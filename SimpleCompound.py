from Compound import Compound
from dataCollection import periodicTable, element_color_map
from Bond import Bond
from math import sqrt, degrees, acos

import numpy as np

import matplotlib.pyplot as plt
from itertools import combinations

def get_perpendicular_vector(v):
    # Return a normalized vector perpendicular to v
    # Choose an arbitrary vector not parallel to v
    if np.allclose(v, [0, 0, 1]):
        other = np.array([1, 0, 0])
    else:
        other = np.array([0, 0, 1])
    perp = np.cross(v, other)
    return perp / np.linalg.norm(perp)

def angle_between(v1, v2):
    # Make sure vectors are float arrays
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)

    # Compute angle in degrees
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    # Clamp to avoid NaN from rounding errors
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return degrees(acos(cos_theta))

# one central atom, only nonmetals
class SimpleCompound(Compound):
    def __init__(self, equation: str):
        super().__init__(equation)
                    
    def totalValenceElectrons(self):
        total = 0
        for element in self.elements: total += periodicTable[element][3] * self.elements[element]
        
        return int(total)
    
    def generateLewisStructure(self):
        # 1. Build atom list (flattened from self.elements)
        atoms = []
        for element, count in self.elements.items():
            atoms.extend([element] * count)

        # 2. Get index of the central atom (first matching element)
        central_symbol = self.getCentralAtom()
        central_idx = atoms.index(central_symbol)

        # 3. Calculate total valence electrons
        totalElectrons = self.totalValenceElectrons()

        # 4. Create single bonds between central atom and all other atoms
        bonds = []
        electronsUsed = 0
        electronCount = [0] * len(atoms)  # Track electrons per atom

        for i, element in enumerate(atoms):
            if i == central_idx:
                continue
            bond = Bond(central_idx, i, atoms, order=1)  
            bonds.append(bond)
            
            electronsUsed += 2  # each bond uses 2 electrons
            electronCount[central_idx] += 2 
            electronCount[i] += 2

        # 6. Distribute remaining electrons as lone pairs
        remainingElectrons = totalElectrons - electronsUsed
        lonePairs = [0] * len(atoms)

        # First pass: satisfy peripheral atoms (non-central atoms)
        for i, element in enumerate(atoms):
            if i == central_idx:
                continue 
            
            maxElectrons = 2 if element == "H" else 8
            needed = max(0, maxElectrons - electronCount[i])
            pairs = min(remainingElectrons, needed) // 2
            lonePairs[i] = pairs
            electronCount[i] += pairs * 2
            remainingElectrons -= pairs * 2

        # Second pass: place remaining electrons on central atom
        if remainingElectrons > 0:
            pairs = remainingElectrons // 2
            lonePairs[central_idx] = pairs
            electronCount[central_idx] += pairs * 2
            remainingElectrons -= pairs * 2

        # 7. Add double/triple bonds if central atom lacks octet
        attempts = 0
        maxAttempts = 10  # Prevent infinite loops

        while electronCount[central_idx] < 8 and attempts < maxAttempts:
            madeChange = False
            attempts += 1
            
            for bond in bonds:
                if not bond.involves(central_idx): continue

                terminal_idx = bond.other(central_idx)

                if atoms[terminal_idx] == "H":
                    continue  # Hydrogen can't form multiple bonds

                if lonePairs[terminal_idx] > 0 and bond.order < 3:
                    # Upgrade bond and reduce lone pair
                    bond.order += 1

                    lonePairs[terminal_idx] -= 1
                    electronCount[central_idx] += 2

                    madeChange = True
                    break

            if not madeChange:
                break

        # 8. Final output
        return {
            "atoms": atoms,              # list of element symbols
            "central": central_idx,      # index of central atom
            "bonds": bonds,              # list of (a_idx, b_idx, bondOrder)
            "lonePairs": lonePairs,      # list of lone pair counts per atom
            "electronCount": electronCount,
            "remainingElectrons": remainingElectrons
        }
    
    def generateVSEPRStructure(self):
        LewisStructure = self.generateLewisStructure()

        central_idx = LewisStructure["central"]

        lone_pairs = LewisStructure["lonePairs"][central_idx]
        steric_number = len(LewisStructure["bonds"]) + lone_pairs

        ideal_positions = self._generateIdealPositions(steric_number)

        ideal_positions = [np.array(pos) for pos in ideal_positions]

        actual_positions = self._adjustForLonePairs(ideal_positions, lone_pairs, bonds=LewisStructure["bonds"])
        # updates each Bond object with the new angleFromCentral vector

        return LewisStructure

    def _generateIdealPositions(self, steric_number):
        if steric_number == 2 or steric_number == 1: 
            # Linear
            return [
                (1.0, 0.0, 0.0),
                (-1.0, 0.0, 0.0)
            ]

        elif steric_number == 3:
            # Trigonal planar
            return [
                (1.0, 0.0, 0.0),
                (-0.5, sqrt(3)/2, 0.0),
                (-0.5, -sqrt(3)/2, 0.0)
            ]
        
        elif steric_number == 4:
            # Tetrahedral
            return [
                (1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3)),
                (-1.0 / sqrt(3), -1.0 / sqrt(3), 1.0 / sqrt(3)),
                (-1.0 / sqrt(3), 1.0 / sqrt(3), -1.0 / sqrt(3)),
                (1.0 / sqrt(3), -1.0 / sqrt(3), -1.0 / sqrt(3))
            ]
        
        elif steric_number == 5:
            # Trigonal bipyramidal
            return [
                (0.0, 0.0, 1.0),      # Axial
                (0.0, 0.0, -1.0),     # Axial
                (1.0, 0.0, 0.0),      # Equatorial
                (-0.5, sqrt(3)/2, 0.0),  # Equatorial
                (-0.5, -sqrt(3)/2, 0.0)  # Equatorial
            ]        
        
        elif steric_number == 6:
            # Octahedral
            return [
                (1.0, 0.0, 0.0),      # +X
                (-1.0, 0.0, 0.0),     # -X
                (0.0, 1.0, 0.0),      # +Y
                (0.0, -1.0, 0.0),     # -Y
                (0.0, 0.0, 1.0),      # +Z
                (0.0, 0.0, -1.0)      # -Z
            ]
        
        else:
            raise ValueError(f"Unsupported steric number: {steric_number}, type: {type(steric_number)}")

    def _adjustForLonePairs(self, ideal_positions, lone_pairs, bonds):
            # returns adjusted positions for lone pairs and bonds as vectors
            lp_vectors = ideal_positions[:lone_pairs]
            bond_vectors = ideal_positions[lone_pairs:]

            adjusted_bonds = []

            bond_order_weight = {
                1: 1.0,
                2: 0.8,
                3: 0.6
            }

            k = 0.1

            for bond_vec, bond in zip(bond_vectors, bonds):
                weight = bond_order_weight.get(int(bond.order), 1)  # Default weight if not found
                repulsion = np.zeros(3)

                for lp in lp_vectors:
                    # Lone pair repels bond vector — stronger when aligned
                    direction = bond_vec - lp
                    if np.linalg.norm(direction) == 0:
                        continue  # Avoid divide-by-zero
                    direction /= np.linalg.norm(direction)

                    repulsion += direction * weight * k  # 0.2 is base strength

                new_vec = bond_vec + repulsion
                new_vec /= np.linalg.norm(new_vec)  # Normalize
                
                bond.angleFromCentral = new_vec

                adjusted_bonds.append(new_vec)

            return [vec for vec in lp_vectors] + adjusted_bonds

    def getBondAngles(self):
        structure = self.generateVSEPRStructure()
        central_idx = structure["central"]
        atoms = structure["atoms"]
        bonds = structure["bonds"]

        angles = []

        for bond1, bond2 in combinations(bonds, 2):
            v1 = np.array(bond1.angleFromCentral)
            if bond1.b == central_idx:
                v1 = -v1
            v1 = v1 * bond1.bondLength()

            v2 = np.array(bond2.angleFromCentral)
            if bond2.b == central_idx:
                v2 = -v2
            v2 = v2 * bond2.bondLength()

            angle = angle_between(v1, v2)

            a1 = atoms[bond1.a if bond1.a != central_idx else bond1.b]
            a2 = atoms[bond2.a if bond2.a != central_idx else bond2.b]
            central = atoms[central_idx]

            angles.append((a1, central, a2, angle))

        return angles

    def displayMolecule(self):
        LewisStructure = self.generateVSEPRStructure()

        bonds = LewisStructure["bonds"]
        atoms = LewisStructure["atoms"]
        central_idx = LewisStructure["central"]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Central atom at origin
        central_pos = np.array([0, 0, 0])
        atom_positions = [central_pos]
        labels = [atoms[central_idx]]

        for bond in bonds:
            # Determine direction and non-central atom
            if bond.a == central_idx:
                direction = np.array(bond.angleFromCentral)
                target_atom = atoms[bond.b]
            elif bond.b == central_idx:
                direction = -np.array(bond.angleFromCentral)
                target_atom = atoms[bond.a]
            else:
                continue  # skip if not bonded to central atom

            # Scale direction vector by bond length
            end = direction * bond.bondLength()
            atom_positions.append(end)
            labels.append(target_atom)

            # Draw single/double/triple bonds as parallel lines
            num_lines = bond.order
            if num_lines == 1:
                ax.plot(
                    [central_pos[0], end[0]],
                    [central_pos[1], end[1]],
                    [central_pos[2], end[2]],
                    color='gray',
                    linewidth=1.5
                )
            else:
                offset_dir = get_perpendicular_vector(direction)
                spacing = 0.05  # Ångstroms of separation

                offsets = np.linspace(-spacing, spacing, num_lines)
                for offset in offsets:
                    shift = offset * offset_dir
                    ax.plot(
                        [central_pos[0] + shift[0], end[0] + shift[0]],
                        [central_pos[1] + shift[1], end[1] + shift[1]],
                        [central_pos[2] + shift[2], end[2] + shift[2]],
                        color='gray',
                        linewidth=1.0
                    )

        # Draw atoms as scatter points
        atom_positions[0], atom_positions[central_idx] = atom_positions[central_idx], atom_positions[0]  # Move central atom to first position
        labels[0], labels[central_idx] = labels[central_idx], labels[0]  # Swap labels

        xs, ys, zs = zip(*atom_positions)

        # Get colors and marker sizes based on covalent radius (index 1)
        colors = [element_color_map.get(elem, element_color_map["default"]) for elem in atoms]
        sizes = [
            periodicTable.get(elem, [None, 0.5])[1] * 500  # scale angstroms to point size; fallback = 0.5 Å
            for elem in atoms
        ]

        ax.scatter(xs, ys, zs, s=sizes, c=colors, alpha = 1)

        # Add text labels
        for (x, y, z), label in zip(atom_positions, labels):
            ax.text(x, y, z, label, fontsize=10, ha='center', color = 'black')

        # Formatting
        ax.set_box_aspect([1, 1, 1])
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ax.set_zlim(-2, 2)
        ax.axis('off')

        plt.tight_layout()
        plt.show()


cmpd = SimpleCompound("PF5")

LewisStructure = cmpd.generateLewisStructure()

print(LewisStructure)
print()
print(cmpd.getBondAngles())

cmpd.displayMolecule()

