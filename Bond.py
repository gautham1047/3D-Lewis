from dataCollection import periodicTable

class Bond:
    def __init__(self, atom1_idx, atom2_idx, atoms, order = 1):
        self.a = atom1_idx  # index of first atom
        self.b = atom2_idx  # index of second atom
        self.order = order  # bond order: 1 = single, 2 = double, 3 = triple

        self.atoms = atoms

        self.angleFromCentral = None  # direction from central atom as a NumPy vector

    def involves(self, idx):
        return self.a == idx or self.b == idx

    def other(self, idx):
        if self.a == idx:
            return self.b
        elif self.b == idx:
            return self.a
        raise ValueError("Atom not in bond")
    
    def bondLength(self):
        firstAtom = periodicTable[self.atoms[self.a]]
        secondAtom = periodicTable[self.atoms[self.b]]

        return (firstAtom[1] + secondAtom[1]) * (1.1 - .1 * self.order)

    def __repr__(self):
        return f"Bond({self.a}, {self.b}, order={self.order})"

bond = Bond(0, 1, ["C", "C"], order = 1)

print(bond.bondLength())

bond = Bond(0, 1, ["C", "C"], order = 2)

print(bond.bondLength())

bond = Bond(0, 1, ["C", "C"], order = 3)

print(bond.bondLength())