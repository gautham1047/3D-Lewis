# from ase.io import read, write

# cif_files = ["1008071.cif"]

# for cif_file in cif_files:
#     atoms = read(cif_file)

#     xyz_file = cif_file.replace(".cif", ".xyz")

#     write(xyz_file, atoms, format="xyz")


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Atom color mapping (extend as needed)
atom_colors = {
    'H': 'white',
    'C': 'black',
    'N': 'blue',
    'O': 'red',
    'F': 'green',
    'Cl': 'green',
    'Br': 'brown',
    'I': 'purple',
    'P': 'orange',
    'S': 'yellow',
}

def read_xyz(filename):
    atoms = []
    with open(filename, 'r') as file:
        lines = file.readlines()[2:]  # Skip first 2 lines
        for line in lines:
            parts = line.split()
            if len(parts) >= 4:
                symbol = parts[0]
                x, y, z = map(float, parts[1:4])
                atoms.append((symbol, x, y, z))
    return atoms

def plot_xyz(atoms):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for atom in atoms:
        symbol, x, y, z = atom
        color = atom_colors.get(symbol, 'gray')
        ax.scatter(x, y, z, color=color, s=200, edgecolors='k', label=symbol)

    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title('Molecule from .xyz file')

    # Avoid duplicate labels in legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())

    plt.show()

# === Usage ===
atoms = read_xyz('1008071.xyz')  # replace with your actual file
plot_xyz(atoms)
