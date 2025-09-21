import csv

# All data is formatted as a dictionary with the element symbol as the key
# Values are lists containing the atomic number, atomic radius, electronegativity, and valence electrons

with open("Nonmetals.csv", "r") as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header row

    nonmetalsData = dict((row[0], [float(r) if r != "" else 0 for r in row[1:]]) for row in reader)

with open("Metals.csv", "r") as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header row

    metalsData = dict((row[0], [float(r) if r != "" else 0 for r in row[1:]]) for row in reader)

with open("Metaloids.csv", "r") as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header row

    metaloidsData = dict((row[0], [float(r) if r != "" else 0 for r in row[1:]]) for row in reader)

with open("TrimmedTable.csv", "r") as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header row

    periodicTable = dict((row[0], [float(r) if r != "" else 0 for r in row[1:]]) for row in reader)

element_color_map = {
    # Nonmetals
    "H": "gray",
    "C": "black",
    "N": "blue",
    "O": "red",
    "F": "green",
    "P": "orange",
    "S": "yellow",
    "Cl": "lime",
    "Br": "maroon",
    "I": "purple",

    # Noble gases
    "He": "cyan",
    "Ne": "cyan",
    "Ar": "cyan",
    "Kr": "cyan",
    "Xe": "cyan",
    "Rn": "cyan",

    # Other common elements
    "B": "salmon",
    "Si": "brown",
    "Se": "gold",
    "Li": "violet",
    "Na": "deepskyblue",
    "K": "orchid",
    "Rb": "plum",
    "Cs": "indigo",
    "Fr": "slateblue",

    "Be": "turquoise",
    "Mg": "lightgreen",
    "Ca": "lightblue",
    "Sr": "steelblue",
    "Ba": "cornflowerblue",
    "Ra": "navy",

    "Ti": "gray",
    "Fe": "darkorange",
    "Zn": "lightgray",
    "Cu": "saddlebrown",
    "Al": "lightcoral",

    # Default (used if element not found)
    "default": "midnightblue"
}