# from Compound import compound

# compounds = []


# with open("testCompounds.txt", "r") as file:
#     for line in file:
#         if "MF:" in line:
#             start = line.find("MF:") + len("MF:")
#             mf_value = line[start:].strip()
#             formula = mf_value
        
#             comp = compound(formula)

#             if comp.totalElements() > 2 and comp.totalElements() < 5 and comp.numUniqueElements() < 3:
#                 compounds.append(formula)

# print(compounds[:100])
# print(len(compounds))

# # save compounds to a comma separated file
# with open("filtered_compounds.txt", "w") as file:
#     for comp in compounds:
#         file.write(comp + "\n")

import csv

nonmetals = []
metals = []
metaloids = []

completeTable = []

with open("bonddata.csv", "r") as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header row
    periodicTable = list(reader)

    for row in periodicTable:
        # data = [row[2], row[0], row[16], row[17], row[27]]

        # if row[12] == "yes": metals.append(data)
        # if row[13] == "yes": nonmetals.append(data)
        # if row[14] == "yes": metaloids.append(data)

        # completeTable.append(data)
        num = row[6]

        if num == "": num = 0
        completeTable.append(float(num) / 100)

import pandas as pd

df = pd.read_csv("TrimmedTable.csv")

df = df.drop('Atomic Radius.1', axis=1)

df.to_csv("TrimmedTable.csv", index=False)

# with open("TrimmedTable.csv", "w", newline='') as file:
#     writer = csv.writer(file)
#     writer.writerow(["Symbol", "Atomic Number", "Atomic Radius", "Electronegativity", "Valence"])
#     writer.writerows(completeTable)

# # save data to their respective CSV files
# with open("Nonmetals.csv", "w", newline='') as file:
#     writer = csv.writer(file)
#     writer.writerow(["Symbol", "Atomic Number", "Atomic Radius", "Electronegativity", "Valence"])
#     writer.writerows(nonmetals)

# with open("Metals.csv", "w", newline='') as file:
#     writer = csv.writer(file)
#     writer.writerow(["Symbol", "Atomic Number", "Atomic Radius", "Electronegativity", "Valence"])
#     writer.writerows(metals)

# with open("Metaloids.csv", "w", newline='') as file:
#     writer = csv.writer(file)
#     writer.writerow(["Symbol", "Atomic Number", "Atomic Radius", "Electronegativity", "Valence"])
#     writer.writerows(metaloids)
