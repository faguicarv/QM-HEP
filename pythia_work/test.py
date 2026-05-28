import numpy as np
import csv
import matplotlib.pyplot as plt

data = []
with open('e+e-_data.csv', mode='r', newline='') as file:
    reader = csv.reader(file)
    next(reader)
    for row in reader:
        numeric_row = list(map(float, row))
        data.append(numeric_row)

# Mucho más fácil trabajar con arrays
data_array = np.array(data)

#creamos un array de solamente las id de las partículas
ids = data_array[:, 1].astype(int)

# Filtramos directamente por id de partícula y creamos un array con toda la información correspondiente
Z_boson, tau, antitau, piplus, piminus, neutau, antineutau = data_array[ids == 23], data_array[ids == 15], data_array[ids == -15], data_array[ids == 211], data_array[ids == -211], data_array[ids == 16], data_array[ids == -16]

print(len(Z_boson))
print(len(tau))
print(len(antitau))
print(len(piplus))
print(len(piminus))
print(len(neutau))
print(len(antineutau))
