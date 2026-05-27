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
Z_boson, tau, antitau, piplus, neutau = data_array[ids == 23], data_array[ids == 15], data_array[ids == -15], data_array[ids == 211], data_array[ids == 16]

# Construimos los cuatrimomentos de las partículas objetivo (Energía, 3momento)
Z_4p, tau_4p, antitau_4p, piplus_4p, neutau_4p = Z_boson[:, 3:7].astype(float), tau[:, 3:7].astype(float), antitau[:, 3:7].astype(float), piplus[:, 3:7].astype(float), neutau[:, 3:7].astype(float)

print(len(tau))
print(len(piplus))
