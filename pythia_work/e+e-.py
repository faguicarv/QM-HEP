import numpy as np
import csv

data = []
with open('e+e-_data.csv', mode='r', newline='') as file:
    reader = csv.reader(file)
    next(reader)
    for row in reader:
        numeric_row = list(map(float, row))
        data.append(numeric_row)

data_array = np.array(data)
Z_boson, tau, antitau, piplus, neutau = [], [], [], [], []

for i in range(len(data_array)):
    if (int(data_array[i][1]) == 23):
        Z_boson.append(data_array[i])
    elif (int(data_array[i][1]) == 15):
        tau.append(data_array[i])
    elif (int(data_array[i][1]) == -15):
        antitau.append(data_array[i])
    elif (int(data_array[i][1]) == 211):
        piplus.append(data_array[i])

    elif (int(data_array[i][1]) == 16):
        neutau.append(data_array[i])


Z_boson_arr, tau_arr, antitau_arr, piplus_arr, neutau_arr = np.array(Z_boson), np.array(tau), np.array(antitau), np.array(piplus), np.array(neutau)

#Ya tenemos las listas de cada partícula. Ahora comenzaremos a obtener solo los momentum y ponerlos en un array para realizar la reconstrucción.





print(len(Z_boson_arr))
print(len(tau_arr))
print(len(antitau_arr))
print(len(piplus_arr))
print(len(neutau_arr))
