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
Z_boson, taus, antitaus = [], [], []

for i in range(len(data_array)):
    if (int(data_array[i][1]) == 23):
        Z_boson.append(data_array[i])
    elif (int(data_array[i][1]) == 15):
        taus.append(data_array[i])
    elif (int(data_array[i][1]) == -15):
        antitaus.append(data_array[i])

Z_boson_arr, taus_arr, antitaus_arr = np.array(Z_boson), np.array(taus), np.array(antitaus)

print(len(Z_boson_arr))
print(len(taus_arr))
print(len(antitaus_arr))
