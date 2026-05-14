import numpy as np
import csv

data = []
with open('tomography_data.csv', mode='r', newline='') as file:
    reader = csv.reader(file)
    next(reader)
    for row in reader:
        numeric_row = list(map(float, row))
        data.append(numeric_row)

data_array = np.array(data)

higgs = []
taus = []
antitaus = []

for i in range(len(data_array)):
    if (int(data_array[i][1]) == 25):
        higgs.append(data_array[i])
    elif (int(data_array[i][1]) == 15):
        taus.append(data_array[i])
    elif (int(data_array[i][1]) == -15):
        antitaus.append(data_array[i])

higgs_arr = np.array(higgs)
taus_arr = np.array(taus)
antitaus_arr = np.array(antitaus)



print(len(higgs_arr))
print(len(taus_arr))
print(len(antitaus_arr))
