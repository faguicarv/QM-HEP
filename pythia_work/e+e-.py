import numpy as np
import csv

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

# Construimos los cuatrimomentos de las partículas objetivo
Z_4p, tau_4p, antitau_4p, piplus_4p, neutau_4p = Z_boson[:, 2:6].astype(float), tau[:, 2:6].astype(float), antitau[:, 2:6].astype(float), piplus[:, 2:6].astype(float), neutau[:, 2:6].astype(float)

# Veamos cuántos eventos vamos a tener
N_events = len(Z_boson)
# Para definir n, r, k, debemos definir el 3-momento del positron (eso dice el paper)
eplus_4p = np.zeros((N_events, 4))
E_beam = 91.1876 / 2 # e+ y e- tienen mitad de energía del CM
eplus_4p[:, 0] = E_beam # Primera columna es energía
eplus_4p[:, 1] = 0.0 # Segunda columna es px
eplus_4p[:, 2] = 0.0 # Tercera columna es py
eplus_4p[:, 3] = E_beam # Cuarta columna es pz (el haz va en dirección z)

# El paper dice que p es la dirección del positron incidente, como es un vector vamos a obtener el 3momento y así obtenemos p
eplus_3p = eplus_4p[:, 1:4].astype(float)
# Ahora obtenemos el vector p para las ecuaciones (2.2)
p_hat = eplus_3p / np.linalg.norm(eplus_3p, axis=1, keepdims=True) # EN REALIDAD BASTABA CON PONER TANTOS 1 EN Z EN EL 3MOMENTO DEL e^+, PERO ASÍ SE ENTIENDE EL SISTEMA DE REFERENCIA DE MEJOR FORMA PARA CONSTRUIR EL RESTO DE VECTORES

# Ahora buscamos k, que es la dirección del momentum del tau
tau_3p = tau_4p[:, 1:4].astype(float) # 3-momento del tau
k_hat = tau_3p / np.linalg.norm(tau_3p, axis=1, keepdims=True) # Normalizamos el vector k

# Como ya tenemos p y k, podemos obtener cos theta
cos_theta = np.sum(p_hat * k_hat, axis=1, keepdims=True)
sin_theta = np.sqrt(1 - cos_theta**2)

# Calculemos n y r
n_hat = np.cross(p_hat, k_hat, axis=1) / sin_theta
r_hat = (p_hat - (k_hat * cos_theta)) / sin_theta

# Calculemos las proyecciones de la dirección del momentum del pion en la base n,r,k
piplus_3p = piplus_4p[:, 1:4].astype(float)
# p_piplus_n =






print(len(Z_boson))















