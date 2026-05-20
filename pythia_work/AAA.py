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

def lorentz_boost(p4_particle, p4_parent):
    """
    Hace el Boost de Lorentz de una partícula al sistema en reposo de su 'parent' (madre).
    p4_particle y p4_parent deben tener la forma (N, 4) -> [E, px, py, pz]
    """
    # Extraer componentes de la madre (el tau)
    E_M, px_M, py_M, pz_M = p4_parent[:, 0], p4_parent[:, 1], p4_parent[:, 2], p4_parent[:, 3]
    M = np.sqrt(E_M**2 - (px_M**2 + py_M**2 + pz_M**2)) # Masa invariante

    # Calcular velocidad beta = p/E de la madre (con signo opuesto para ir a su reposo)
    bx = -px_M / E_M
    by = -py_M / E_M
    bz = -pz_M / E_M

    b2 = bx**2 + by**2 + bz**2
    gamma = 1.0 / np.sqrt(1.0 - b2)

    # Para evitar divisiones por cero si beta es cero
    gamma2 = np.where(b2 > 0, (gamma - 1.0) / b2, 0.0)

    # Extraer componentes de la hija (el pión)
    E, px, py, pz = p4_particle[:, 0], p4_particle[:, 1], p4_particle[:, 2], p4_particle[:, 3]

    # Producto punto escalar de beta y el trimomento de la hija
    b_dot_p = bx*px + by*py + bz*pz

    # Transformaciones de Lorentz vectoriales generales
    E_prime  = gamma * (E + b_dot_p)
    px_prime = px + gamma2 * b_dot_p * bx + gamma * E * bx
    py_prime = py + gamma2 * b_dot_p * by + gamma * E * by
    pz_prime = pz + gamma2 * b_dot_p * bz + gamma * E * bz

    return np.column_stack((E_prime, px_prime, py_prime, pz_prime))


piplus_4p_rest = lorentz_boost(piplus_4p, tau_4p)

# 2. Extraemos el NUEVO trimomento en reposo
piplus_3p_rest = piplus_4p_rest[:, 1:4]

# 3. Calculamos los cosenos directores usando el momento en REPOSO
mag_pi_rest = np.linalg.norm(piplus_3p_rest, axis=1, keepdims=True)

p_piplus_n = np.sum(piplus_3p_rest * n_hat, axis=1, keepdims=True) / mag_pi_rest
p_piplus_r = np.sum(piplus_3p_rest * r_hat, axis=1, keepdims=True) / mag_pi_rest
p_piplus_k = np.sum(piplus_3p_rest * k_hat, axis=1, keepdims=True) / mag_pi_rest

# 4. Graficar (usamos density=True para ver frecuencias relativas normalizadas)
plt.hist(p_piplus_n.flatten(), bins=40, density=True, alpha=0.5, label=r'$\cos\theta_n$')
plt.hist(p_piplus_r.flatten(), bins=40, density=True, alpha=0.5, label=r'$\cos\theta_r$')
plt.hist(p_piplus_k.flatten(), bins=40, density=True, alpha=0.5, label=r'$\cos\theta_k$')

plt.xlim(-1, 1)
plt.ylim(0, 1) # Como es uniforme de -1 a 1, la densidad debería ser plana cerca de 0.5
plt.legend()
plt.show()





print(len(Z_boson))
print(len(tau))
print(len(piplus))
print(p_piplus_n)

