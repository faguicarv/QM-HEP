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

# Veamos cuántos eventos vamos a tener
N_events = len(Z_boson)
# Para definir n, r, k, debemos definir el 3-momento del positron (eso dice el paper)
eplus_4p = np.zeros((N_events, 4))
E_beam = 91.1876 / 2 # e+ y e- tienen mitad de energía del CM
eplus_4p[:, 0], eplus_4p[:, 1], eplus_4p[:, 2], eplus_4p[:, 3] = E_beam, 0.0, 0.0, E_beam # Primera columna es energía, Segunda columna es px, Tercera columna es py, Cuarta columna es pz (el haz va en dirección z)

# El paper dice que p es la dirección del positron incidente, como es un vector vamos a obtener el 3momento y así obtenemos p (Que sólo debería contener una componente en z)
eplus_3p = eplus_4p[:, 1:4].astype(float)
# Ahora obtenemos el vector p para las ecuaciones (2.2)
p_hat = eplus_3p / np.linalg.norm(eplus_3p, axis=1, keepdims=True) # EN REALIDAD BASTABA CON PONER TANTOS 1 EN Z EN EL 3MOMENTO DEL e^+, PERO ASÍ SE ENTIENDE EL SISTEMA DE REFERENCIA DE MEJOR FORMA PARA CONSTRUIR EL RESTO DE VECTORES

# Ahora buscamos k, que es la dirección del momentum del tau^+
tau_3p = tau_4p[:, 1:4].astype(float) # 3-momento del tau
k_hat = tau_3p / np.linalg.norm(tau_3p, axis=1, keepdims=True) # Normalizamos el vector k

# Como ya tenemos p y k, podemos obtener cos theta
cos_theta = np.sum(p_hat * k_hat, axis=1, keepdims=True)
sin_theta = np.sqrt(1 - cos_theta**2)

# Calculemos n y r
n_hat = np.cross(p_hat, k_hat, axis=1) / sin_theta
r_hat = (p_hat - (k_hat * cos_theta)) / sin_theta

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# HASTA ACÁ SE DEFINEN TODAS LAS COMPONENTES DE LA BASE EN EL FERMION-PAIR CENTER OF MASS FRAME Y ENCONTRAMOS TODAS LAS COMPONENTES DE LA BASE (n, r, k) DE LAS ECUACIONES (2.2)
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------


# Calculemos las proyecciones de la dirección del momentum del pion en la base n,r,k
piplus_3p = piplus_4p[:, 1:4].astype(float)

p_piplus_n = np.sum(piplus_3p * n_hat, axis=1, keepdims=True) / np.linalg.norm(piplus_3p, axis=1, keepdims=True)
p_piplus_r = np.sum(piplus_3p * r_hat, axis=1, keepdims=True) / np.linalg.norm(piplus_3p, axis=1, keepdims=True)
p_piplus_k = np.sum(piplus_3p * k_hat, axis=1, keepdims=True) / np.linalg.norm(piplus_3p, axis=1, keepdims=True)

# Graficar
plt.figure(figsize=(8, 6))
plt.hist(p_piplus_n, bins=100, histtype='step', linewidth=2, label=r'$\cos\theta_n^+$')
plt.hist(p_piplus_r, bins=100, histtype='step', linewidth=2, label=r'$\cos\theta_r^+$')
plt.hist(p_piplus_k, bins=100, histtype='step', linewidth=2, label=r'$\cos\theta_k^+$')
plt.xlim(-1.1, 1.1)
plt.xlabel(r"$\cos\theta_n$")
plt.ylabel("Eventos")
plt.title("Proyección $\pi^+$-momentum en la base ($\hat{n}, \hat{r}, \hat{k}$)")
plt.legend()
plt.show()

# Como todavía los gráficos son planos, vamos a trabajar en realizar un boost del 3-momento del pion+ hacia el marco de referencia de la desintegración del tau

# 3-momento y Energía del tau
tau_energy = tau_4p[:, 0, np.newaxis].astype(float)

# Factores de Lorentz beta y gamma
beta, gamma = tau_3p / tau_energy, tau_energy / tau[:, 2, np.newaxis].astype(float)

print(beta)
print(len(beta))
print(gamma)
print(len(gamma))

#def boost_to_rest_frame(p4, beta, gamma):
    #"""
    #Aplica un boost de Lorentz a un arreglo de 4-momentos para llevarlos
    #al sistema de reposo definido por la velocidad beta.
    #"""
    #E = p4[:, 0:1]
    #p3 = p4[:, 1:4]
    #gamma = gamma.reshape(-1, 1)

    ## ¡AQUÍ ESTABA EL ERROR! No hay que invertir beta.
    ## Las ecuaciones ya asumen la transformación de CM a Rest Frame.
    #v = beta

    #v_dot_p3 = np.sum(v * p3, axis=1, keepdims=True)
    #v2 = np.sum(v**2, axis=1, keepdims=True)

    #factor = np.zeros_like(v2)
    #mask = v2[:, 0] > 0
    #factor[mask] = (gamma[mask] - 1.0) / v2[mask]

    ## Ecuaciones correctas del boost de Lorentz
    #E_prime = gamma * (E - v_dot_p3)
    #p3_prime = p3 + factor * v_dot_p3 * v - gamma * E * v

    #p4_prime = np.hstack((E_prime, p3_prime))

    #return p4_prime

## Realizar boost a todos los productos (ZRF = Z-boson rest frame )
#Z_ZRF = boost_to_rest_frame(Z_4p, beta, gamma)
#tau_ZRF = boost_to_rest_frame(tau_4p, beta, gamma)
#antitau_ZRF = boost_to_rest_frame(antitau_4p, beta, gamma)
#piplus_ZRF = boost_to_rest_frame(piplus_4p, beta, gamma)
#neutau_ZRF = boost_to_rest_frame(neutau_4p, beta, gamma)
#eplus_ZRF = boost_to_rest_frame(eplus_4p, beta, gamma)

## Todos los productos de desintegración boosteados, vamos a buscar las coordenadas también boosteadas (n, r, k)
#p_hat_ZRF = eplus_ZRF[:, 1:4] / np.linalg.norm(eplus_ZRF[:, 1:4], axis=1, keepdims=True)

#k_hat_ZRF = tau_ZRF[:, 1:4] / np.linalg.norm(tau_ZRF[:, 1:4], axis=1, keepdims=True)

#cos_theta_ZRF = np.sum(p_hat_ZRF * k_hat_ZRF, axis=1, keepdims=True)
#sin_theta_ZRF = np.sqrt(1 - cos_theta_ZRF**2)

#n_hat_ZRF = np.cross(p_hat_ZRF, k_hat_ZRF, axis=1) / sin_theta_ZRF
#r_hat_ZRF = (p_hat_ZRF - (k_hat_ZRF * cos_theta_ZRF)) / sin_theta_ZRF

## Proyecciones del pi+ en esta nueva base boosteada
#p_piplus_n_ZRF = np.sum(piplus_ZRF[:,1:4] * n_hat_ZRF, axis=1, keepdims=True) / np.linalg.norm(piplus_ZRF[:,1:4], axis=1, keepdims=True)
#p_piplus_r_ZRF = np.sum(piplus_ZRF[:,1:4] * r_hat_ZRF, axis=1, keepdims=True) / np.linalg.norm(piplus_ZRF[:,1:4], axis=1, keepdims=True)
#p_piplus_k_ZRF = np.sum(piplus_ZRF[:,1:4] * k_hat_ZRF, axis=1, keepdims=True) / np.linalg.norm(piplus_ZRF[:,1:4], axis=1, keepdims=True)


##print(p_hat_ZRF)
##print(k_hat_ZRF)
##print(p_hat_ZRF)

##plt.hist(p_piplus_n_ZRF, bins=100)
##plt.xlim(-1, 1)
##plt.show()
## plt.hist(p_piplus_r)
## plt.hist(p_piplus_k)



#p_hat_CM = eplus_4p[:, 1:4] / np.linalg.norm(eplus_4p[:, 1:4], axis=1, keepdims=True)

## Momento del tau original
#k_hat_CM = tau_4p[:, 1:4] / np.linalg.norm(tau_4p[:, 1:4], axis=1, keepdims=True)

#cos_theta_CM = np.sum(p_hat_CM * k_hat_CM, axis=1, keepdims=True)
#sin_theta_CM = np.sqrt(1 - cos_theta_CM**2)

#n_hat_CM = np.cross(p_hat_CM, k_hat_CM, axis=1) / sin_theta_CM
#r_hat_CM = (p_hat_CM - (k_hat_CM * cos_theta_CM)) / sin_theta_CM


## -------------------------------------------------------------------
## 2. LLEVAR LOS PRODUCTOS DE DECAIMIENTO AL TAU REST FRAME (TRF)
## Usando el beta y gamma que calculaste a partir del tau_4p
## -------------------------------------------------------------------

## Solo necesitas hacerle el boost a los productos, no al tau ni al haz inicial
#piplus_TRF = boost_to_rest_frame(piplus_4p, beta, gamma)
## neutau_TRF = boost_to_rest_frame(neutau_4p, beta, gamma) # etc...

## Obtenemos la dirección del pion en el rest frame del tau
#p_piplus_hat_TRF = piplus_TRF[:, 1:4] / np.linalg.norm(piplus_TRF[:, 1:4], axis=1, keepdims=True)


## -------------------------------------------------------------------
## 3. PROYECTAR EL PION (TRF) SOBRE LA BASE (CM)
## -------------------------------------------------------------------

## Cosenos directores (proyecciones)
#cos_theta_pi_n = np.sum(p_piplus_hat_TRF * n_hat_CM, axis=1, keepdims=True)
#cos_theta_pi_r = np.sum(p_piplus_hat_TRF * r_hat_CM, axis=1, keepdims=True)
#cos_theta_pi_k = np.sum(p_piplus_hat_TRF * k_hat_CM, axis=1, keepdims=True)

