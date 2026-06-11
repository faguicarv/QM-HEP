import numpy as np
import csv
import matplotlib.pyplot as plt

files = [f"output_Z_decay_{i}.csv" for i in range(1,401)] # Creamos una lista con todos los nombres de los archivos extraídos de pythia
data = []
for file_name in files:
    with open(file_name, mode='r', newline='') as file:
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

# ----- NOTA: EL TAU ES tau^- (id=15)
# Construimos los cuatrimomentos de las partículas objetivo (Energía, 3momento)
Z_4p, tau_4p, antitau_4p, piplus_4p, piminus_4p, neutau_4p, antineutau_4p = Z_boson[:, 3:7].astype(float), tau[:, 3:7].astype(float), antitau[:, 3:7].astype(float), piplus[:, 3:7].astype(float), piminus[:, 3:7].astype(float), neutau[:, 3:7].astype(float), antineutau[:, 3:7].astype(float)

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
antitau_3p = antitau_4p[:, 1:4].astype(float) # 3-momento del tau^+
k_hat = antitau_3p / np.linalg.norm(antitau_3p, axis=1, keepdims=True) # Normalizamos el vector k

# Como ya tenemos p y k, podemos obtener cos theta
cos_theta = np.sum(p_hat * k_hat, axis=1, keepdims=True)
sin_theta = np.sqrt(1 - cos_theta**2)

# Calculemos n y r
n_hat = np.cross(p_hat, k_hat, axis=1) / sin_theta
k_cos = (k_hat * cos_theta)
r_hat = (p_hat - k_cos) / sin_theta

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# HASTA ACÁ SE DEFINEN TODAS LAS COMPONENTES DE LA BASE EN EL FERMION-PAIR CENTER OF MASS FRAME Y ENCONTRAMOS TODAS LAS COMPONENTES DE LA BASE (n, r, k) DE LAS ECUACIONES (2.2)
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

# Como el pion se encuentra en el CM, vamos a boostearlo: (tRF:= tau rest frame)
# 3-momento y Energía del tau
antitau_energy = antitau_4p[:, 0, np.newaxis].astype(float)
# Factores de Lorentz beta y gamma
beta, gamma = antitau_3p / antitau_energy, antitau_energy / antitau[:, 2, np.newaxis].astype(float)

def boost_to_rest_frame(p4, beta, gamma, aux): # Aplica un boost de Lorentz a un arreglo de 4-momentos para llevarlos al sistema de reposo definido por la velocidad beta.
    E = p4[:, 0:1]
    p3 = p4[:, 1:4]
    gamma = gamma.reshape(-1, 1)

    v_dot_p3 = np.sum(beta * p3, axis=1, keepdims=True)
    v2 = np.sum(beta**2, axis=1, keepdims=True)

    factor, mask = np.zeros_like(v2), v2[:, 0] > 0
    factor[mask] = (gamma[mask] - 1.0) / v2[mask]

    # Ecuaciones boost de Lorentz
    E_prime = gamma * (E - v_dot_p3)
    if (aux == +1):
        p3_prime = p3 + factor * v_dot_p3 * beta + gamma * E * beta
    elif (aux == -1):
        p3_prime = p3 + factor * v_dot_p3 * beta - gamma * E * beta

    p4_prime = np.hstack((E_prime, p3_prime))

    return p4_prime

# La función anterior nos ayuda a realizar las operaciones del boost.
# Realicemos un boost al 4-momento del pion^+ hacia el rest-frame del tau
piplus_4p_tRF = boost_to_rest_frame(piplus_4p, beta, gamma, -1) #4-momento boosteado
piplus_3p_tRF = piplus_4p_tRF[:, 1:4] #3-momento boosteado
piplus_tRF = piplus_3p_tRF / np.linalg.norm(piplus_3p_tRF, axis=1, keepdims=True) # 3-momento normalizado

# Obtenemos las proyecciones de la dirección del 3-momento del pión en la base n,r,k
cos_n_plus = np.sum(piplus_tRF * n_hat, axis=1, keepdims=True)
cos_r_plus = np.sum(piplus_tRF * r_hat, axis=1, keepdims=True)
cos_k_plus = np.sum(piplus_tRF * k_hat, axis=1, keepdims=True)


# Boost al 4-momento del pion^- en la base n,r,k
piminus_4p_tRF = boost_to_rest_frame(piminus_4p, beta, gamma, 1)
piminus_3p_tRF = piminus_4p_tRF[:, 1:4]
piminus_tRF = piminus_3p_tRF / np.linalg.norm(piminus_3p_tRF, axis=1, keepdims=True)

# Proyecciones del 3-momento del pi^- en la base n,r,k
cos_n_minus = np.sum(piminus_tRF * n_hat, axis=1, keepdims=True)
cos_r_minus = np.sum(piminus_tRF * r_hat, axis=1, keepdims=True)
cos_k_minus = np.sum(piminus_tRF * k_hat, axis=1, keepdims=True)


# Vamos a calcular los productos de estos cosenos para obtener la FIG 2
np_nm, np_rm, np_km = cos_n_plus * cos_n_minus, cos_n_plus * cos_r_minus, cos_n_plus * cos_k_minus
rp_nm, rp_rm, rp_km = cos_r_plus * cos_n_minus, cos_r_plus * cos_r_minus, cos_r_plus * cos_k_minus
kp_nm, kp_rm, kp_km = cos_k_plus * cos_n_minus, cos_k_plus * cos_r_minus, cos_k_plus * cos_k_minus


fig, axs = plt.subplots(2, 3, figsize=(3, 14))
for ax in axs.flat:
    # 1. Fijar los ticks
    ax.set_xticks(np.arange(-1.0, 1.5, 0.5))
    ax.set_yticks(np.arange(0.0, 0.8, 0.1))

    # 2. Ajustar los límites
    ax.set_xlim(-1.0, 1.0)

    # 3. Encender la grilla suave que seguirá a estos ticks
    ax.grid(True, linestyle='--', alpha=0.3)

bins = 150

# Gráfico 1: Proyección en n
axs[0, 0].hist(cos_n_plus, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='darksalmon', label=r'$\cos\theta_n^+$')
axs[0, 0].axhline(y=0.5, color='red', linestyle='--', alpha=0.7, label=r'$y=0.5$')
axs[0, 0].set_title(r"Proyección $\pi^+$-momentum en $\hat{n}$")
axs[0, 0].set_xlabel(r"$\cos\theta_n^+$")
axs[0, 0].legend(loc="lower center")

# Gráfico 2: Proyección en r
axs[0, 1].hist(cos_r_plus, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='lightpink', label=r'$\cos\theta_r^+$')
axs[0, 1].axhline(y=0.5, color='red', linestyle='--', alpha=0.7, label=r'$y=0.5$')
axs[0, 1].set_title(r"Proyección $\pi^+$-momentum en $\hat{r}$")
axs[0, 1].set_xlabel(r"$\cos\theta_r^+$")
axs[0, 1].legend(loc="lower center")

# Gráfico 3: Proyección en k
axs[0, 2].hist(cos_k_plus, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='thistle', label=r'$\cos\theta_k^+$')
axs[0, 2].axhline(y=0.5, color='red', linestyle='--', alpha=0.7, label=r'$y=0.5$')
axs[0, 2].set_title(r"Proyección $\pi^+$-momentum en $\hat{k}$")
axs[0, 2].set_xlabel(r"$\cos\theta_k^+$")
axs[0, 2].legend(loc="lower center")

# Gráfico 4: Proyección cos_n^-
axs[1, 0].hist(cos_n_minus, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='darksalmon', label=r'$\cos\theta_n^-$')
axs[1, 0].axhline(y=0.5, color='red', linestyle='--', alpha=0.7, label=r'$y=0.5$')
axs[1, 0].set_title(r"Proyección $\pi^-$-momentum en $\hat{n}$")
axs[1, 0].set_xlabel(r"$\cos\theta_n^-$")
axs[1, 0].legend(loc="lower center")

# Gráfico 5: Proyección cos_r^-
axs[1, 1].hist(cos_r_minus, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='lightpink', label=r'$\cos\theta_n^-$')
axs[1, 1].axhline(y=0.5, color='red', linestyle='--', alpha=0.7, label=r'$y=0.5$')
axs[1, 1].set_title(r"Proyección $\pi^-$-momentum en $\hat{r}$")
axs[1, 1].set_xlabel(r"$\cos\theta_r^-$")
axs[1, 1].legend(loc="lower center")

# Gráfico 6: Proyección cos_r^-
axs[1, 2].hist(cos_k_minus, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='thistle', label=r'$\cos\theta_n^-$')
axs[1, 2].axhline(y=0.5, color='red', linestyle='--', alpha=0.7, label=r'$y=0.5$')
axs[1, 2].set_title(r"Proyección $\pi^-$-momentum en $\hat{k}$")
axs[1, 2].set_xlabel(r"$\cos\theta_k^-$")
axs[1, 2].legend(loc="lower center")


fig2, axs2 = plt.subplots(3, 3, figsize=(3, 14))
for ax in axs2.flat:
    # 1. Fijar los ticks
    ax.set_xticks(np.arange(-1.0, 1.5, 0.5))
    ax.set_yticks(np.arange(0.0, 2.5, 0.5))

    # 2. Ajustar los límites
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(0.0, 2.0)

    # 3. Encender la grilla suave que seguirá a estos ticks
    ax.grid(True, linestyle='--', alpha=0.3)

bins = 40

# Gráfico 1: cos_n^+ cos_n^-
axs2[0, 0].hist(np_nm, bins=bins, density=True, histtype='stepfilled', linewidth=1, color='lightsalmon')
axs2[0, 0].set_xlabel(r"$\cos\theta_n^+ \cos\theta_n^-$")
axs2[0, 0].set_xlim(-1.1, 1.1)

# Gráfico 2: cos_n^+ cos_r^-
axs2[0, 1].hist(np_rm, bins=bins, density=True, histtype='stepfilled', linewidth=1, color='sandybrown')
axs2[0, 1].set_xlabel(r"$\cos\theta_n^+ \cos\theta_r^-$")
axs2[0, 1].set_xlim(-1.1, 1.1)
# axs2[0, 1].grid(True)

# Gráfico 3: cos_n^+ cos_k^-
axs2[0, 2].hist(np_km, bins=bins, density=True, histtype='stepfilled', linewidth=1, color='plum')
axs2[0, 2].set_xlabel(r"$\cos\theta_n^+ \cos\theta_k^-$")
axs2[0, 2].set_xlim(-1.1, 1.1)
# axs2[0, 2].grid(True)

# Gráfico 4: cos_r^+ cos_n^-
axs2[1, 0].hist(rp_nm, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='lightsalmon')
axs2[1, 0].set_xlabel(r"$\cos\theta_r^+ \cos\theta_n^-$")
axs2[1, 0].set_xlim(-1.1, 1.1)
# axs2[1, 0].grid(True)

# Gráfico 5: cos_r^+ cos_r^-
axs2[1, 1].hist(rp_rm, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='sandybrown')
axs2[1, 1].set_xlabel(r"$\cos\theta_r^+ \cos\theta_r^-$")
axs2[1, 1].set_xlim(-1.1, 1.1)
# axs2[1, 1].grid(True)

# Gráfico 6: cos_r^+ cos_k^-
axs2[1, 2].hist(rp_km, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='plum')
axs2[1, 2].set_xlabel(r"$\cos\theta_r^+ \cos\theta_k^-$")
axs2[1, 2].set_xlim(-1.1, 1.1)
# axs2[1, 2].grid(True)

# Gráfico 7: cos_k^+ cos_n^-
axs2[2, 0].hist(kp_nm, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='lightsalmon')
axs2[2, 0].set_xlabel(r"$\cos\theta_k^+ \cos\theta_n^-$")
axs2[2, 0].set_xlim(-1.1, 1.1)
# axs2[2, 0].grid(True)

# Gráfico 8: cos_k^+ cos_r^-
axs2[2, 1].hist(kp_rm, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='sandybrown')
axs2[2, 1].set_xlabel(r"$\cos\theta_k^+ \cos\theta_r^-$")
axs2[2, 1].set_xlim(-1.1, 1.1)
# axs2[2, 1].grid(True)

# Gráfico 9: cos_k^+ cos_k^-
axs2[2, 2].hist(kp_km, bins=bins, density=True, histtype='stepfilled', linewidth=2, color='plum')
axs2[2, 2].set_xlabel(r"$\cos\theta_k^+ \cos\theta_k^-$")
axs2[2, 2].set_xlim(-1.1, 1.1)
# axs2[2, 2].grid(True)







# Ajustar automáticamente los márgenes para que los títulos y etiquetas no se superpongan
plt.tight_layout()

# Mostrar la figura con los 3 subgráficos
plt.show()






