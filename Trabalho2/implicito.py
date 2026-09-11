import numpy as np
import matplotlib.pyplot as plt

from dados import case_data, physical_data, numerical_data

# Propriedades
mu = physical_data["mu"]
ct = physical_data["ct"]
k = physical_data["k"]
phi = physical_data["phi"]
Bo = case_data["Bo"]
L = physical_data["L"]
A = physical_data["A"]

# Dados iniciais
tempo_horas = case_data["tempo_horas"]
tempo = tempo_horas * 3600

Po = case_data["Po"]
pe = case_data["pe"]
pw = case_data["pw"]

qo = case_data["qo"]
q = qo * Bo / 86400

# Malha
nx = numerical_data["nx"]
nt = numerical_data["nt"]

# Cálculos iniciais
eta = k / (phi * mu * ct)
dx = L / nx
dt = tempo / nt
beta = eta * dt / dx**2

# Escolha da condição de contorno
tipo_contorno = case_data["boundary_type"]

# Condição inicial
P = np.zeros((nt + 1, nx))

for i in range(nx):
    P[0, i] = Po

for n in range(nt):
    A_matriz = np.zeros((nx, nx))
    b = np.zeros(nx)

    for i in range(1, nx - 1):
        A_matriz[i, i-1] = -beta
        A_matriz[i, i] = 1 + 2 * beta
        A_matriz[i, i+1] = -beta
        b[i] = P[n,i]

    if tipo_contorno == "pessure":
        # Fronteira esquerda
        A_matriz[0,0] = 1 + 4 * beta
        A_matriz[0,1] = -4/3 * beta
        b[0] = P[n,0] + 8/3 * beta * pw
        # Fronteira direita
        A_matriz[-1,-1] = 1 + 4 * beta
        A_matriz[-1,-2] = -4/3 * beta
        b[-1] = P[n,-1] + 8/3 * beta * pe

    elif tipo_contorno == "flow":
        # Fronteira esquerda
        A_matriz[0,0]= 1 + beta
        A_matriz[0,1] = -beta
        b[0] = P[n,0] +  ((dt * q) / (phi * ct * A * dx))
        # Fronteira direita 
        A_matriz[-1,-1] = 1 + 4*beta
        A_matriz[-1,-2] = -4/3 * beta
        b[-1] = P[n,-1] + 8/3 * beta * pe

    P[n+1, :] = np.linalg.solve(A, b)