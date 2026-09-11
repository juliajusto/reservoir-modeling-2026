import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.special as sp

# Propriedades do reservatório
k = 100 * 9.869233e-16 # 100 md -> m²
mu = 3.0e-3  # 3 cp -> Pa.s
phi = 0.20
ct = 130e-6 / 98066.5 # (kgf/cm²)^-1 -> Pa^-1
h = 4       # m
Bo = 1.25   # m³/m³ std
qo = 35   # m³ std/d
qw = qo * Bo / 86400  # m³/s
rw = 0.10  # m

# Pressão inicial adotada
p0_bar = 200
p0 = p0_bar * 1e5            # Pa

# Cálculos Iniciais
eta = k/(phi*mu*ct)

def pressao_radial(r, t):
    u = r**2 / (4 * eta * t)
    E1 = -sp.expi(-u)
    p = p0 - (qw * mu / (4 * np.pi * k * h)) * E1

    return p

# Define o dominio
r = np.linspace(rw, 1000, 100)
t = np.linspace(1, 30 * 24 * 3600, 500)
R, T = np.meshgrid(r, t)
Pressao = pressao_radial(R, T)
Pressao_bar = Pressao / 1e5

# Gráfico 1
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(R,T,Pressao_bar,rstride=2,cstride=1,cmap=plt.cm.viridis,linewidth=0.2,alpha=1)
ax.set_xlabel('r (m)')
ax.set_ylabel('t (segundos)')
ax.set_zlabel('Pressão(r,t) (bar)')
plt.title('Solução Analítica Radial da Equação da Difusividade Hidráulica')
fig.colorbar(surf, shrink=0.5, aspect=10)
ax.view_init(30, 30)

# Gráfico 2
plt.figure(figsize=(10, 6))
tempos_dias = [0.1, 1, 10, 30]
for tempo_dias in tempos_dias:

    tempo_segundos = tempo_dias * 24 * 3600
    P = pressao_radial(r, tempo_segundos)
    P_bar = P / 1e5
    plt.plot(r, P_bar, label=f'{tempo_dias} dias')
plt.xlabel('r (m)')
plt.ylabel('Pressão (bar)')
plt.title('Perfis de Pressão em Diferentes Tempos')
plt.legend()
plt.grid()
plt.show()