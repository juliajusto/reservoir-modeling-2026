import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

mu = 8e-4          # Pa.s
ct = 1.53e-9       # Pa^-1
k = 1.974e-14      # m^2
phi = 0.18
h = 10             # m
xf = 100            # m
A = 2*h*xf         # m^2
qo = 30             # m^3std/d
Bo = 1.2            # m^3/m^3std
qw = qo*Bo          # m^3/d

# Conversão da vazão para SI
qw = qw/(24*3600)  # m^3/s

# Pressão inicial
p0 = 200e5         # Pa

# Cálculos iniciais
eta = k/(phi*mu*ct)

# Coeficiente da solução
coef = qw*mu/(k*A)

print(f"eta = {eta:.4f} m²/s")
print(f"qw = {qw:.6e} m³/s")
print(f"A = {A:.0f} m²")
print(f"coef = {coef:.2f} Pa/m")


def pressao(x, t):
    termo_A = np.sqrt((4*eta*t)/(np.pi))*np.exp(-x**2/(4*eta*t))
    termo_B = x*erfc(x/(np.sqrt(4*eta*t)))
    p = p0 - coef*(termo_A - termo_B)
    return p


#Gráfico 1 - Perfis de pressão em diferentes tempos

x = np.linspace(0, 2000, 200)

tempos_perfis = [0.1, 1, 5, 10, 30]

plt.figure(figsize=(10, 6))

for tempo_dias in tempos_perfis:
    tempo_segundos = tempo_dias * 24 * 3600

    P = pressao(x, tempo_segundos)

    if tempo_dias == 0.1:
        legenda = '0,1 dia'
    else:
        legenda = f'{tempo_dias} dias'

    plt.plot(x, P/1e5, label=legenda)

plt.xlabel('x (m)')
plt.ylabel('Pressão (bar)')
plt.title('Perfis de Pressão em Diferentes Tempos')
plt.legend()
plt.grid()
plt.show()


#Gráfico 2 - Queda de pressão em x = 0

tempos_log = np.logspace(-1, np.log10(30), 100)

delta_p = []

for tempo_dias in tempos_log:
    tempo_segundos = tempo_dias * 24 * 3600

    P = pressao(0, tempo_segundos)

    delta_p.append((p0 - P)/1e5)

plt.figure(figsize=(10, 6))

plt.loglog(tempos_log, delta_p)

plt.xlabel('Tempo (dias)')
plt.ylabel('Δp (bar)')
plt.title('Queda de Pressão em x = 0 em Função do Tempo')
plt.grid()
plt.show()


#Gráfico 3 - Superfície 3D da pressão

t = np.linspace(3600, 30*24*3600, 100)

X, T = np.meshgrid(x, t)

P = pressao(X, T)/1e5

fig = plt.figure(figsize=(10, 6))

ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(X, T/(24*3600), P,
                       rstride=2, cstride=1,
                       cmap=plt.cm.viridis,
                       linewidth=0.2,
                       alpha=1)

ax.set_xlabel('x (m)')
ax.set_ylabel('Tempo (dias)')
ax.set_zlabel('Pressão (bar)')

plt.title('Superfície de Pressão em Função do Espaço e Tempo')

fig.colorbar(surf, shrink=0.5, aspect=10)

ax.view_init(30, 30)
plt.show()

