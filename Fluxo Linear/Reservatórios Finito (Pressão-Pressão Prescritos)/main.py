import numpy as np
import math
import matplotlib.pyplot as plt

# Propriedades do reservatório
mu = 8e-4 # Pa.s
ct = 1.53e-9 # Pa^-1
k = 1.974e-14 # m^2
phi = 0.18
pe = 200 # bar
pw = 150 # bar
L = 1000 # m
N = 200

# Cálculos Iniciais
eta = k/(phi*mu*ct)

# Define a função desejada:
def f(x):
  return (pe-pw)*(1- (x/L))

# Coeficiente da Série de Fourier
def fourier(f, L, N):
  Bn = np.zeros(N)

  for n in range(1, N):
    def integrando(x):
      return f(x)*np.sin(n*np.pi*x/L)

    Bn[n] = 2/L*np.trapezoid(integrando(np.linspace(0,L,1000)),np.linspace(0,L,1000))

  return Bn

# Solução da Equação hidraulica:
def pressao(x, t, L, N):
  transiente = np.zeros_like(x) # np.zeros: cria uma matriz de zeros com o mesmo formato e tipo de dados de uma matriz de entrada fornecida.
  coef = fourier(f, L, N)

  for n in range(1,N):

    transiente += coef[n]*np.sin(n*np.pi*x/L)*np.exp(-(n*np.pi/L)**2*eta*t)

  permanente = pw + (pe-pw)*(x/L)

  pressao_total = permanente + transiente

  return pressao_total


# Define os valores de x e t
x = np.linspace(0, L, 100)
t = np.linspace(0, 864000, 1000)
X, T = np.meshgrid(x,t)

Pressao = pressao(X, T, L, N)

# Gráfico 1
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, T, Pressao, rstride=2, cstride=1, cmap=plt.cm.viridis, linewidth=0.2, alpha=1)
ax.set_xlabel('x (cm)')
ax.set_ylabel('t (segundos)')
ax.set_zlabel('Pressão(x,t) (bar)')
plt.title('Solução Analítica 1D da Equação da Difusividade Hidráulica')
fig.colorbar(surf, shrink=0.5, aspect=10)
ax.view_init(30, 30)


# Gráfico 2
plt.figure(figsize=(10, 6))
tempos_dias = [0, 0.1, 1, 5, 10]

for tempo_dias in tempos_dias:

    tempo_segundos = tempo_dias * 24 * 3600
    P = pressao(x, tempo_segundos, L, N)
    plt.plot(x, P, label=f'{tempo_dias} dias')
plt.xlabel('x (m)')
plt.ylabel('Pressão (bar)')
plt.title('Perfis de Pressão em Diferentes Tempos')
plt.legend()
plt.grid()
plt.show()