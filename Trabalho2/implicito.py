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

# Nome legível da condição de contorno (usado nos títulos dos gráficos)
if tipo_contorno == "pressure":
    nome_contorno = "Pressão Prescrita"
elif tipo_contorno == "flow":
    nome_contorno = "Vazão Prescrita"

# Condição inicial
P = np.zeros((nt + 1, nx))

for i in range(nx):
    P[0, i] = Po

# Equação discretizada (totalmente implícito)
for n in range(nt):
    A_matriz = np.zeros((nx, nx))
    b = np.zeros(nx)

    for i in range(1, nx - 1):
        A_matriz[i, i-1] = -beta
        A_matriz[i, i] = 1 + 2 * beta
        A_matriz[i, i+1] = -beta
        b[i] = P[n, i]

    if tipo_contorno == "pressure":
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
        A_matriz[0,0] = 1 + beta
        A_matriz[0,1] = -beta
        b[0] = P[n,0] - ((dt * q) / (phi * ct * A * dx * 1e5))  # conversão p/ bar
        # Fronteira direita
        A_matriz[-1,-1] = 1 + 4*beta
        A_matriz[-1,-2] = -4/3 * beta
        b[-1] = P[n,-1] + 8/3 * beta * Po

    P[n+1, :] = np.linalg.solve(A_matriz, b)


# Plotagem

# Gráfico 1
X = np.linspace(0, L, nx)
Y = np.linspace(0, tempo, nt + 1)

X, Y = np.meshgrid(X, Y)

fig = plt.figure(figsize=(10, 10))

ax = fig.add_subplot(
    111,
    projection='3d'
)

surf = ax.plot_surface(
    X,
    Y,
    P,
    rstride=2,
    cstride=1,
    cmap=plt.cm.viridis,
    linewidth=0.2,
    alpha=1
)

ax.set_xlabel('x (m)')
ax.set_ylabel('t (segundos)')
ax.set_zlabel('P(x,t) (bar)')

plt.title(f'Método Implícito - {nome_contorno}')

fig.colorbar(
    surf,
    shrink=0.5,
    aspect=10
)

ax.view_init(30, 30)


# Gráfico 2
plt.figure(figsize=(10, 6))

tempos_horas = [1, 6, 12, 24]

x = np.linspace(0, L, nx)

for tempo_hora in tempos_horas:

    tempo_segundos = tempo_hora * 3600

    indice = int(tempo_segundos / dt)

    P_tempo = P[indice, :]

    plt.plot(
        x,
        P_tempo,
        label=f'{tempo_hora} horas'
    )

plt.xlabel('x (m)')
plt.ylabel('Pressão (bar)')
plt.title(
    f'Perfis de Pressão em Diferentes Tempos - Método Implícito - {nome_contorno}'
)

plt.legend()
plt.grid()

plt.show()