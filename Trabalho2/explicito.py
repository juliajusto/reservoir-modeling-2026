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


# Análise de estabilidade e convergência
if dt <= (dx**2) / (4 * eta):
    print("Método Convergente.")
else:
    print("Método Não Convergente. Refinar a malha do tempo")


# Solução da Equação

# Condição inicial
P = np.zeros((nt + 1, nx))

for i in range(nx):
    P[0, i] = Po


# Equação discretizada
for n in range(nt):

    for i in range(nx):

        if tipo_contorno == "pressure":

            if i == 0:
                # Fronteira esquerda
                P[n+1, i] = (
                    8/3 * beta * pw
                    + (1 - 4*beta) * P[n, i]
                    + 4/3 * beta * P[n, i+1]
                )

            elif i == nx - 1:
                # Fronteira direita
                P[n+1, i] = (
                    4/3 * beta * P[n, i-1]
                    + (1 - 4*beta) * P[n, i]
                    + 8/3 * beta * pe
                )

            else:
                P[n+1, i] = (
                    beta * P[n, i-1]
                    + (1 - 2*beta) * P[n, i]
                    + beta * P[n, i+1]
                )

        elif tipo_contorno == "flow":

            if i == 0:
                # Fronteira esquerda
                P[n+1, i] = (
                    (1 - beta) * P[n, i]
                    + beta * P[n, i+1]
                    - beta * ((q * mu * dx) / (k * A) / 1e5)
                )

            elif i == nx - 1:
                # Fronteira direita
                P[n+1, i] = (
                    4/3 * beta * P[n, i-1]
                    + (1 - 4*beta) * P[n, i]
                    + 8/3 * beta * Po
                )

            else:
                # Nós internos
                P[n+1, i] = (
                    beta * P[n, i-1]
                    + (1 - 2*beta) * P[n, i]
                    + beta * P[n, i+1]
                )


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

plt.title(f'Método Explícito - {nome_contorno}')

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
    f'Perfis de Pressão em Diferentes Tempos - Método Explícito - {nome_contorno}'
)

plt.legend()
plt.grid()

plt.show()