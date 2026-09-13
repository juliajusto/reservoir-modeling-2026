import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

from dados import case_data, physical_data, numerical_data


# Propriedades
mu = physical_data["mu"]
ct = physical_data["ct"]
k = physical_data["k"]
phi = physical_data["phi"]
L = physical_data["L"]
A = physical_data["A"]

# Dados
Po = case_data["Po"]
pe = case_data["pe"]
pw = case_data["pw"]
qo = case_data["qo"]
Bo = case_data["Bo"]

tipo_contorno = case_data["boundary_type"]

N = numerical_data["N_fourier"]

# Cálculo da difusividade
eta = k / (phi * mu * ct)


# ============================================================
# PRESSÃO PRESCRITA
# ============================================================

def f(x):
    return (pe - pw) * (1 - x / L)


def fourier(f, L, N):

    Bn = np.zeros(N)

    x_integracao = np.linspace(0, L, 1000)

    for n in range(1, N):

        integrando = (
            f(x_integracao)
            * np.sin(n * np.pi * x_integracao / L)
        )

        Bn[n] = 2 / L * np.trapezoid(
            integrando,
            x_integracao
        )

    return Bn


def pressao_prescrita(x, t):

    transiente = np.zeros_like(x)

    coef = fourier(f, L, N)

    for n in range(1, N):

        transiente += (
            coef[n]
            * np.sin(n * np.pi * x / L)
            * np.exp(
                -(n * np.pi / L)**2
                * eta * t
            )
        )

    permanente = pw + (pe - pw) * (x / L)

    return permanente + transiente


# ============================================================
# VAZÃO PRESCRITA
# ============================================================

def pressao_vazao(x, t):

    if t == 0:
        return np.full_like(x, float(Po))

    q = qo * Bo / 86400

    delta_p_pa = (
        (q * mu * L) / (k * A)
        * (
            np.sqrt(
                (4 * eta * t) / (np.pi * L**2)
            )
            * np.exp(
                -x**2 / (4 * eta * t)
            )
            - (x / L)
            * erfc(
                x / np.sqrt(4 * eta * t)
            )
        )
    )

    # delta_p_pa está em Pascal; convertendo para bar
    return Po - delta_p_pa / 1e5


# solução
def pressao(x, t):

    if tipo_contorno == "pressure":
        return pressao_prescrita(x, t)

    elif tipo_contorno == "flow":
        return pressao_vazao(x, t)
#inicio dos graficos

x = np.linspace(0, L, 100)

tempo_horas = case_data["tempo_horas"]

t = np.linspace(
    0,
    tempo_horas * 3600,
    1000
)

X, T = np.meshgrid(x, t)

Pressao = np.zeros_like(X)

for i in range(len(t)):

    Pressao[i, :] = pressao(x, t[i])


# Gráfico 3D

fig = plt.figure(figsize=(10, 10))

ax = fig.add_subplot(111, projection="3d")

surf = ax.plot_surface(
    X,
    T,
    Pressao,
    rstride=2,
    cstride=1,
    cmap=plt.cm.viridis,
    linewidth=0.2
)

ax.set_xlabel("x (m)")
ax.set_ylabel("t (segundos)")
ax.set_zlabel("Pressão (bar)")

plt.title("Solução Analítica")

fig.colorbar(
    surf,
    shrink=0.5,
    aspect=10
)
ax.view_init(30, 30)


# Perfis de pressão

plt.figure(figsize=(10, 6))

tempos_horas = [1, 6, 12, 24]

for tempo_hora in tempos_horas:

    tempo_segundos = tempo_hora * 3600

    P = pressao(
        x,
        tempo_segundos
    )

    plt.plot(
        x,
        P,
        label=f"{tempo_hora} horas"
    )

plt.xlabel("x (m)")
plt.ylabel("Pressão (bar)")
plt.title("Perfis de Pressão - Solução Analítica")
plt.legend()
plt.grid()

plt.show()