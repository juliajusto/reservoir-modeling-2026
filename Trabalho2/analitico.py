import numpy as np
import matplotlib.pyplot as plt

from dados import case_data, physical_data, numerical_data


# ============================================================
# Cálculo da difusividade hidráulica
# ============================================================

mu = physical_data["mu"]
ct = physical_data["ct"]
k = physical_data["k"]
phi = physical_data["phi"]

eta = k / (phi * mu * ct)


# ============================================================
# Função inicial
# ============================================================

def f(x):

    pe = case_data["pe"]
    pw = case_data["pw"]
    L = physical_data["L"]

    return (pe - pw) * (1 - (x / L))


# ============================================================
# Coeficientes da Série de Fourier
# ============================================================

def fourier(f, L, N):

    Bn = np.zeros(N)

    x_integracao = np.linspace(0, L, 1000)

    for n in range(1, N):

        integrando = (
            f(x_integracao)
            * np.sin(n * np.pi * x_integracao / L)
        )

        Bn[n] = (
            2 / L
            * np.trapezoid(integrando, x_integracao)
        )

    return Bn


# ============================================================
# Solução Analítica
# ============================================================

def pressao(x, t):

    L = physical_data["L"]
    N = numerical_data["N_fourier"]

    pw = case_data["pw"]
    pe = case_data["pe"]

    transiente = np.zeros_like(x)

    coef = fourier(f, L, N)

    for n in range(1, N):

        transiente += (
            coef[n]
            * np.sin(n * np.pi * x / L)
            * np.exp(
                -(n * np.pi / L)**2
                * eta
                * t
            )
        )

    permanente = pw + (pe - pw) * (x / L)

    pressao_total = permanente + transiente

    return pressao_total


# ============================================================
# Gráficos da solução analítica
# ============================================================

L = physical_data["L"]
tempo_horas_total = case_data["tempo_horas"]

x = np.linspace(0, L, 100)

t = np.linspace(
    0,
    tempo_horas_total * 3600,
    1000
)

X, T = np.meshgrid(x, t)

Pressao = pressao(X, T)


# ============================================================
# Gráfico 1 - Solução Analítica 3D
# ============================================================

fig = plt.figure(figsize=(10, 10))

ax = fig.add_subplot(
    111,
    projection="3d"
)

surf = ax.plot_surface(
    X,
    T,
    Pressao,
    rstride=2,
    cstride=1,
    cmap=plt.cm.viridis,
    linewidth=0.2,
    alpha=1
)

ax.set_xlabel("x (m)")
ax.set_ylabel("t (segundos)")
ax.set_zlabel("Pressão(x,t) (bar)")

plt.title(
    "Solução Analítica 1D da Equação da Difusividade Hidráulica"
)

fig.colorbar(
    surf,
    shrink=0.5,
    aspect=10
)

ax.view_init(30, 30)


# ============================================================
# Gráfico 2 - Perfis de Pressão
# ============================================================

plt.figure(figsize=(10, 6))

tempos_horas = [1, 6, 12, 24]

for tempo_horas in tempos_horas:

    tempo_segundos = tempo_horas * 3600

    P = pressao(
        x,
        tempo_segundos
    )

    plt.plot(
        x,
        P,
        label=f"{tempo_horas} horas"
    )

plt.xlabel("x (m)")
plt.ylabel("Pressão (bar)")

plt.title(
    "Perfis de Pressão em Diferentes Tempos"
)

plt.legend()
plt.grid()

plt.show()