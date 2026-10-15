import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from dados import case_data, physical_data, numerical_data, boundary_data

# Propriedades
mu = physical_data["mu"]
ct = physical_data["ct"]
k = physical_data["k_ref"]
phi = physical_data["phi"]
L = physical_data["Lx"]

# Dados (pressão prescrita nas bordas oeste e leste)
Po = case_data["Po"]
pe = boundary_data["p_leste"]
pw = boundary_data["p_oeste"]

nome_contorno = "Pressão Prescrita"

N = numerical_data.get("N_fourier", 200)

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
# GRÁFICOS (só rodam quando este arquivo é executado direto,
# não quando ele é importado pelo verificacao.py)
# ============================================================

if __name__ == "__main__":

    PASTA_PLOTS = Path(__file__).parent / "plots" / "analitico"
    PASTA_PLOTS.mkdir(parents=True, exist_ok=True)

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

        Pressao[i, :] = pressao_prescrita(x, t[i])


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

    plt.title(f"Solução Analítica - {nome_contorno}")

    fig.colorbar(
        surf,
        shrink=0.5,
        aspect=10
    )
    ax.view_init(30, 30)
    plt.savefig(
        PASTA_PLOTS / "analitico3d_pressure.png",
        dpi=150,
        bbox_inches="tight")


    # Perfis de pressão

    plt.figure(figsize=(10, 6))

    tempos_horas = [1, 6, 12, 24]

    for tempo_hora in tempos_horas:

        tempo_segundos = tempo_hora * 3600

        P = pressao_prescrita(
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
    plt.title(f"Perfis de Pressão - Solução Analítica - {nome_contorno}")
    plt.legend()
    plt.grid()
    plt.savefig(
        PASTA_PLOTS / "analitico_pressure.png",
        dpi=150,
        bbox_inches="tight")

    plt.show()