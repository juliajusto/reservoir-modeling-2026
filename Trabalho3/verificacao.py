import importlib
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import dados
from dados import (
    case_data,
    physical_data,
    numerical_data,
    boundary_data,
    fault_data,
    well_data
)

# Caso de verificacao 1D
physical_data["k_homogeneo"] = True
well_data["produtor_ativo"] = False
well_data["injetor_ativo"] = False
fault_data["ativa"] = False
boundary_data["oeste"] = "pressao"
boundary_data["leste"] = "pressao"
case_data["tempo_horas"] = 24
case_data["tempos_snapshot"] = [1, 6, 12, 24]
numerical_data["solver"] = "direto"

from analitico import pressao_prescrita

PASTA_PLOTS = Path(__file__).parent / "plots" / "verificacao"
PASTA_PLOTS.mkdir(parents=True, exist_ok=True)


def comparar_pressao(x_numerico, pressao_numerica, tempo_horas):

    tempo_segundos = tempo_horas * 3600

    pressao_analitica = pressao_prescrita(
        x_numerico,
        tempo_segundos
    )

    erro_absoluto = np.abs(
        pressao_numerica - pressao_analitica
    )

    return pressao_analitica, erro_absoluto


def rodar_simulador(nx, ny, nt):

    numerical_data["nx"] = nx
    numerical_data["ny"] = ny
    numerical_data["nt"] = nt

    import simulador
    importlib.reload(simulador)

    return simulador


# Comparacao numerico x analitico

nx_base = dados.numerical_data["nx"]
ny_base = dados.numerical_data["ny"]
nt_base = dados.numerical_data["nt"]

sim = rodar_simulador(nx_base, ny_base, nt_base)

x_numerico = (np.arange(sim.nx) + 0.5) * sim.dx   # centro dos blocos
x_analitico = np.linspace(0, sim.Lx, 500)

fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(9, 8.5),
    sharex=True,
    gridspec_kw={"height_ratios": [2, 1]}
)

cores = plt.cm.viridis(np.linspace(0, 0.85, len(sim.horas_snapshot)))

print("\nNumerico x analitico")

for P, hora, cor in zip(sim.pressoes_snapshot, sim.horas_snapshot, cores):

    # linha do meio
    pressao_numerica = P[sim.ny // 2, :] / 1e5

    pressao_analitica, erro = comparar_pressao(
        x_numerico,
        pressao_numerica,
        hora
    )

    print(f"{hora} h - erro medio: {erro.mean():.4f} bar, erro maximo: {erro.max():.4f} bar")

    ax1.plot(x_analitico, pressao_prescrita(x_analitico, hora * 3600),
             color=cor, linewidth=1.5)
    ax1.plot(x_numerico, pressao_numerica, 'o',
             color=cor, markersize=4, label=f'{hora} h')

    ax2.plot(x_numerico, erro, 'o-', color=cor, markersize=3, label=f'{hora} h')

dif_y = np.ptp(sim.pressoes_snapshot[-1], axis=0).max() / 1e5
print(f"Diferenca entre as linhas em y: {dif_y:.2e} bar")

ax1.plot([], [], 'k-', label='Analítico')
ax1.plot([], [], 'ko', markersize=4, label='Numérico')
ax1.set_ylabel("Pressão (bar)")
ax1.set_title(f"Numérico x Analítico (malha {sim.nx}x{sim.ny})")
ax1.legend()
ax1.grid(True, alpha=0.3)

ax2.set_xlabel("x (m)")
ax2.set_ylabel("Erro absoluto (bar)")
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(PASTA_PLOTS / "comparacao_analitico.png", dpi=150, bbox_inches="tight")
plt.show()


# Convergencia no espaco (ny = 3 pq o caso e 1D)

lista_nx = [11, 21, 41, 81]
nt_conv = 2000
lista_dx = []
erros_rms_x = []

print("\nConvergencia no espaco")

for nx in lista_nx:

    sim = rodar_simulador(nx, 3, nt_conv)

    x_numerico = (np.arange(sim.nx) + 0.5) * sim.dx
    pressao_numerica = sim.pressoes_snapshot[-1][1, :] / 1e5

    _, erro = comparar_pressao(x_numerico, pressao_numerica, 24)

    lista_dx.append(sim.dx)
    erros_rms_x.append(np.sqrt(np.mean(erro**2)))

    print(f"nx = {nx}, dx = {sim.dx:.1f} m, erro RMS = {erros_rms_x[-1]:.4f} bar")

lista_dx = np.array(lista_dx)
erros_rms_x = np.array(erros_rms_x)

ordem_x = np.log(erros_rms_x[:-1] / erros_rms_x[1:]) / np.log(lista_dx[:-1] / lista_dx[1:])
print("Ordem no espaco:", np.round(ordem_x, 2))


# Convergencia no tempo (referencia com dt bem pequeno)

nx_t = 41
lista_nt = [24, 48, 96, 192, 384]

sim_ref = rodar_simulador(nx_t, 3, 6144)
P_ref = sim_ref.pressoes_snapshot[-1][1, :] / 1e5

lista_dt = []
erros_rms_t = []

print("\nConvergencia no tempo")

for nt in lista_nt:

    sim = rodar_simulador(nx_t, 3, nt)
    erro = np.abs(sim.pressoes_snapshot[-1][1, :] / 1e5 - P_ref)

    lista_dt.append(sim.dt)
    erros_rms_t.append(np.sqrt(np.mean(erro**2)))

    print(f"nt = {nt}, dt = {sim.dt:.1f} s, erro RMS = {erros_rms_t[-1]:.4f} bar")

lista_dt = np.array(lista_dt)
erros_rms_t = np.array(erros_rms_t)

ordem_t = np.log(erros_rms_t[:-1] / erros_rms_t[1:]) / np.log(lista_dt[:-1] / lista_dt[1:])
print("Ordem no tempo:", np.round(ordem_t, 2))


# Estabilidade - dt grande (12 h)

sim_est = rodar_simulador(81, 3, 2)

eta_dif = physical_data["k_ref"] / (physical_data["phi"] * physical_data["mu"] * physical_data["ct"])
dt_explicito = sim_est.dx**2 / (4 * eta_dif)   # criterio do explicito 2D

P_est = sim_est.pressoes_snapshot[-1][1, :] / 1e5
x_est = (np.arange(sim_est.nx) + 0.5) * sim_est.dx
_, erro_est = comparar_pressao(x_est, P_est, 24)

dentro_limites = (P_est.min() >= 150 - 1e-6) and (P_est.max() <= 200 + 1e-6)
monotono = np.all(np.diff(P_est) >= -1e-9)

print("\nEstabilidade (nx = 81, dt = 12 h)")
print(f"Limite do explicito: dt <= {dt_explicito:.1f} s")
print(f"dt usado: {sim_est.dt:.0f} s ({sim_est.dt/dt_explicito:.0f} vezes o limite)")
print(f"Pressao entre 150 e 200 bar: {dentro_limites}, sem oscilacao: {monotono}")
print(f"Erro RMS em 24 h: {np.sqrt(np.mean(erro_est**2)):.3f} bar")


# Grafico de convergencia

fig, (axa, axb) = plt.subplots(1, 2, figsize=(12, 5))

# Espaco
axa.loglog(lista_dx, erros_rms_x, 'o-', linewidth=2, label='Erro numérico')
ref2 = erros_rms_x[-1] * (lista_dx / lista_dx[-1])**2
axa.loglog(lista_dx, ref2, 'k--', label='Inclinação 2')
axa.set_xlabel("Δx (m)")
axa.set_ylabel("Erro RMS (bar)")
axa.set_title(f"Convergência no espaço (ordem ≈ {ordem_x.mean():.2f})")
axa.legend()
axa.grid(True, which='both', alpha=0.3)

# Tempo
axb.loglog(lista_dt, erros_rms_t, 's-', color='tab:orange', linewidth=2, label='Erro numérico')
ref1 = erros_rms_t[-1] * (lista_dt / lista_dt[-1])**1
axb.loglog(lista_dt, ref1, 'k--', label='Inclinação 1')
axb.set_xlabel("Δt (s)")
axb.set_ylabel("Erro RMS (bar)")
axb.set_title(f"Convergência no tempo (ordem ≈ {ordem_t.mean():.2f})")
axb.legend()
axb.grid(True, which='both', alpha=0.3)

plt.tight_layout()
plt.savefig(PASTA_PLOTS / "convergencia.png", dpi=150, bbox_inches="tight")
plt.show()