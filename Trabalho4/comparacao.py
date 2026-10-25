import copy
import runpy
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
warnings.filterwarnings("ignore", message=".*non-interactive.*")
import matplotlib.pyplot as plt

# Pastas
PASTA = Path(__file__).parent
SIM_MDF = PASTA.parent / "Trabalho3" / "simulador.py"   # MDF (trabalho 3)
SIM_MVF = PASTA / "simulador.py"                         # MVF (trabalho 4)

PASTA_PLOTS = PASTA / "plots" / "comparacao"
PASTA_PLOTS.mkdir(parents=True, exist_ok=True)

# Apaga as imagens antigas da comparacao
for arquivo in PASTA_PLOTS.glob("*.png"):
    arquivo.unlink()

sys.path.insert(0, str(PASTA))
import dados

# Copia dos dados originais
NOMES = ["case_data", "physical_data", "numerical_data",
         "fault_data", "boundary_data", "well_data"]
ORIGINAL = {nome: copy.deepcopy(getattr(dados, nome)) for nome in NOMES}


def restaurar_dados():
    for nome, valores in ORIGINAL.items():
        dic = getattr(dados, nome)
        dic.clear()
        dic.update(copy.deepcopy(valores))


def configurar_caso(caso):
    nx = dados.numerical_data["nx"]
    ny = dados.numerical_data["ny"]

    dados.physical_data["k_homogeneo"] = False
    dados.boundary_data["oeste"] = "fluxo_nulo"
    dados.boundary_data["leste"] = "fluxo_nulo"
    dados.well_data["produtor_ativo"] = True

    if caso == "produtor":
        dados.well_data["i"] = ny // 2
        dados.well_data["j"] = nx // 2
        dados.well_data["injetor_ativo"] = False
        dados.fault_data["ativa"] = False

    elif caso == "falha":
        dados.well_data["i"] = ny // 2
        dados.well_data["j"] = nx // 2
        dados.well_data["injetor_ativo"] = False
        dados.fault_data["ativa"] = True

    elif caso == "fivespot":
        dados.well_data["i"] = ny - 1
        dados.well_data["j"] = nx - 1
        dados.well_data["i_inj"] = 0
        dados.well_data["j_inj"] = 0
        dados.well_data["injetor_ativo"] = True
        dados.fault_data["ativa"] = False


def rodar(arquivo):
    # roda o simulador e devolve as variaveis + tempo
    t0 = time.time()
    res = runpy.run_path(str(arquivo))
    return res, time.time() - t0


titulos = {
    "produtor": "Caso 1",
    "falha": "Caso 2",
    "fivespot": "Caso 3",
}

resultados = {}

for caso in titulos:

    print(f"\n{titulos[caso]}")

    restaurar_dados()
    configurar_caso(caso)
    mdf, t_mdf = rodar(SIM_MDF)

    restaurar_dados()
    configurar_caso(caso)
    mvf, t_mvf = rodar(SIM_MVF)

    # diferencas entre os metodos
    dP = np.abs(mdf["P_new"] - mvf["P_new"]).max() / 1e5           # bar
    dq = np.abs(np.array(mdf["vazoes"]) - np.array(mvf["vazoes"])).max()
    du = np.abs(mdf["u_mod"] - mvf["u_mod"]).max() * 86400           # m/dia

    m = mdf["m"]
    resultados[caso] = {
        "mdf": mdf, "mvf": mvf,
        "t_mdf": t_mdf, "t_mvf": t_mvf,
        "dP": dP, "dq": dq, "du": du,
        "Pbloco_mdf": mdf["Pnew"][m] / 1e5, "Pbloco_mvf": mvf["Pnew"][m] / 1e5,
    }

    print(f"max |P_MDF - P_MVF| = {dP:.2e} bar")
    print(f"max |q_MDF - q_MVF| = {dq:.2e} m3 std/dia")
    print(f"max |u_MDF - u_MVF| = {du:.2e} m/dia")

restaurar_dados()


# Tabela resumo
print("\n" + "=" * 92)
print(f"{'Caso':<8}{'Np MDF':>11}{'Np MVF':>11}{'q final MDF':>13}{'q final MVF':>13}"
      f"{'max dP (bar)':>14}{'t MDF (s)':>11}{'t MVF (s)':>11}")
print("=" * 92)
for caso, r in resultados.items():
    print(f"{titulos[caso]:<8}{r['mdf']['Np']:>11.3f}{r['mvf']['Np']:>11.3f}"
          f"{r['mdf']['vazoes'][-1]:>13.3f}{r['mvf']['vazoes'][-1]:>13.3f}"
          f"{r['dP']:>14.1e}{r['t_mdf']:>11.1f}{r['t_mvf']:>11.1f}")

r3 = resultados["fivespot"]
print(f"\npwf do injetor (Caso 3): MDF {r3['mdf']['pwf_injetor'][-1]:.2f} bar | "
      f"MVF {r3['mvf']['pwf_injetor'][-1]:.2f} bar")


# Grafico 1: vazao MDF x MVF nos 3 casos
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)

for ax, (caso, r) in zip(axes, resultados.items()):
    ax.plot(r["mdf"]["tempos"], r["mdf"]["vazoes"], linewidth=3, label="MDF")
    ax.plot(r["mvf"]["tempos"], r["mvf"]["vazoes"], "--", linewidth=2, label="MVF")
    ax.set_title(titulos[caso])
    ax.set_xlabel("Tempo (h)")
    ax.set_ylabel("Vazão (m³ std/dia)")
    ax.grid(True)
    ax.legend()

fig.suptitle("Vazão de produção: MDF x MVF")
plt.savefig(PASTA_PLOTS / "vazao_mdf_mvf.png", dpi=150, bbox_inches="tight")
plt.close()


# Grafico 2: producao acumulada MDF x MVF nos 3 casos
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)

for ax, (caso, r) in zip(axes, resultados.items()):
    ax.plot(r["mdf"]["tempos"], r["mdf"]["producoes_acumuladas"], linewidth=3, label="MDF")
    ax.plot(r["mvf"]["tempos"], r["mvf"]["producoes_acumuladas"], "--", linewidth=2, label="MVF")
    ax.set_title(titulos[caso])
    ax.set_xlabel("Tempo (h)")
    ax.set_ylabel("Produção acumulada (m³ std)")
    ax.grid(True)
    ax.legend()

fig.suptitle("Produção acumulada: MDF x MVF")
plt.savefig(PASTA_PLOTS / "producao_acumulada_mdf_mvf.png", dpi=150, bbox_inches="tight")
plt.close()


# Grafico 3: perfil de pressao na diagonal em 720 h (MDF x MVF)
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)

nx = dados.numerical_data["nx"]
dx = dados.physical_data["Lx"] / nx
s = np.sqrt(2) * (np.arange(nx) + 0.5) * dx   # distancia ao longo da diagonal

for ax, (caso, r) in zip(axes, resultados.items()):

    if caso == "fivespot":
        # injetor -> produtor
        perfil_mdf = r["mdf"]["P_new"].diagonal()
        perfil_mvf = r["mvf"]["P_new"].diagonal()
    else:
        # canto superior esquerdo -> canto inferior direito (passa pelo poco)
        perfil_mdf = np.flipud(r["mdf"]["P_new"]).diagonal()
        perfil_mvf = np.flipud(r["mvf"]["P_new"]).diagonal()

    ax.plot(s, perfil_mdf / 1e5, linewidth=2, label="MDF")
    ax.plot(s, perfil_mvf / 1e5, "o", markersize=5, markerfacecolor="none", label="MVF")
    ax.set_title(f"{titulos[caso]}  (máx. |ΔP| = {r['dP']:.1e} bar)")
    ax.set_xlabel("Distância ao longo da diagonal (m)")
    ax.set_ylabel("Pressão (bar)")
    ax.grid(True)
    ax.legend()

fig.suptitle("Perfil de pressão na diagonal em 720 h: MDF x MVF")
plt.savefig(PASTA_PLOTS / "perfil_pressao_mdf_mvf.png", dpi=150, bbox_inches="tight")
plt.close()

print(f"\nGraficos salvos em {PASTA_PLOTS}")