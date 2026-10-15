import copy
import runpy
import sys
import time
from pathlib import Path

import matplotlib

# Opções
MOSTRAR_GRAFICOS = False   # False = so salva os png
RODAR_CASOS = ["produtor", "falha", "fivespot"]
RODAR_VERIFICACAO = True

if not MOSTRAR_GRAFICOS:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt

PASTA = Path(__file__).parent
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


def rodar_script(arquivo):
    for modulo in ["simulador", "graficos", "analitico", "verificacao",
                   "gerador_permeabilidade"]:
        sys.modules.pop(modulo, None)
    runpy.run_path(str(PASTA / arquivo), run_name="__main__")
    plt.close("all")


def configurar_caso(caso):
    nx = dados.numerical_data["nx"]
    ny = dados.numerical_data["ny"]

    dados.physical_data["k_homogeneo"] = False
    dados.boundary_data["oeste"] = "fluxo_nulo"
    dados.boundary_data["leste"] = "fluxo_nulo"
    dados.well_data["produtor_ativo"] = True

    # Caso 1
    if caso == "produtor":
        dados.well_data["i"] = ny // 2
        dados.well_data["j"] = nx // 2
        dados.well_data["injetor_ativo"] = False
        dados.fault_data["ativa"] = False

    # Caso 2
    elif caso == "falha":
        dados.well_data["i"] = ny // 2
        dados.well_data["j"] = nx // 2
        dados.well_data["injetor_ativo"] = False
        dados.fault_data["ativa"] = True

    # Caso 3
    elif caso == "fivespot":
        dados.well_data["i"] = ny - 1
        dados.well_data["j"] = nx - 1
        dados.well_data["i_inj"] = 0
        dados.well_data["j_inj"] = 0
        dados.well_data["injetor_ativo"] = True
        dados.fault_data["ativa"] = False

    else:
        raise ValueError(f"Caso desconhecido: {caso}")


inicio_total = time.time()

# Campo de permeabilidade
arquivo_k = PASTA / dados.numerical_data["arquivo_permeabilidade"]

if not arquivo_k.exists():
    print("Gerando campo de permeabilidade")
    rodar_script("gerador_permeabilidade.py")

# Casos
titulos = {
    "produtor": "Caso 1 - Produtor central",
    "falha": "Caso 2 - Produtor central com falha",
    "fivespot": "Caso 3 - Five-spot",
}

for caso in RODAR_CASOS:

    print(f"\n{titulos[caso]}")

    restaurar_dados()
    configurar_caso(caso)
    rodar_script("graficos.py")

# Verificacao com o analitico
if RODAR_VERIFICACAO:

    print("\nVerificacao")

    restaurar_dados()
    rodar_script("verificacao.py")

restaurar_dados()
print(f"\nTempo total: {time.time() - inicio_total:.1f} segundos")