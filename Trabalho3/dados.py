# ============================================================
# GUIA DE CASOS - o que mudar para rodar cada caso
# ============================================================
#
# CASO 1 - Produtor central
#   well_data:     'i' = ny // 2, 'j' = nx // 2
#                  'produtor_ativo' = True, 'injetor_ativo' = False
#   fault_data:    'ativa' = False
#   boundary_data: 'oeste' = 'fluxo_nulo', 'leste' = 'fluxo_nulo'
#   case_data:     'tempo_horas' = 720, 'tempos_snapshot' = [72, 240, 480, 720]
#
# CASO 2 - Produtor central com falha selante
#   igual ao caso 1, mas fault_data: 'ativa' = True
#
# CASO 3 - Five-spot (1/4): injetor e produtor em cantos opostos
#   well_data:     'i' = ny - 1, 'j' = nx - 1
#                  'produtor_ativo' = True, 'injetor_ativo' = True
#                  'i_inj' = 0, 'j_inj' = 0
#   fault_data:    'ativa' = False
#   boundary_data: 'oeste' = 'fluxo_nulo', 'leste' = 'fluxo_nulo'
#   case_data:     'tempo_horas' = 720, 'tempos_snapshot' = [72, 240, 480, 720]
#
# VERIFICAÇÃO 1D - comparação com a solução analítica (Fourier)
#   physical_data: 'k_homogeneo' = True
#   well_data:     'produtor_ativo' = False, 'injetor_ativo' = False
#   fault_data:    'ativa' = False
#   boundary_data: 'oeste' = 'pressao' (150 bar), 'leste' = 'pressao' (200 bar)
#   case_data:     'tempo_horas' = 24, 'tempos_snapshot' = [1, 6, 12, 24]
#   numerical_data: 'solver' = 'direto'
#
# Obs.: rodar gerador_permeabilidade.py antes, se o CSV não existir.
# ============================================================


# ============================================================
# DADOS DO CASO
# ============================================================

case_data = {

    # Condição inicial
    "Po": 200,                  # Pressão inicial [bar]

    # Poço produtor (pressão prescrita)
    "pw": 150,                  # Pressão de fundo do poço [bar]

    "Bo": 1.25,                 # Fator volume de formação [m³/m³ std]

    # Tempo total da simulação
    "tempo_horas": 720,                # Tempo total da simulação [h]

    'tempos_snapshot': [72, 240, 480, 720]  # Tempos para snapshots [h]
}


# ============================================================
# PROPRIEDADES FÍSICAS
# ============================================================

physical_data = {

    "mu": 8e-4,                 # Viscosidade [Pa.s]
    "ct": 1.53e-9,              # Compressibilidade total [Pa^-1]
    "k_ref": 1.974e-14,         # Permeabilidade de referência [m²]
    "phi": 0.18,                # Porosidade [-]
    "h": 4,                     # Espessura [m]
    "rw": 0.10,                 # Raio do poço [m]

    # True = k_ref em todos os blocos (verificação)
    # False = campo heterogêneo lido do CSV
    "k_homogeneo": False,

    "L": 1000,                  # [m]
    "Lx": 1000,                 # Comprimento em x [m]
    "Ly": 1000                  # Comprimento em y [m]
}


# ============================================================
# DADOS NUMÉRICOS
# ============================================================

numerical_data = {

    "nx": 21,
    "ny": 21,
    "nt": 2000,

    # Solver do sistema linear: "gauss_seidel" ou "direto"
    "solver": "gauss_seidel",

    # Gauss-Seidel
    "tolerancia": 1e-6,
    "max_iter": 10000,

    # Solução analítica (número de termos da série de Fourier)
    "N_fourier": 200,

    # Gerador de permeabilidade
    "semente": 42,
    "arquivo_permeabilidade": "permeabilidade.csv"
}


# ============================================================
# FALHA SELANTE
# ============================================================

fault_data = {

    # False = caso 1 (sem falha) | True = caso 2 (com falha)
    "ativa": False,

    "coluna": 13,
    "inicio": 4,
    "fim": 16
}


# ============================================================
# CONDIÇÕES DE CONTORNO (bordas oeste e leste)
# ============================================================

boundary_data = {

    # "fluxo_nulo" = fronteira selada | "pressao" = pressão prescrita
    "oeste": "fluxo_nulo",
    "p_oeste": 150,             # Pressão na borda oeste [bar] (se "pressao")

    "leste": "fluxo_nulo",
    "p_leste": 200              # Pressão na borda leste [bar] (se "pressao")
}


# ============================================================
# LOCALIZAÇÃO DO POÇO
# ============================================================

well_data = {

    'i': numerical_data['ny'] //  2,  # Índice da linha do poço
    'j': numerical_data['nx'] // 2,  # Índice da coluna do poço
    'produtor_ativo': True,          # False = sem produtor (verificação)

    'i_inj': 0,  # Índice da linha do poço injetor
    'j_inj': 0,   # Índice da coluna do poço injetor
    'q_inj': 35, # Vazão de injeção [m³ std/dia]
    'injetor_ativo': False     # True = five-spot | False = só produtor

}