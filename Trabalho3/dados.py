# Casos (o main.py ja configura sozinho)
# Caso 1: produtor no centro, sem falha, sem injetor
# Caso 2: produtor no centro, com falha
# Caso 3: five-spot, injetor em (0,0) e produtor em (ny-1, nx-1)
# Verificacao: k homogeneo, sem pocos, pressao prescrita oeste/leste

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

    "k_homogeneo": False,       # True = k_ref em todos os blocos

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

    # "gauss_seidel" ou "direto"
    "solver": "gauss_seidel",

    # Gauss-Seidel
    "tolerancia": 1e-6,
    "max_iter": 10000,

    # Termos da serie de Fourier
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
# CONDIÇÕES DE CONTORNO
# ============================================================

boundary_data = {

    # "fluxo_nulo" ou "pressao"
    "oeste": "fluxo_nulo",
    "p_oeste": 150,             # [bar]

    "leste": "fluxo_nulo",
    "p_leste": 200              # [bar]
}


# ============================================================
# LOCALIZAÇÃO DO POÇO
# ============================================================

well_data = {

    'i': numerical_data['ny'] //  2,  # Índice da linha do poço
    'j': numerical_data['nx'] // 2,  # Índice da coluna do poço
    'produtor_ativo': True,

    'i_inj': 0,  # Índice da linha do poço injetor
    'j_inj': 0,   # Índice da coluna do poço injetor
    'q_inj': 35, # Vazão de injeção [m³ std/dia]
    'injetor_ativo': False     # True = five-spot | False = só produtor

}