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

    # Gauss-Seidel
    "tolerancia": 1e-6,
    "max_iter": 10000,

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
# LOCALIZAÇÃO DO POÇO 
# ============================================================

well_data = {

    'i': numerical_data['ny'] - 1,  # Índice da linha do poço
    'j': numerical_data['nx'] - 1,  # Índice da coluna do poço

    'i_inj': 0,  # Índice da linha do poço injetor
    'j_inj': 0,   # Índice da coluna do poço injetor
    'q_inj': 35, # Vazão de injeção [m³ std/dia]
    'injetor_ativo': True

}