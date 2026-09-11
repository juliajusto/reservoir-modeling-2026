# ============================================================
# DADOS DO CASO
# ============================================================

case_data = {

    # Condição inicial e condições de contorno
    "Po": 200,                  # Pressão inicial [bar]
    "pe": 200,                  # Pressão na extremidade direita [bar]
    "pw": 150,                  # Pressão na extremidade esquerda [bar]

    # Dados de produção
    "qo": 35,                   # Vazão de óleo [m³ std/d]
    "Bo": 1.25,                 # Fator volume de formação [m³/m³ std]

    # Condição de contorno
    # "pressure" = pressão prescrita
    # "flow" = vazão prescrita
    "boundary_type": "pressure",

    # Tempo total da simulação
    "tempo_horas": 24
}


# ============================================================
# PROPRIEDADES FÍSICAS
# ============================================================

physical_data = {

    "mu": 8e-4,                 # Viscosidade [Pa.s]
    "ct": 1.53e-9,              # Compressibilidade total [Pa^-1]
    "k": 1.974e-14,             # Permeabilidade [m²]
    "phi": 0.18,                # Porosidade

    "L": 1000,                  # Comprimento do reservatório [m]
    "A": 92.9                   # Área da seção transversal [m²]
}


# ============================================================
# DADOS NUMÉRICOS
# ============================================================

numerical_data = {

    "nx": 100,                  # Número de nós na direção x
    "nt": 2000,                 # Número de passos de tempo
    "N_fourier": 200            # Número de termos da série de Fourier
}