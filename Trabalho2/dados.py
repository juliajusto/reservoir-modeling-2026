# ============================================================
# DADOS DO CASO
# ============================================================

case_data = {

    # Condição inicial
    "Po": 200,                  # Pressão inicial [bar]

    # Pressão prescrita
    "pe": 200,                  # Pressão na extremidade direita [bar]
    "pw": 150,                  # Pressão na extremidade esquerda [bar]

    # Vazão prescrita
    "qo": 35,                   # Vazão de óleo [m³ std/d]
    "Bo": 1.25,                 # Fator volume de formação [m³/m³ std]

    # Escolha da condição de contorno
    # "pressure" = pressão prescrita
    # "flow" = vazão prescrita
    "boundary_type": "flow",

    # Tempo total da simulação
    "tempo_horas": 24
}


# ============================================================
# PROPRIEDADES FÍSICAS
# ============================================================

physical_data = {

    "mu": 8e-4,
    "ct": 1.53e-9,
    "k": 1.974e-14,
    "phi": 0.18,

    "L": 1000,
    "A": 92.9
}


# ============================================================
# DADOS NUMÉRICOS
# ============================================================

numerical_data = {

    "nx": 100,
    "nt": 2000,
    "N_fourier": 200
}