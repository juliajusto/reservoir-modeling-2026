import numpy as np
from pathlib import Path

from dados import physical_data, numerical_data

# Dimensões da malha
nx = numerical_data["nx"]
ny = numerical_data["ny"]

# Permeabilidade de referência
k_ref = physical_data["k_ref"]  # m²

# Gerador de números aleatórios
rng = np.random.default_rng(numerical_data["semente"])

# Campo heterogêneo: 50% a 150% do valor de referência
k = k_ref * rng.uniform(0.5, 1.5, size=(ny, nx))

# Salvar o campo de permeabilidade em CSV
np.savetxt(
    Path(__file__).parent / numerical_data["arquivo_permeabilidade"],
    k,
    delimiter=','
)

print('Arquivo CSV gerado')