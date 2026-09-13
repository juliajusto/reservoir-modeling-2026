import matplotlib.pyplot as plt
import numpy as np

from dados import case_data, physical_data, numerical_data
from analitico import pressao_prescrita, pressao_vazao

import explicito
import implicito

L = physical_data["L"]
nx = numerical_data["nx"]

tempo_horas = case_data["tempo_horas"]
tempo = tempo_horas * 3600

tipo_contorno = case_data["boundary_type"]

x = np.linspace(0, L, nx)

# último passo de tempo de cada método
P_explicito = explicito.P[-1, :]
P_implicito = implicito.P[-1, :]

# comparando analítico x implícito x explícito

if tipo_contorno == "pressure":
    P_analitico = pressao_prescrita(x, tempo)
    titulo = "Pressão Prescrita"
    nome_arquivo = "comparacao_pressure.png"

elif tipo_contorno == "flow":
    P_analitico = pressao_vazao(x, tempo)
    titulo = "Vazão Prescrita"
    nome_arquivo = "comparacao_flow.png"

erro_explicito = np.abs(P_explicito - P_analitico)
erro_implicito = np.abs(P_implicito - P_analitico)

erro_medio_exp = np.mean(erro_explicito)
erro_max_exp = np.max(erro_explicito)
erro_medio_imp = np.mean(erro_implicito)
erro_max_imp = np.max(erro_implicito)

fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(9, 8.5),
    sharex=True,
    gridspec_kw={"height_ratios": [2, 1]}
)

# --- Painel 1: perfis de pressão ---
ax1.plot(x, P_explicito, label="Explícito", linewidth=2)
ax1.plot(x, P_implicito, label="Implícito", linewidth=2)
ax1.plot(x, P_analitico, '--', label="Analítico", linewidth=2, color='black')

# Linha de referência em P=0 (útil no caso de vazão, onde a pressão pode ficar negativa)
if tipo_contorno == "flow":
    ax1.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.6)

ax1.set_ylabel("Pressão (bar)")
ax1.set_title(f"Comparação dos Métodos - {titulo}")
ax1.legend()
ax1.grid(True, alpha=0.3)

# --- Painel 2: erro absoluto ao longo de x ---
ax2.plot(x, erro_explicito, label="Erro Explícito", linewidth=1.5)
ax2.plot(x, erro_implicito, label="Erro Implícito", linewidth=1.5)

ax2.set_xlabel("x (m)")
ax2.set_ylabel("Erro absoluto (bar)")
ax2.legend()
ax2.grid(True, alpha=0.3)

# Texto de erro médio e máximo abaixo dos gráficos, dentro da figura
texto_erro = (
    f"Erro médio (Explícito): {erro_medio_exp:.4f} bar   |   "
    f"Erro máximo (Explícito): {erro_max_exp:.4f} bar\n"
    f"Erro médio (Implícito): {erro_medio_imp:.4f} bar   |   "
    f"Erro máximo (Implícito): {erro_max_imp:.4f} bar"
)

fig.text(
    0.5, 0.005, texto_erro,
    ha='center', va='bottom',
    fontsize=9,
    bbox=dict(boxstyle="round", facecolor="whitesmoke", alpha=0.9, edgecolor="gray")
)

plt.tight_layout(rect=[0, 0.07, 1, 1])  # reserva espaço embaixo para o texto
plt.savefig(nome_arquivo, dpi=150)
plt.show()
plt.close(fig)


print(f"=== Caso: {tipo_contorno} ===")
print(f"Erro médio  (Explícito): {erro_medio_exp:.4f} bar")
print(f"Erro máximo (Explícito): {erro_max_exp:.4f} bar")
print(f"Erro médio  (Implícito): {erro_medio_imp:.4f} bar")
print(f"Erro máximo (Implícito): {erro_max_imp:.4f} bar")