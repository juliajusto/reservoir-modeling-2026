from pathlib import Path

# Roda a simulação e traz as variáveis (k, P_new, tempos, vazoes, ...)
from simulador import *

# Pasta dos gráficos (separada por caso)
if falha_ativa:
    nome_caso = "falha"
elif injetor_ativo:
    nome_caso = "injetor"
else:
    nome_caso = "produtor"

PASTA_PLOTS = Path(__file__).parent / "plots" / nome_caso
PASTA_PLOTS.mkdir(parents=True, exist_ok=True)


# Gráfico do campo de permeabilidade

plt.figure(figsize=(8,6))

plt.imshow(
    k/9.869233e-16,
    origin='lower',
    extent=[0, Lx, 0, Ly],
    aspect='equal',
    cmap='viridis'
)
plt.colorbar(label='Permeabilidade (mD)')

# Poço produtor
plt.scatter(
    x=(well_data['j'] +0.5)* dx, y=(well_data['i'] +0.5)* dy,
    color='red',
    marker='x',
    s=100,
    label='Poço produtor'
)

# Poço injetor
if injetor_ativo:
    plt.scatter(
        (j_inj + 0.5)*dx, (i_inj + 0.5)*dy,
        color='white',
        marker='^',
        s=100,
        edgecolors='black',
        label='Poço injetor'
    )

plt.xlabel('Distância x (m)')
plt.ylabel('Distância y (m)')
plt.title('Campo de Permeabilidade do reservatório')
plt.legend()
plt.tight_layout()
plt.savefig(PASTA_PLOTS / "permeabilidade.png", dpi=150, bbox_inches="tight")
plt.show()


# Gráfico da Distribuição de Pressão
plt.figure(figsize=(8,6))

plt.imshow(
    P_new/1e5,
    origin='lower',
    extent=[0, Lx, 0, Ly],
    aspect='equal',
    cmap='viridis'
)

# Posição da falha (só no caso com falha)
if falha_ativa:
    x_falha = (coluna_falha + 1)* dx

    y_inicio = inicio_falha * dy
    y_fim = (fim_falha + 1) * dy

    plt.plot(
        [x_falha, x_falha],
        [y_inicio, y_fim],
        color='red',
        linewidth=3,
        label='Falha selante'
    )

plt.colorbar(label='Pressão (bar)')

# Poço produtor
plt.scatter(
    x=(well_data['j'] +0.5)* dx, y=(well_data['i'] +0.5)* dy,
    color='red',
    marker='x',
    s=100,
    label='Poço produtor'
)

# Poço injetor
if injetor_ativo:
    plt.scatter(
        (j_inj + 0.5)*dx, (i_inj + 0.5)*dy,
        color='white',
        marker='^',
        s=100,
        edgecolors='black',
        label='Poço injetor'
    )

plt.xlabel('Distância x (m)')
plt.ylabel('Distância y (m)')
plt.title(f'Distribuição de Pressão no Reservatório')
plt.legend()
plt.tight_layout()
plt.savefig(PASTA_PLOTS / "pressao.png", dpi=150, bbox_inches="tight")
plt.show()


# Gráfico da pressão em diferentes tempos (2x2)
fig, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)

# Mesma escala de cores para todos os mapas (para poder comparar)
p_min = min(P.min() for P in pressoes_snapshot) / 1e5
p_max = max(P.max() for P in pressoes_snapshot) / 1e5

for ax, P, hora in zip(axes.flat, pressoes_snapshot, horas_snapshot):

    im = ax.imshow(
        P/1e5,
        origin='lower',
        extent=[0, Lx, 0, Ly],
        aspect='equal',
        cmap='viridis',
        vmin=p_min,
        vmax=p_max
    )

    # Falha (só no caso com falha)
    if falha_ativa:
        ax.plot(
            [(coluna_falha + 1)*dx, (coluna_falha + 1)*dx],
            [inicio_falha*dy, (fim_falha + 1)*dy],
            color='red',
            linewidth=3
        )

    # Poço produtor
    ax.scatter(
        (well_data['j'] + 0.5)*dx, (well_data['i'] + 0.5)*dy,
        color='red',
        marker='x',
        s=80
    )

    # Poço injetor
    if injetor_ativo:
        ax.scatter(
            (j_inj + 0.5)*dx, (i_inj + 0.5)*dy,
            color='white',
            marker='^',
            s=80,
            edgecolors='black'
        )

    ax.set_title(f'{hora} h')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')

fig.colorbar(im, ax=axes, label='Pressão (bar)', shrink=0.8)
fig.suptitle('Evolução da Pressão no Reservatório')
plt.savefig(PASTA_PLOTS / "pressao_snapshots.png", dpi=150, bbox_inches="tight")
plt.show()


# Gráfico da vazão de produção
plt.figure(figsize=(8,6))

plt.plot(tempos, vazoes, linewidth=2)

plt.xlabel('Tempo (h)')
plt.ylabel('Vazão de Produção (m³ std/dia)')
plt.title('Vazão de Produção ao Longo do Tempo')

plt.grid(True)
plt.tight_layout()
plt.savefig(PASTA_PLOTS / "vazao.png", dpi=150, bbox_inches="tight")
plt.show()


# Gráfico da produção acumulada
plt.figure(figsize=(8,5))

plt.plot(tempos, producoes_acumuladas, linewidth=2)

plt.xlabel('Tempo (h)')
plt.ylabel('Produção Acumulada (m³ std)')
plt.title('Produção Acumulada ao longo do tempo')

plt.grid(True)
plt.tight_layout()
plt.savefig(PASTA_PLOTS / "producao_acumulada.png", dpi=150, bbox_inches="tight")
plt.show()