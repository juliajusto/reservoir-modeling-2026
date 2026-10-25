import numpy as np
import matplotlib.pyplot as plt
import time
from pathlib import Path

from dados import case_data, physical_data, numerical_data, fault_data, well_data, boundary_data

# Propriedades
mu = physical_data["mu"]
ct = physical_data["ct"]
# k = 1.974e-14  # m² - referência para o campo de permeabilidades
phi = physical_data["phi"]
Bo = case_data["Bo"]
L = physical_data["L"]
h = physical_data["h"]
rw = physical_data["rw"]

# Dados Iniciais
tempo_horas = case_data["tempo_horas"]
tempo = tempo_horas * 3600
Po = case_data["Po"] * 1e5  # Pa
pw = case_data["pw"] * 1e5  # Pa
tempos_snapshot = case_data["tempos_snapshot"]

# Dados do Injetor
i_inj = well_data['i_inj']
j_inj = well_data['j_inj']
q_inj = well_data['q_inj']
injetor_ativo = well_data['injetor_ativo']

# Produtor ligado/desligado
produtor_ativo = well_data['produtor_ativo']

# Permeabilidade homogênea ou heterogênea
k_ref = physical_data["k_ref"]
k_homogeneo = physical_data["k_homogeneo"]

# Solver do sistema linear
solver = numerical_data["solver"]

# Condições de contorno oeste e leste
contorno_oeste = boundary_data["oeste"]
contorno_leste = boundary_data["leste"]
p_oeste = boundary_data["p_oeste"] * 1e5  # Pa
p_leste = boundary_data["p_leste"] * 1e5  # Pa
q_oeste = boundary_data["q_oeste"] 
q_leste = boundary_data["q_leste"] 

#Conversão da vazão de injeção para m³/s
q_inj = well_data['q_inj'] * Bo / 86400  # m³/s

# Parâmetros da malha
nx = numerical_data["nx"]
ny = numerical_data["ny"]
N = nx * ny
nt = numerical_data["nt"]

#vazao das bordas por bloco (m3/s no reservaotiro
q_oeste = q_oeste * Bo / 86400 / ny
q_leste = q_leste * Bo / 86400 / ny

# Dimensão do reservatório
Lx = physical_data["Lx"]
Ly = physical_data["Ly"]

# Discretização espacial e temporal
dx = Lx / nx
dy = Ly / ny
dt = tempo / nt

# Cálculos Iniciais
V = dx * dy * h  # Volume do bloco [m³]
Ax = dy * h  # Área da face em x [m²]
Ay = dx * h  # Área da face em y [m²]
Az = dx * dy  # Área da face em z [m²]

acumulacao = (phi * ct* V) / dt

# Importa o campo de permeabilidade
if k_homogeneo:
    k = np.full((ny, nx), k_ref)
else:
    k = np.loadtxt(
        Path(__file__).parent / numerical_data["arquivo_permeabilidade"],
        delimiter=','
    )

# verifica dimensões
assert k.shape == (ny, nx)
# Condição inicial
P_old = np.full((ny, nx), Po)
P_new = np.full((ny, nx), Po)
Pold = np.full(N, Po)

# Paramêtros do método de Gauss-Seidel
tolerancia = numerical_data["tolerancia"]
max_iter = numerical_data["max_iter"]

# Armazenamento dos resultados
tempos = []
vazoes = []
producoes_acumuladas = []
pressoes_snapshot = []
horas_snapshot = []
pwf_injetor = []

Np = 0 # Produção acumulada em m³ std

# Inicializando as Matrizes
A = np.zeros((N, N))
b = np.zeros(N)
tempo_atual = 0
passo = 0

# Índices da matriz 2D para 1D
ind = np.arange(N).reshape(ny, nx)

# Falha selante Vertical
falha_ativa = fault_data["ativa"]
coluna_falha = fault_data["coluna"]
inicio_falha = fault_data["inicio"]
fim_falha = fault_data["fim"]

# Método de Gauss-Seidel
def Gauss_Seidel(A, b, x0, Eppara, maxit):

    ne = len(b)
    x = x0.copy()
    iter = 0
    Epest = np.full(ne, 100.0)

    while max(Epest) >= Eppara and iter < maxit:

        x_old = x.copy()

        for i in range(ne):

            sum1 = np.dot(A[i,:i], x[:i])
            sum2 = np.dot(A[i,i+1:], x_old[i+1:])

            x[i] = (b[i] - sum1 - sum2) / A[i,i]

        # Erro percentual estimado
        Epest = np.abs(
            (x - x_old) / np.maximum(np.abs(x), 1e-30)
        ) * 100

        iter += 1

    if max(Epest) >= Eppara:
        raise RuntimeError("Gauss-Seidel não convergiu.")

    return x

#permeabilidade harmônica nas fases em x 
k_harm_x = (
    2*k[:, :-1] * k[:,1:] 
    / (k[:, :-1] + k[:, 1:])
)

#permeabilidade harmônica nas fases em y
k_harm_y = (
    2*k[:-1, :] * k[1:, :] 
    / (k[:-1, :] + k[1:, :])
)

#transmissibilidades 
Tx = k_harm_x * Ax / (mu * dx)
Ty = k_harm_y * Ay / (mu * dy)

#falha selante: transmissibilidade hulas nas faces da falha 
if falha_ativa:
    Tx[inicio_falha:fim_falha + 1, coluna_falha] = 0.0

# Inicio da simulação
inicio = time.time()

for passo in range(1, nt + 1):

    #zerar matriz e vetor a cada passo de tempo
    A.fill(0)
    b.fill(0)

    #pt 1: acumulação
    A[np.arange(N), np.arange(N)] = acumulacao
    b[:] = acumulacao * Pold

    #pt 2: Faces em x 
    for i in range(ny):
        for j in range(nx - 1):
            m1 = i * nx + j
            m2 = m1 + 1
            T = Tx[i, j]

            A[m1, m1] += T
            A[m2, m2] += T
            A[m1, m2] -= T
            A[m2, m1] -= T

    #pt 3: Faces em y
    for i in range(ny - 1):
        for j in range(nx):
            m1 = i * nx + j
            m2 = m1 + nx
            T = Ty[i, j]

            A[m1, m1] += T
            A[m2, m2] += T
            A[m1, m2] -= T
            A[m2, m1] -= T
           


    # Contorno oeste - pressão prescrita
    if contorno_oeste == "pressao":
        for i in range(ny):
            m_b = i * nx
            Tb = (k[i, 0] * Ax) / (mu * dx/2)  # dx/2 ate a parede
            A[m_b, m_b] += Tb
            b[m_b] += Tb * p_oeste

    # Contorno leste - pressão prescrita
    if contorno_leste == "pressao":
        for i in range(ny):
            m_b = i * nx + (nx - 1)
            Tb = (k[i, nx - 1] * Ax) / (mu * dx/2)  # dx/2 ate a parede
            A[m_b, m_b] += Tb
            b[m_b] += Tb * p_leste

   
    # Contorno oeste - fluxo prescrito
    if contorno_oeste == "fluxo":
        for i in range(ny):
            m_b = i * nx
            b[m_b] += q_oeste

   
    # Contorno leste - fluxo prescrito
    if contorno_leste == "fluxo":
        for i in range(ny):
            m_b = i * nx + (nx - 1)
            b[m_b] += q_leste


   
    # Poço produtor central - Modelo de Peaceman
    if produtor_ativo:
        i = well_data['i']
        j = well_data['j']
        m = i * nx + j

        kp = k[i, j]
        re = 0.208 * dx

        WI = (2 * np.pi * kp * h) / (mu * np.log(re / rw))

        A[m, m] += WI
        b[m] += WI * pw


   
    
    # Poço injetor - vazão prescrita
    if injetor_ativo:
        m_inj = i_inj * nx + j_inj

        # Permeabilidade do bloco do injetor
        k_inj = k[i_inj, j_inj]

        # Raio equivalente de Peaceman
        re_inj = 0.208 * dx

        # Índice de produtividade do injetor
        WI_inj = (
            2 * np.pi * k_inj * h
            / (mu * np.log(re_inj / rw))
        )

        # Termo-fonte da vazão de injeção
        b[m_inj] += q_inj


    # Resolução do sistema linear
    if solver == "direto":
        Pnew = np.linalg.solve(A, b)
    else:
        Pnew = Gauss_Seidel(A, b, Pold, tolerancia, max_iter)

    # Vazão de produção pelo modelo de Peaceman
    if produtor_ativo:
        q_res = WI * (Pnew[m] - pw)  # m³/s no reservatório
    else:
        q_res = 0
    q_std = q_res / Bo   # m³ std/s

    # Pressão de fundo do injetor pelo modelo de Peaceman
    if injetor_ativo:
        pwf_inj = Pnew[m_inj] + q_inj / WI_inj
        pwf_injetor.append(pwf_inj / 1e5)  # Pa -> bar
    else:
        pwf_injetor.append(np.nan)

    # Produção acumulada
    Np += q_std * dt

    # Atualização do tempo
    tempo_atual += dt

    # Armazenamento dos resultados
    tempos.append(tempo_atual/3600)
    vazoes.append(q_std * 86400)
    producoes_acumuladas.append(Np)

    # Atualização das pressões
    Pold = Pnew.copy()
    P_new = Pnew.reshape(ny, nx).copy()

     # Salvar mapas de pressão nos tempos escolhidos
    for hora in tempos_snapshot:
        if passo == round(hora * 3600 / dt):
            pressoes_snapshot.append(P_new.copy())
            horas_snapshot.append(hora)

    P_old = P_new.copy()

fim = time.time()


def calcular_velocidades(P_new):

    # Velocidade de Darcy nas faces em x
    F_x = Tx * (P_new[:, :-1] - P_new[:, 1:])
    ux = F_x / Ax

    # Velocidade de Darcy nas faces em y
    F_y = Ty * (P_new[:-1, :] - P_new[1:, :])
    uy = F_y / Ay

    # Velocidade de Darcy no centro dos blocos em x
    ux_c = np.zeros((ny, nx))

    for j in range(nx):
        oeste = ux[:, j-1] if j > 0 else 0
        leste = ux[:, j] if j < nx - 1 else 0
        ux_c[:, j] = (oeste + leste) / 2

    # Velocidade de Darcy no centro dos blocos em y
    uy_c = np.zeros((ny, nx))

    for i in range(ny):
        norte = uy[i-1, :] if i > 0 else 0
        sul = uy[i, :] if i < ny - 1 else 0
        uy_c[i, :] = (norte + sul) / 2

    # Módulo da velocidade
    u_mod = np.sqrt(ux_c**2 + uy_c**2)

    return ux_c, uy_c, u_mod



# Velocidade no tempo final
ux_c, uy_c, u_mod = calcular_velocidades(P_new)

# Velocidade nos tempos dos snapshots
velocidades_snapshot = []
for P in pressoes_snapshot:
    velocidades_snapshot.append(calcular_velocidades(P))

print(f"Tempo de execução: {fim - inicio:.2f} segundos")
print(f"Pressão final no poço: {Pnew[m]/1e5:.2f} bar")
print(f"Vazão final: {vazoes[-1]:.4f} m³ std/dia")
print(f"Produção acumulada: {Np:.4f} m³ std")
print(f"Snapshots salvos: {horas_snapshot}")

if injetor_ativo:
    print(
        f"Pressão de fundo final do injetor: "
        f"{pwf_injetor[-1]:.2f} bar"
    )