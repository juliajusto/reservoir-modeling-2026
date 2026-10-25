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
eta = 1 / (phi * mu * ct)
rx = dt / dx**2
ry = dt / dy**2

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

# Inicio da simulação
inicio = time.time()

for passo in range(1, nt + 1):

    #zerar matriz e vetor a cada passo de tempo
    A.fill(0)
    b.fill(0)

    # Parte central
    for i in range(1, ny-1):
        for m in range(i*nx+1, (i+1)*nx-1):

            # Índice da coluna
            j = m % nx

            # Permeabilidades equivalentes nas faces
            kE = 2*k[i,j]*k[i,j+1] / (k[i,j]+k[i,j+1])
            kW = 2*k[i,j]*k[i,j-1] / (k[i,j]+k[i,j-1])
            kN = 2*k[i,j]*k[i-1,j] / (k[i,j]+k[i-1,j])
            kS = 2*k[i,j]*k[i+1,j] / (k[i,j]+k[i+1,j])

            if falha_ativa and inicio_falha <= i <= fim_falha:
              if j == coluna_falha:
                kE = 0
              if j == coluna_falha + 1:
                kW = 0

            # Coeficientes
            Ae = -eta * rx * kE
            Aw = -eta * rx * kW
            An = -eta * ry * kN
            As = -eta * ry * kS
            Ap = 1 - Ae - Aw - An - As

            # Matriz
            A[m,m] = Ap
            A[m,m-1] = Aw
            A[m,m+1] = Ae
            A[m,m-nx] = An
            A[m,m+nx] = As

            b[m] = Pold[m]

    # Fronteira norte - fluxo nulo
    for m in range(1, nx-1):
        i = 0
        j = m % nx

        # Permeabilidades equivalentes
        kE = 2*k[i,j]*k[i,j+1] / (k[i,j]+k[i,j+1])
        kW = 2*k[i,j]*k[i,j-1] / (k[i,j]+k[i,j-1])
        kS = 2*k[i,j]*k[i+1,j] / (k[i,j]+k[i+1,j])

        # Coeficientes
        Ae = -eta * rx * kE
        Aw = -eta * rx * kW
        An = 0
        As = -eta * ry * kS
        Ap = 1 - Ae - Aw - An - As

        # Matriz
        A[m,m] = Ap
        A[m,m-1] = Aw
        A[m,m+1] = Ae
        A[m,m+nx] = As

        b[m] = Pold[m]

    # Fronteira sul - fluxo nulo
    for m in range((ny-1)*nx+1, N-1):
        i = ny-1
        j = m % nx

        # Permeabilidades equivalentes nas faces
        kE = 2*k[i,j]*k[i,j+1] / (k[i,j]+k[i,j+1])
        kW = 2*k[i,j]*k[i,j-1] / (k[i,j]+k[i,j-1])
        kN = 2*k[i,j]*k[i-1,j] / (k[i,j]+k[i-1,j])

        # Coeficientes
        Ae = -eta * rx * kE
        Aw = -eta * rx * kW
        An = -eta * ry * kN
        As = 0
        Ap = 1 - Ae - Aw - An - As

        # Matriz
        A[m,m] = Ap
        A[m,m-1] = Aw
        A[m,m+1] = Ae
        A[m,m-nx] = An

        b[m] = Pold[m]

    # Fronteira oeste - fluxo nulo
    for i in range(1, ny-1):
        m = i * nx
        j = 0

        # Permeabilidades equivalentes
        kE = 2*k[i,j]*k[i,j+1] / (k[i,j]+k[i,j+1])
        kN = 2*k[i,j]*k[i-1,j] / (k[i,j]+k[i-1,j])
        kS = 2*k[i,j]*k[i+1,j] / (k[i,j]+k[i+1,j])

        # Coeficientes
        Ae = -eta * rx * kE
        Aw = 0
        An = -eta * ry * kN
        As = -eta * ry * kS
        Ap = 1 - Ae - Aw - An - As

        # Matriz
        A[m,m] = Ap
        A[m,m+1] = Ae
        A[m,m-nx] = An
        A[m,m+nx] = As

        b[m] = Pold[m]

    # Fronteira leste - fluxo nulo
    for i in range(1, ny-1):
        m = i * nx + (nx-1)
        j = nx-1

        # Permeabilidades equivalentes
        kW = 2*k[i,j]*k[i,j-1] / (k[i,j]+k[i,j-1])
        kN = 2*k[i,j]*k[i-1,j] / (k[i,j]+k[i-1,j])
        kS = 2*k[i,j]*k[i+1,j] / (k[i,j]+k[i+1,j])

        # Coeficientes
        Ae = 0
        Aw = -eta * rx * kW
        An = -eta * ry * kN
        As = -eta * ry * kS
        Ap = 1 - Ae - Aw - An - As

        # Matriz
        A[m,m] = Ap
        A[m,m-1] = Aw
        A[m,m-nx] = An
        A[m,m+nx] = As

        b[m] = Pold[m]

    # Canto Noroeste
    i = 0
    j = 0
    m = 0

    # Permeabilidades equivalentes
    kE = 2*k[i,j]*k[i,j+1] / (k[i,j]+k[i,j+1])
    kS = 2*k[i,j]*k[i+1,j] / (k[i,j]+k[i+1,j])

    # Coeficientes
    Ae = -eta * rx * kE
    Aw = 0
    An = 0
    As = -eta * ry * kS
    Ap = 1 - Ae - Aw - An - As

    # Matriz
    A[m,m] = Ap
    A[m,m+1] = Ae
    A[m,m+nx] = As

    b[m] = Pold[m]

    # Canto Nordeste
    i = 0
    j = nx-1
    m = nx-1

    # Permeabilidades equivalentes
    kW = 2*k[i,j]*k[i,j-1] / (k[i,j]+k[i,j-1])
    kS = 2*k[i,j]*k[i+1,j] / (k[i,j]+k[i+1,j])

    # Coeficientes
    Ae = 0
    Aw = -eta * rx * kW
    An = 0
    As = -eta * ry * kS
    Ap = 1 - Ae - Aw - An - As

    # Matriz
    A[m,m] = Ap
    A[m,m-1] = Aw
    A[m,m+nx] = As

    b[m] = Pold[m]

    # Canto Sudoeste
    i = ny-1
    j = 0
    m = (ny-1)*nx

    # Permeabilidades equivalentes
    kE = 2*k[i,j]*k[i,j+1] / (k[i,j]+k[i,j+1])
    kN = 2*k[i,j]*k[i-1,j] / (k[i,j]+k[i-1,j])

    # Coeficientes
    Ae = -eta * rx * kE
    Aw = 0
    An = -eta * ry * kN
    As = 0
    Ap = 1 - Ae - Aw - An - As

    # Matriz
    A[m,m] = Ap
    A[m,m+1] = Ae
    A[m,m-nx] = An

    b[m] = Pold[m]

    # Canto Sudeste
    i = ny-1
    j = nx-1
    m = N-1

    # Permeabilidades equivalentes
    kW = 2*k[i,j]*k[i,j-1] / (k[i,j]+k[i,j-1])
    kN = 2*k[i,j]*k[i-1,j] / (k[i,j]+k[i-1,j])

    # Coeficientes
    Ae = 0
    Aw = -eta * rx * kW
    An = -eta * ry * kN
    As = 0
    Ap = 1 - Ae - Aw - An - As

    # Matriz
    A[m,m] = Ap
    A[m,m-1] = Aw
    A[m,m-nx] = An

    b[m] = Pold[m]

    # Contorno oeste - pressão prescrita
    if contorno_oeste == "pressao":
        for i in range(ny):
            m_b = i * nx
            Cb = 2 * eta * rx * k[i, 0]  # dx/2 ate a parede
            A[m_b, m_b] += Cb
            b[m_b] += Cb * p_oeste

    # Contorno leste - pressão prescrita
    if contorno_leste == "pressao":
        for i in range(ny):
            m_b = i * nx + (nx - 1)
            Cb = 2 * eta * rx * k[i, nx - 1]
            A[m_b, m_b] += Cb
            b[m_b] += Cb * p_leste

    # Contorno oeste - fluxo prescrito
    if contorno_oeste == "fluxo":
        for i in range(ny):
            m_b = i * nx
            b[m_b] += (
                dt * q_oeste
                / (phi * ct * dx * dy * h)
            )
    # Contorno leste - fluxo prescrito
    if contorno_leste == "fluxo":
        for i in range(ny):
            m_b = i * nx + (nx - 1)
            b[m_b] += (
                dt * q_leste
                / (phi * ct * dx * dy * h)
            )

    # Poço produtos Central - Modelo de Peaceman
    if produtor_ativo:
        i = well_data['i']
        j = well_data['j']
        m = i * nx + j

        # Permeabilidade
        kp = k[i,j]

        # Raio equivalente de Peaceman
        re = 0.208 * dx

        # indice de produtividade
        WI = (2 * np.pi * kp * h) / (mu * np.log(re/rw))

        # Coeficiente do poço na equação implícita
        Cw = dt * WI / (phi * ct * dx * dy * h)

        # Modificação na matriz e no vetor
        A[m,m] += Cw
        b[m] += Cw * pw

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
        b[m_inj] += (
            dt * q_inj
            / (phi * ct * dx * dy * h)
        )

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

# Campo de velocidades
def calcular_velocidades(P_new):

    # Velocidade de Darcy nas faces em x
    ux = np.zeros((ny, nx - 1))

    for i in range(ny):
        for j in range(nx - 1):
            k_face = 2 * k[i,j] * k[i,j+1] / (k[i,j] + k[i,j+1])

            ux[i,j] = -(k_face / mu) * (
                P_new[i,j+1] - P_new[i,j]
            ) / dx


    #velocidade de Darcy nas faces em y
    uy = np.zeros((ny - 1, nx))

    for i in range(ny - 1):
        for j in range(nx):
            k_face = 2 * k[i,j] * k[i+1,j] / (k[i,j] + k[i+1,j])

            uy[i,j] = -(k_face / mu) * (
                P_new[i+1,j] - P_new[i,j]
            ) / dy


    # Falha selante: velocidade nula nas faces da falha
    if falha_ativa:
        for i in range(inicio_falha, fim_falha + 1):
            ux[i, coluna_falha] = 0


    # Velocidade em x no centro dos blocos
    ux_c = np.zeros((ny, nx))

    for i in range(ny):
        for j in range(nx):
            oeste = ux[i, j-1] if j > 0 else 0
            leste = ux[i, j] if j < nx - 1 else 0

            ux_c[i, j] = (oeste + leste) / 2


    # Velocidade em y no centro dos blocos
    uy_c = np.zeros((ny, nx))

    for i in range(ny):
        for j in range(nx):
            baixo = uy[i-1, j] if i > 0 else 0
            cima = uy[i, j] if i < ny - 1 else 0

            uy_c[i, j] = (baixo + cima) / 2


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