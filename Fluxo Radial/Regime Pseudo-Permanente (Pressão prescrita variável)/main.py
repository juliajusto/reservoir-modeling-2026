import numpy as np
import matplotlib.pyplot as plt

mu = 3e-3         # Pa.s
ct = 130e-6 / 98066.5 # Pa^-1
k = 100 * 9.869233e-16 # m^2
phi = 0.20
h = 4             # m
B_o= 1.25          # m^3/m^3 std
qo = 35           # m^3 std/d
r_w = 0.10         # m
p0 = 200 * 1e5      # Pa
r_e = 500          # m
q_w = qo * B_o/ 86400  # m^3/s

eta = k/(phi*mu*ct)

print(f"eta = {eta:.4f} m²/s")
print(f"q_w = {q_w:.6e} m³/s")

def pressao_radial(r, t):
    coef = q_w * mu / (2 * np.pi * k * h)
    termo_tempo = 2*eta*t/r_e**2
    termo_geometrico = -np.log(r/r_w) + 1/2*(r/r_e)**2 + np.log(r_e/r_w) - 3/4
    p = p0 - coef * (termo_tempo + termo_geometrico)
    return p 


#Gráfico 1 - Perfis de pressão em diferentes tempos

r = np.linspace(r_w, r_e, 200)

tempos_perfis = [10, 30, 60, 100]

plt.figure(figsize=(10, 6))

for tempo_dias in tempos_perfis:
    tempo_segundos = tempo_dias * 24 * 3600
    
    P = pressao_radial(r, tempo_segundos)
    
    plt.plot(r, P/1e5, label=f'{tempo_dias} dias')

plt.xlabel('r (m)')
plt.ylabel('Pressão (bar)')
plt.title('Perfis de Pressão em Diferentes Tempos')
plt.legend()
plt.grid()

plt.savefig('perfis_pressao_radial.png', dpi=150)
plt.show()

#Gráfico 2 - Superfície 3D da pressão

t = np.linspace(10*24*3600, 100*24*3600, 100)

R, T = np.meshgrid(r, t)

P = pressao_radial(R, T)/1e5

fig = plt.figure(figsize=(10, 6))

ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(R, T/(24*3600), P,
                       rstride=2, cstride=1,
                       cmap=plt.cm.viridis,
                       linewidth=0.2,
                       alpha=1)

ax.set_xlabel('r (m)')
ax.set_ylabel('Tempo (dias)')
ax.set_zlabel('Pressão (bar)')

plt.title('Superfície de Pressão em Função do Raio e Tempo')

fig.colorbar(surf, shrink=0.5, aspect=10)

ax.view_init(30, 30)

plt.savefig('superficie_pressao_radial.png', dpi=150)
plt.show()