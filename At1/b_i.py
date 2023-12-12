import numpy as np
import matplotlib.pyplot as plt

# Parâmetros iniciais
N = 500000
t0, dt, tMax = 0, 0.002, 365
S0, I0, R0 = 0.99 * N, 0.01 * N, 0
listS, listI, listR = [S0], [I0], [R0]


# define a função do modelo SIR
def SIR(S, I, R, beta, gamma):
    ds_dt = -beta * I * S / N
    di_dt = beta * I * S / N - (gamma * I)
    dr_dt = gamma * I
    return ds_dt, di_dt, dr_dt


# define a interpolação em Runge-Kutta de 2nd ordem
def RK2(SIR, S0, I0, R0, beta, gamma):
    k1s = SIR(S0, I0, R0, beta, gamma)[0] * dt
    k2s = SIR(S0 + k1s, I0, R0, beta, gamma)[0] * dt

    k1i = SIR(S0, I0, R0, beta, gamma)[1] * dt
    k2i = SIR(S0, I0 + k1i, R0, beta, gamma)[1] * dt

    k1r = SIR(S0, I0, R0, beta, gamma)[2] * dt
    k2r = SIR(S0, I0, R0 + k1r, beta, gamma)[2] * dt

    Sn = S0 + (1 / 2) * (k1s + k2s)
    In = I0 + (1 / 2) * (k1i + k2i)
    Rn = R0 + (1 / 2) * (k1r + k2r)

    return Sn, In, Rn


# Cria pontos conforme o avanço do tempo, sendo 16 o dMax para os parâmetros iniciais beta 0.5 e gamma 0.2
while t0 < 16:
    Snew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[0]
    Inew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[1]
    Rnew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[2]

    listS.append(Snew)
    listI.append(Inew)
    listR.append(Rnew)

    S0 = Snew
    I0 = Inew
    R0 = Rnew

    t0 += dt

# gamma de 0.2 p/ 0.1 quando t0 > dMax
while 16 < t0 <= tMax:
    Snew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.1)[0]
    Inew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.1)[1]
    Rnew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.1)[2]

    listS.append(Snew)
    listI.append(Inew)
    listR.append(Rnew)

    S0 = Snew
    I0 = Inew
    R0 = Rnew

    t0 += dt

# Entrega o dia de infecção máxima
maxI = max(listI)
indexMaxI = listI.index(maxI)
maxR = max(listR)
indexMaxR = listR.index(maxR)
print(f"Número máximo de infectados: {int(maxI)}")
print(f"Número máximo de Recuperados: {int(maxR)}")
tempo = np.linspace(t0-tMax, tMax, int(tMax/dt)+1)
print(len(tempo))
print(len(listS))
print(f"Dia de máxima infecção: {int(tempo[indexMaxI])}")
print(f"Dia final: {int(tempo[indexMaxR])}")

# configura o gráfico
plt.scatter(tempo, listS, label='Suscetíveis', color="green", s=1)
plt.scatter(tempo, listI, label='Infectados', color="orange", s=1)
plt.scatter(tempo, listR, label='Recuperados', color="blue", s=1)
plt.scatter(tempo[indexMaxI], maxI, label='Máxima Infecção', color="red", s=5)
plt.axvline(x=tempo[indexMaxR], ymin=0, ymax=maxR, linestyle=':', linewidth=1)
plt.legend()
plt.title("SIR Model", size=18)
plt.xlabel("Tempo", size=14)
plt.ylabel("População", size=14)
plt.show()
