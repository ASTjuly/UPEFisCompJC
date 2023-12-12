import numpy as np
import matplotlib.pyplot as plt

# Parâmetros iniciais
N = 500000
t0, dt, tMax = 0, 0.002, 365
S0, I0, R0 = 0.99 * N, 0.01 * N, 0
Sbi0, Ibi0, Rbi0 = 0.99 * N, 0.01 * N, 0
Sbii0, Ibii0, Rbii0 = 0.99 * N, 0.01 * N, 0
listS, listI, listR = [S0], [I0], [R0]
listSbi, listIbi, listRbi = [Sbi0], [Ibi0], [Rbi0]
listSbii, listIbii, listRbii = [Sbii0], [Ibii0], [Rbii0]


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


while t0 < 16:
    Snew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[0]
    Inew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[1]
    Rnew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[2]

    Sbinew = RK2(SIR, Sbi0, Ibi0, Rbi0, beta=0.5, gamma=0.2)[0]
    Ibinew = RK2(SIR, Sbi0, Ibi0, Rbi0, beta=0.5, gamma=0.2)[1]
    Rbinew = RK2(SIR, Sbi0, Ibi0, Rbi0, beta=0.5, gamma=0.2)[2]

    Sbiinew = RK2(SIR, Sbii0, Ibii0, Rbii0, beta=0.5, gamma=0.2)[0]
    Ibiinew = RK2(SIR, Sbii0, Ibii0, Rbii0, beta=0.5, gamma=0.2)[1]
    Rbiinew = RK2(SIR, Sbii0, Ibii0, Rbii0, beta=0.5, gamma=0.2)[2]

    listS.append(Snew)
    listI.append(Inew)
    listR.append(Rnew)

    listSbi.append(Sbinew)
    listIbi.append(Ibinew)
    listRbi.append(Rbinew)

    listSbii.append(Sbiinew)
    listIbii.append(Ibiinew)
    listRbii.append(Rbiinew)

    S0 = Snew
    I0 = Inew
    R0 = Rnew

    Sbi0 = Sbinew
    Ibi0 = Ibinew
    Rbi0 = Rbinew

    Sbii0 = Sbiinew
    Ibii0 = Ibiinew
    Rbii0 = Rbiinew

    t0 += dt

# gammas alterados para b)i 0.1 e b)ii 0.4
while 16 < t0 <= tMax:
    Snew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[0]
    Inew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[1]
    Rnew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[2]

    Sbinew = RK2(SIR, Sbi0, Ibi0, Rbi0, beta=0.5, gamma=0.1)[0]
    Ibinew = RK2(SIR, Sbi0, Ibi0, Rbi0, beta=0.5, gamma=0.1)[1]
    Rbinew = RK2(SIR, Sbi0, Ibi0, Rbi0, beta=0.5, gamma=0.1)[2]

    Sbiinew = RK2(SIR, Sbii0, Ibii0, Rbii0, beta=0.5, gamma=0.4)[0]
    Ibiinew = RK2(SIR, Sbii0, Ibii0, Rbii0, beta=0.5, gamma=0.4)[1]
    Rbiinew = RK2(SIR, Sbii0, Ibii0, Rbii0, beta=0.5, gamma=0.4)[2]

    listS.append(Snew)
    listI.append(Inew)
    listR.append(Rnew)

    listSbi.append(Sbinew)
    listIbi.append(Ibinew)
    listRbi.append(Rbinew)

    listSbii.append(Sbiinew)
    listIbii.append(Ibiinew)
    listRbii.append(Rbiinew)

    S0 = Snew
    I0 = Inew
    R0 = Rnew

    Sbi0 = Sbinew
    Ibi0 = Ibinew
    Rbi0 = Rbinew

    Sbii0 = Sbiinew
    Ibii0 = Ibiinew
    Rbii0 = Rbiinew

    t0 += dt

# Entrega o dia de infecção máxima
maxI = max(listI)
indexMaxI = listI.index(maxI)

maxIbi = max(listIbi)
indexMaxIbi = listIbi.index(maxIbi)

maxIbii = max(listIbii)
indexMaxIbii = listIbii.index(maxIbii)

maxR = max(listR)
indexMaxR = listR.index(maxR)

maxRbi = max(listRbi)
indexMaxRbi = listRbi.index(maxRbi)

maxRbii = max(listRbii)
indexMaxRbii = listRbii.index(maxRbii)

print(f"Número máximo de infectados a): {int(maxI)}")
print(f"Número máximo de infectados b)i: {int(maxIbi)}")
print(f"Número máximo de infectados b)ii: {int(maxIbii)}")
tempo = np.linspace(t0-tMax, tMax, int(tMax/dt)+1)
print(f"Dia de máxima infecção a): {int(tempo[indexMaxI])}")
print(f"Dia de máxima infecção b)i: {int(tempo[indexMaxIbi])}")
print(f"Dia de máxima infecção b)ii: {int(tempo[indexMaxIbii])}")
print(f"Dia final a): {int(tempo[indexMaxR])}")
print(f"Dia final b)i: {int(tempo[indexMaxRbi])}")
print(f"Dia final b)ii: {int(tempo[indexMaxRbii])}")

# configura o gráfico
plt.scatter(tempo, listI, label='I - a)', color="green", s=1)
plt.scatter(tempo, listIbi, label='I - b)i', color="blue", s=1)
plt.scatter(tempo, listIbii, label='I - b)ii', color="orange", s=1)
plt.scatter(tempo[indexMaxI], maxI, label='Máxima Infecção a)', color="red", s=5)
plt.scatter(tempo[indexMaxIbi], maxIbi, label='Máxima Infecção b)i', color="purple", s=5)
plt.scatter(tempo[indexMaxIbii], maxIbii, label='Máxima Infecção b)ii', color="pink", s=5)
plt.axvline(x=tempo[indexMaxR], ymin=0, ymax=maxR, color="green", linestyle=':', linewidth=1)
plt.axvline(x=tempo[indexMaxRbi], ymin=0, ymax=maxR, color="blue", linestyle=':', linewidth=1)
plt.axvline(x=tempo[indexMaxRbii], ymin=0, ymax=maxR, color="orange", linestyle=':', linewidth=1)
plt.legend()
plt.title("SIR Model", size=18)
plt.xlabel("Tempo", size=14)
plt.ylabel("População", size=14)
plt.show()
