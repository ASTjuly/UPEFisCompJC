import numpy as np
import matplotlib.pyplot as plt

# Parâmetros iniciais
N = 500000
t0, dt, tMax = 0, 0.002, 365
S0, I0, R0 = 0.99 * N, 0.01 * N, 0
Sci0, Ici0, Rci0 = 0.99 * N, 0.01 * N, 0
Scii0, Icii0, Rcii0 = 0.99 * N, 0.01 * N, 0
listS, listI, listR = [S0], [I0], [R0]
listSci, listIci, listRci = [Sci0], [Ici0], [Rci0]
listScii, listIcii, listRcii = [Scii0], [Icii0], [Rcii0]


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

    Scinew = RK2(SIR, Sci0, Ici0, Rci0, beta=0.5, gamma=0.2)[0]
    Icinew = RK2(SIR, Sci0, Ici0, Rci0, beta=0.5, gamma=0.2)[1]
    Rcinew = RK2(SIR, Sci0, Ici0, Rci0, beta=0.5, gamma=0.2)[2]

    Sciinew = RK2(SIR, Scii0, Icii0, Rcii0, beta=0.5, gamma=0.2)[0]
    Iciinew = RK2(SIR, Scii0, Icii0, Rcii0, beta=0.5, gamma=0.2)[1]
    Rciinew = RK2(SIR, Scii0, Icii0, Rcii0, beta=0.5, gamma=0.2)[2]

    listS.append(Snew)
    listI.append(Inew)
    listR.append(Rnew)

    listSci.append(Scinew)
    listIci.append(Icinew)
    listRci.append(Rcinew)

    listScii.append(Sciinew)
    listIcii.append(Iciinew)
    listRcii.append(Rciinew)

    S0 = Snew
    I0 = Inew
    R0 = Rnew

    Sci0 = Scinew
    Ici0 = Icinew
    Rci0 = Rcinew

    Scii0 = Sciinew
    Icii0 = Iciinew
    Rcii0 = Rciinew

    t0 += dt

# beta variam c)i 0.25, c)ii 0.75
while 16 < t0 <= tMax:
    Snew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[0]
    Inew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[1]
    Rnew = RK2(SIR, S0, I0, R0, beta=0.5, gamma=0.2)[2]

    Scinew = RK2(SIR, Sci0, Ici0, Rci0, beta=0.25, gamma=0.2)[0]
    Icinew = RK2(SIR, Sci0, Ici0, Rci0, beta=0.25, gamma=0.2)[1]
    Rcinew = RK2(SIR, Sci0, Ici0, Rci0, beta=0.25, gamma=0.2)[2]

    Sciinew = RK2(SIR, Scii0, Icii0, Rcii0, beta=0.75, gamma=0.2)[0]
    Iciinew = RK2(SIR, Scii0, Icii0, Rcii0, beta=0.75, gamma=0.2)[1]
    Rciinew = RK2(SIR, Scii0, Icii0, Rcii0, beta=0.75, gamma=0.2)[2]

    listS.append(Snew)
    listI.append(Inew)
    listR.append(Rnew)

    listSci.append(Scinew)
    listIci.append(Icinew)
    listRci.append(Rcinew)

    listScii.append(Sciinew)
    listIcii.append(Iciinew)
    listRcii.append(Rciinew)

    S0 = Snew
    I0 = Inew
    R0 = Rnew

    Sci0 = Scinew
    Ici0 = Icinew
    Rci0 = Rcinew

    Scii0 = Sciinew
    Icii0 = Iciinew
    Rcii0 = Rciinew

    t0 += dt

# Entrega o dia de infecção máxima
maxI = max(listI)
indexMaxI = listI.index(maxI)

maxIci = max(listIci)
indexMaxIci = listIci.index(maxIci)

maxIcii = max(listIcii)
indexMaxIcii = listIcii.index(maxIcii)

maxR = max(listR)
indexMaxR = listR.index(maxR)

maxRci = max(listRci)
indexMaxRci = listRci.index(maxRci)

maxRcii = max(listRcii)
indexMaxRcii = listRcii.index(maxRcii)

print(f"Número máximo de infectados a): {int(maxI)}")
print(f"Número máximo de infectados c)i: {int(maxIci)}")
print(f"Número máximo de infectados c)ii: {int(maxIcii)}")
tempo = np.linspace(t0-tMax, tMax, int(tMax/dt)+1)
print(f"Dia de máxima infecção a): {int(tempo[indexMaxI])}")
print(f"Dia de máxima infecção c)i: {int(tempo[indexMaxIci])}")
print(f"Dia de máxima infecção c)ii: {int(tempo[indexMaxIcii])}")
print(f"Dia final a): {int(tempo[indexMaxR])}")
print(f"Dia final c)i: {int(tempo[indexMaxRci])}")
print(f"Dia final c)ii: {int(tempo[indexMaxRcii])}")

# configura o gráfico
plt.scatter(tempo, listI, label='I - a)', color="green", s=1)
plt.scatter(tempo, listIci, label='I - c)i', color="blue", s=1)
plt.scatter(tempo, listIcii, label='I - c)ii', color="orange", s=1)
plt.scatter(tempo[indexMaxI], maxI, label='Máxima Infecção a)', color="red", s=5)
plt.scatter(tempo[indexMaxIci], maxIci, label='Máxima Infecção c)i', color="purple", s=5)
plt.scatter(tempo[indexMaxIcii], maxIcii, label='Máxima Infecção c)ii', color="pink", s=5)
plt.axvline(x=tempo[indexMaxR], ymin=0, ymax=maxR, color="green", linestyle=':', linewidth=1)
plt.axvline(x=tempo[indexMaxRci], ymin=0, ymax=maxR, color="blue", linestyle=':', linewidth=1)
plt.axvline(x=tempo[indexMaxRcii], ymin=0, ymax=maxR, color="orange", linestyle=':', linewidth=1)
plt.legend()
plt.title("SIR Model", size=18)
plt.xlabel("Tempo", size=14)
plt.ylabel("População", size=14)
plt.show()
