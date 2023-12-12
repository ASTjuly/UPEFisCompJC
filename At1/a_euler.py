import matplotlib.pyplot as plt
import numpy as np

N = 500000
s0, i0, r0 = 0.99*N, 0.01*N, 0
t0, tMax, dt = 0, 365, 0.002
beta, gamma = 0.5, 0.2
sList, iList, rList = [s0], [i0], [r0]


def sir(s, i, r, beta, gamma):
    ds_dt = -beta * s/N * i
    di_dt = beta * s/N * i - gamma * i
    dr_dt = gamma * i
    return ds_dt, di_dt, dr_dt


while t0 < tMax:
    ds, di, dr = sir(s0, i0, r0, beta, gamma)

    sNew = s0 + ds * dt
    iNew = i0 + di * dt
    rNew = r0 + dr * dt

    sList.append(sNew)
    iList.append(iNew)
    rList.append(rNew)

    s0 = sNew
    i0 = iNew
    r0 = rNew
    t0 += dt

# Entrega o dia de infecção máxima
maxI = max(iList)
indexMaxI = iList.index(maxI)
maxR = max(rList)
indexMaxR = rList.index(maxR)
print(f"Número máximo de infectados: {int(maxI)}")
print(f"Número máximo de Recuperados: {int(maxR)}")
tempo = np.linspace(t0-tMax, tMax, int(tMax/dt)+1)
print(len(tempo))
print(len(sList))
print(f"Dia de máxima infecção: {int(tempo[indexMaxI])}")
print(f"Dia final: {int(tempo[indexMaxR])}")

plt.plot(tempo, sList, label="Suscetíveis")
plt.plot(tempo, iList, label="Infectados")
plt.plot(tempo, rList, label="Removidos")
plt.scatter(tempo[indexMaxI], maxI, label='Máxima Infecção', color="red", s=5)
plt.axvline(x=tempo[indexMaxR], ymin=0, ymax=maxR, linestyle=':', linewidth=1)
plt.legend()
plt.show()
