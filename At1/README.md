# Primeira Atividade - 04/12/2023
## Atividade sobre Método de Runge-Kutta e Modelo SIR

### Imagine que as taxas de infecção e de recuperação para a dinâmica de uma virose sobre a população de 500.000 habitantes de uma cidade sejam de **&beta; = 0.5** e **&gamma; = 0.2**. **Usando o método de Runge-Kutta de 2ª ordem** e considerando a unidade de tempo como um dia, encontre:<br /> <br />a) As curvas S(t), I(t) e R(t), usando as equações **(II)** para a dinâmica de infecção dessa cidade. Em qual dia o número de infectados será máximo? (Chame esse dia de d<sub>max</sub>) Em aproximadamente qual dia o surto da doença terá sido controlado?

### Resposta

Primeiramente, é preciso definir a função SIR no programa. Sendo esta:

```
# define a função do modelo SIR
def SIR(S, I, R, beta, gamma):
    ds_dt = -beta * I * S / N
    di_dt = beta * I * S / N - (gamma * I)
    dr_dt = gamma * I
    return ds_dt, di_dt, dr_dt
```

Logo após, é necessário a implementação do Método de Runge-Kutta de 2ª Ordem — que será chamado ao longo do trabalho de RK2. Para tal, usaremos o seguinte trecho que define a função RK2 no código:

```
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
```
Em que *Sn* equivale ao novo número de Suscetíveis obtido a partir da função SIR com parâmetro inicial *S0*, o mesmo equivale para o número de Infectados *In* e o número de Recuperados *Rn*, respectivamente calculados a partir de *I0* e *R0*.

Perceba que foram calculados diferentes valores para cada *k1* e *k2* dependendo de qual parâmetro estava em questão.

Logo após a definição de RK2, cria-se um laço *while* que irá definir os pontos de *S*, *I* e *R* conforme o avançar do tempo em variação de *dt*.

```
while t0 < tMax:
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
```

Agora, para obter o dia de máxima infecção *d<sub>max</sub>* foi necessário criar um local para armazenar o ponto de máximo da função que define *I*. Este local foi chamado de "maxI" e aparece no programa da seguinte forma:

```
# Entrega o dia de infecção máxima
maxI = max(listI)
indexMaxI = listI.index(maxI)

# Imprime o número máximo de infectados e o dia em que aconteceu o d_max
print(f"Número máximo de infectados a): {int(maxI)}")
print(f"Dia de máxima infecção a): {int(tempo[indexMaxI])}")

# Gera um ponto de máximo na função I(t) no gráfico
plt.scatter(tempo[indexMaxI], maxI, label='Máxima Infecção a)', color="red", s=5)
```

Já para o dia de controle do surto (chamado no código como "dia final"), a ideia foi um pouco mais elaborada, visto que o que define um controle no surto não está localizado somente na função *I(t)*, mas também no ponto máximo da função *R(t)*. De acordo com o modelo SIR utilizado na atividade:

$\frac{ds}{dt}=-\beta * I * \frac{S}{N}$<br />
$\frac{di}{dt}=\beta * I * \frac{S}{N} - \gamma * I$<br />
$\frac{dr}{dt}=\gamma * I$

É possível observar, por fim, que o valor de *I(t)* continuará diminuindo conforme o avanço do tempo. Já o valor de *R(t)*, por outro lado, terá um limite máximo melhor definido na computação, e este valor definirá o dia de controle do surto de infecção porque quando o valor de Recuperados se torna máximo, logo o valor de Infectados se torna insuficiente para infectar habitantes Suscetíveis e criar novos Recuperados.

Partindo disto, a parcela do código que irá determinar o dia final será:

```
# Entrega o dia de máxima dos recuperados
maxR = max(listR)
indexMaxR = listR.index(maxR)

# Imprime o dia final do surto a partir do momento em que R atinge o pico
print(f"Dia final: {int(tempo[indexMaxR])}")

# Cria uma linha vertical que demonstra o dia final no eixo temporal do gráfico
plt.axvline(x=tempo[indexMaxR], ymin=0, ymax=maxR, color="green", linestyle=':', linewidth=1)
```
De acordo com o **[programa do item a)](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/a.py)** e os parâmetros iniciais entregues no item, o gráfico que contem as curvas *S(t)*, *I(t), e *R(t)* é:

![UPEFisCompJC Preview](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/Graph/a.png)
Sendo seus dados:

Número máximo de infectados: 118794;

Número máximo de Recuperados: 447146;

Dia de máxima infecção: 16;

Dia final: 213.


Assim, finalizando o item a).
