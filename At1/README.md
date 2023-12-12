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

Número máximo de infectados: 118794; Número máximo de Recuperados: 447146; Dia de máxima infecção: 16; Dia final: 213.


Assim, finalizando o item a).

### b) Considere que, por algum motivo, a taxa de infecção seja constante e igual a **&beta;=0.5** e a taxa de recuperação *dependa do tempo* das seguintes formas: <br /> **i) &gamma;=0.2, para d < d<sub>max</sub> e &gamma;=0.1, para d >= d<sub>max</sub>** <br /> **ii) &gamma;=0.2, para d < d<sub>max</sub> e &gamma;=0.4, para d >= d<sub>max</sub>** <br /> Plote em um mesmo gráfico as curvas para **I(t)** para os três casos acima, **a), b-i) e b-ii)**.

Para que fosse possível variar o &gamma; conforme o indicado na questão, é necessário a separação do laço *while* em dois: um para o momento antes do pico dos parâmetros iniciais, e outro após o pico dos parâmetros iniciais. Desta forma, o trecho do código fica:
```
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
```

Observe como o número 16 dita onde será quebrado o *while*, isto porque o valor 16 é respectivo ao dia de infecção máxima nos parâmetros iniciais da questão, e este valores mudam somente após este pico.

Desta forma, obtemos os seguintes gráficos:

Gráfico obtido a partir do **[programa do item b-i)](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/b_i.py)**
![UPEFisCompJC Preview](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/Graph/b_i.png)
Seus dados:

Número máximo de infectados: 155527; Número máximo de Recuperados: 491632; Dia de máxima infecção: 21; Dia final: 330.

É possível observar que a função criou um novo pico conforme a mudança da taxa de recuperação &gamma; para um valor menor. Isto é possível pelo papel da taxa de recuperação estar interligado com o tempo em que o indivíduo estará infectado desta maneira: $\frac{1}{\gamma}=$Dias ifectados, neste caso o valor passa de $\frac{1}{0.2}=5$ para $\frac{1}{0.1}=10$, dobrando o dia em que os Infectados terão que transitar para o número de Recuperados e, assim, criando um novo acúmulo de Infectados gerando um novo pico no dia 21. Além da consequência de um atraso no dia final do surto, que agora está no dia 330, quase ao final do tempo máximo.

Gráfico obtido a partir do **[programa do item b-ii)](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/b_ii.py)**
![UPEFisCompJC Preview](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/Graph/b_ii.png)
Seus dados:

Número máximo de infectados: 118422; Número máximo de Recuperados: 372193; Dia de máxima infecção: 16; Dia final: 122.

Nesse caso, é possível observar que a alteração da taxa de recuperação &gamma; para um valor maior, resultou na diminuição do caso de infectados e em uma aceleração para o fim do surto, já que os contaminados irão passar $\frac{1}{0.4}=2,5$ dias para saírem do estado de Infectados para Recuperados — tempo menor que o inicial de 5 dias.

Agora, o gráfico obtido a partir do **[programa do item b)](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/b.py)** pedido na questão que compara a função *I(t)* para as situações dos itens a), b-i) e b-ii):
![UPEFisCompJC Preview](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/Graph/b.png)

Dados em comparação:
|  | a) | b-i) | b-ii) |       
| -- | -- | -- | -- |
| d<sub>max</sub> | 16 | 21 | 16 |
| d<sub>final</sub> | 213 | 330 | 122 |

### c) Considere agora que a taxa de recuperação seja constante e igual a **&gamma;=0.2** e a taxa de transmissão dependa do tempo das seguintes formas: <br /> **i) &beta;=0.5, para d < d<sub>max</sub> e &beta;=0.25, para d >= d<sub>max</sub>** <br /> **ii) &beta;=0.5, para d < d<sub>max</sub> e &beta;=0.75, para d >= d<sub>max</sub>** <br /> Plote em um mesmo gráfico as curvas para **I(t)** para os três casos acima, **a), c-i) e c-ii)**.

O mesmo passo inicial do item b foi utilizado no trecho correspondente aos laços *while* do programa, sendo sua única diferença o parâmetro em que foram aplicadas as mudanças (exemplo do item c-i) abaixo):
```
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

# beta de 0.5 p/ 0.25 quando t0 > dMax
while 16 < t0 <= tMax:
    Snew = RK2(SIR, S0, I0, R0, beta=0.25, gamma=0.2)[0]
    Inew = RK2(SIR, S0, I0, R0, beta=0.25, gamma=0.2)[1]
    Rnew = RK2(SIR, S0, I0, R0, beta=0.25, gamma=0.2)[2]

    listS.append(Snew)
    listI.append(Inew)
    listR.append(Rnew)

    S0 = Snew
    I0 = Inew
    R0 = Rnew

    t0 += dt
```

Desta forma, obtemos os seguintes gráficos:

Gráfico obtido a partir do **[programa do item c-i)](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/c_i.py)**
![UPEFisCompJC Preview](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/Graph/c_i.png)
Seus dados:

Número máximo de infectados: 118422; Número máximo de Recuperados: 372168; Dia de máxima infecção: 16; Dia final: 224.

Já nesta situação, observa-se que o pico é mantido no dia 16, embora ocorra um atraso para o controle do surto, que chega apenas dia 224 — 11 dias depois em comparação caso os valores permancessem os mesmos, como no item a). Essas ocasionalidades são oriundas da diminuição da taxa de transmissão **&beta;**, que é responsável por ditar a possibilidade de um indivíduo Infectado transmitir para outro Suscetível. Conforme o número diminui, menos chance de infectar alguém que não contaminou-se com o vírus, embora, a doença ainda será lentamente espalhada até chegar no dia em que for controlada.

Gráfico obtido a partir do **[programa do item c-ii)](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/c_ii.py)**
![UPEFisCompJC Preview](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/Graph/c_ii.png)
Seus dados:

Número máximo de infectados: 135435; Número máximo de Recuperados: 479345; Dia de máxima infecção: 18; Dia final: 190.

Já nesse caso, a variação temporal no **&beta;** ocasionou em um novo pico 2 dias depois do primeiro e no adiantamento do dia de controle do surto. Isto ocorre porque o vírus possui uma taxa de transmissão maior, ou seja, uma doença altamente contagiosa e que dissemina com bastante verocidade. Embora, com a mesma taxa de recuperação **&gamma;=0.2**, os Infectados logo irão ser transferidos para os Recuperados, o que resultará no dia do controle chegar mais cedo que a situação do item a).

Agora, o gráfico obtido a partir do **[programa do item c)](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/c.py)** pedido na questão que compara a função *I(t)* para as situações dos itens a), c-i) e c-ii):
![UPEFisCompJC Preview](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/Graph/c.png)

Dados em comparação:
|  | a) | c-i) | c-ii) |
| -- | -- | -- | -- |
| d<sub>max</sub> | 16 | 16 | 18 |
| d<sub>final</sub> | 213 | 224 | 190 |

## Conclusões da Atividade

Por fim, é possível entender como o Método de Runge-Kutta contribui de forma mais precisa na interpolação dos dados de uma EDO. Em comparativo, refazendo a situação do item a) com o seguinte **[programa)](https://github.com/ASTjuly/UPEFisCompJC/blob/main/At1/a_euler.py)** utilizando o Método de Euler, os valores são ainda próximos porém com uma diferença pequena:

|  | I<sub>max</sub> | R<sub>max</sub> | d<sub>max</sub> | d<sub>final</sub> |
| -- | -- | -- | -- | -- |
| Euler | 118767 | 447065 | 16 | 213 |
| RK2 | 118794 | 447146 | 16 | 213 |

Visto a utilização do Método de Runge-Kutta de 2ª Ordem ser o nível mais básico e não oferecer tanta precisão, a diferença entre os métodos matemáticos são mínimas mas estão presentes. Conforme escalasse o número de ordem em Runge-Kutta, maior a exatidão na qual serão observados os dados.

Também foram estudados na atividade as alterações devidas as variações de **&beta;(t)** e **&gamma;(t)**. Que mudaram os gráficos em seus pontos de infecção máxima e momento de controle do surto epidemiológico conforme o previsto pelo modelo teórico SIR. Assim, comprovando a funcionalidade do modelo no método computacional.
