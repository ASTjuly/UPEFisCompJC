import numpy as np
import matplotlib.pyplot as plt


def rule_index(quintuplet):
    # Ll = esquerda distante, L = esquerda, C = centro, R = direita, Rr= direita distante
    Ll, L, C, R, Rr = quintuplet
    index = 31 - (16 * Ll + 8 * L + 4 * C + 2 * R + Rr)
    return int(index)


def CA_run(initial_state, n_steps, rule_number):
    rule_string = np.binary_repr(rule_number, 32)
    rule = np.array([int(bit) for bit in rule_string])

    m_cells = len(initial_state)
    CA_run = np.zeros((n_steps, m_cells))
    CA_run[0, :] = initial_state

    for step in range(1, n_steps):
        all_quintuplets = np.stack(
            [
                np.roll(CA_run[step-1, :], 1),
                np.roll(CA_run[step - 1, :], 1),
                CA_run[step - 1, :],
                np.roll(CA_run[step - 1, :], -1),
                np.roll(CA_run[step - 1, :], -1),
            ]
        )

        CA_run[step, :] = rule[np.apply_along_axis(rule_index, 0, all_quintuplets)]

    return CA_run


plt.rcParams["image.cmap"] = "binary"

# cria um ponto central isolado
IS = np.linspace(0, 0, 301)

# indice 150 é o ponto central
IS[150] = 1

# rng = np.random.RandomState(0)
data = CA_run(IS, 150, 4246798)

fig, ax = plt.subplots(figsize=(16, 9))
ax.matshow(data)
ax.axis(False)

plt.show()
