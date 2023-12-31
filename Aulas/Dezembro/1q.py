# Começa com um ponto único no centro

import numpy as np
import matplotlib.pyplot as plt


def rule_index(triplet):
    L, C, R = triplet
    index = 7 - (4 * L + 2 * C + R)
    return int(index)


def CA_run(initial_state, n_steps, rule_number):
    rule_string = np.binary_repr(rule_number, 8)
    rule = np.array([int(bit) for bit in rule_string])

    m_cells = len(initial_state)
    CA_run = np.zeros((n_steps, m_cells))
    CA_run[0, :] = initial_state
# laço for faz a leitura da vizinhança de cada ponto
    for step in range(1, n_steps):
        all_triplets = np.stack(
            [
                np.roll(CA_run[step-1, :], 1),
                CA_run[step - 1, :],
                np.roll(CA_run[step - 1, :], -1),
            ]
        )

        CA_run[step, :] = rule[np.apply_along_axis(rule_index, 0, all_triplets)]

    return CA_run


plt.rcParams["image.cmap"] = "binary"

# cria um ponto central isolado
IS = np.linspace(0, 0, 301)

# indice 150 é o ponto central
IS[150] = 1

"""rng = np.random.RandomState(20)"""
# o terceiro termo define a regra
data = CA_run(IS, 150, 30)

fig, ax = plt.subplots(figsize=(16, 9))
ax.matshow(data)
ax.axis(False)

plt.show()
