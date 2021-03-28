import numpy as np

np.seterr(divide='ignore', invalid='ignore')


def get_hmm_params(x, sigma, p_i, states):
    sigma_ind = {val: i for i, val in enumerate(sigma)}
    states_ind = {val: i for i, val in enumerate(states)}

    transition_matrix = np.zeros((len(states), len(states)))
    emission_matrix = np.zeros((len(states), len(sigma)))
    for x_i, state, next_state in zip(x, p_i, p_i[1:]):
        transition_matrix[states_ind[state], states_ind[next_state]] += 1
        emission_matrix[states_ind[state], sigma_ind[x_i]] += 1

    emission_matrix[states_ind[p_i[-1]], sigma_ind[x[-1]]] += 1
    e_sum = emission_matrix.sum(axis=1).reshape(len(states), 1)
    emission_matrix = np.where(e_sum, emission_matrix / e_sum, 1 / len(sigma))

    t_sum = transition_matrix.sum(axis=1).reshape(len(states), 1)
    transition_matrix = np.where(t_sum, transition_matrix / t_sum, 1 / len(sigma))
    return transition_matrix, emission_matrix


with open('../data/rosalind_ba10h.txt') as data:
    lines = [line.strip() for line in data.readlines()]
x = lines[0]
sigma = lines[2].split()
p_i = lines[4]
states = lines[6].split()

transition, emission = get_hmm_params(x, sigma, p_i, states)
transition = np.around(transition, decimals=3)
emission = np.around(emission, decimals=3)

print('', *states, sep='\t')
for i, row in zip(states, transition):
    print(i, *row, sep='\t')
print('--------')
print('', *sigma, sep='\t')
for i, row in zip(states, emission):
    print(i, *row, sep='\t')
