import numpy as np


def viterbi(observations, states, start_proba, transition_matrix, emition_matrix):
    V = [{}]
    for state in states:
        V[0][state] = {"prob": start_proba[state] * emition_matrix[state][observations[0]], "prev": None}
    for t in range(1, len(observations)):
        V.append({})
        for state in states:
            max_proba = max(V[t - 1][prev_state]["prob"] * transition_matrix[prev_state][state] for prev_state in states)
            for prev_state in states:
                if V[t - 1][prev_state]["prob"] * transition_matrix[prev_state][state] == max_proba:
                    max_prob = max_proba * emition_matrix[state][observations[t]]
                    V[t][state] = {"prob": max_prob, "prev": prev_state}
                    break
    optimal = []
    max_prob = max(value["prob"] for value in V[-1].values())
    previous = None
    for state, data in V[-1].items():
        if data["prob"] == max_prob:
            optimal.append(state)
            previous = state
            break
    for t in range(len(V) - 2, -1, -1):
        optimal.insert(0, V[t + 1][previous]["prev"])
        previous = V[t + 1][previous]["prev"]
    return optimal, max_prob


with open("../data/rosalind_ba10c.txt") as data:
    x = data.readline().strip()
    data.readline()

    alphabet = data.readline().strip().split()
    alphabet = dict(zip(alphabet, range(len(alphabet))))
    observations = [alphabet[i] for i in x]
    data.readline()

    states = data.readline().strip().split()
    n = len(states)
    states = dict(zip(range(n), states))
    data.readline()
    data.readline()

    lines = []
    for _ in range(n):
        lines.append(data.readline().strip())
    transition = np.array([line.split()[1:] for line in lines], dtype=float)
    data.readline()
    data.readline()

    lines = []
    for _ in range(n):
        lines.append(data.readline().strip())
    emission = np.array([line.split()[1:] for line in lines], dtype=float)

path, p = viterbi(observations, states.keys(), [1 / n] * n, transition, emission)
print(''.join([states[i] for i in path]))
