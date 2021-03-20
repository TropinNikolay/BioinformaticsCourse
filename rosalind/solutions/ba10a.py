transition_matrix = {}
with open("../data/rosalind_ba10a.txt") as data:
    for ind, line in enumerate(data):
        line = line.strip()
        if ind == 0:
            sequence = line
        elif ind == 2:
            states = line.split()
            number_of_states = len(states)
        elif ind >= 5:
            probabilities = line.split()
            letter = probabilities.pop(0)
            for j, state in enumerate(states):
                transition_matrix[f"{letter}{state}"] = float(probabilities[j])

probability = 1 / len(states)
for i in range(len(sequence) - 1):
    probability *= transition_matrix[f"{sequence[i]}{sequence[i + 1]}"]
print(probability)
