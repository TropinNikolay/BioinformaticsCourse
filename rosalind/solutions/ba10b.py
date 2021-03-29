emission_matrix = {}
with open("../data/rosalind_ba10b.txt") as data:
    for ind, line in enumerate(data):
        line = line.strip()
        if ind == 0:
            sequence = line
        elif ind == 2:
            alphabet = line.split()
        elif ind == 4:
            hidden_path = line
        elif ind == 6:
            states = line.split()
        elif ind >= 9:
            probabilities = line.split()
            state = probabilities.pop(0)
            for j, char in enumerate(alphabet):
                emission_matrix[f"{state}{char}"] = float(probabilities[j])

probability = 1
for i, state in enumerate(hidden_path):
    probability *= emission_matrix[f"{state}{sequence[i]}"]
print(probability)
