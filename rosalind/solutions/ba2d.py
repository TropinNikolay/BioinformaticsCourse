import numpy as np

k_mer_transform = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1]}
matrix_transform = {(1, 0, 0, 0): 'A', (0, 1, 0, 0): 'C', (0, 0, 1, 0): 'G', (0, 0, 0, 1): 'T'}


def get_all_k_mers(seq, k):
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]


def greedy_motif_search(dna, k, t):
    all_dna = np.array([[np.array([k_mer_transform[i] for i in j]).T for j in get_all_k_mers(s, k)] for s in dna])
    best_motifs = all_dna[:, 0]
    for k_mer in all_dna[0]:
        motifs = np.array([k_mer])
        for i in range(1, t):
            profile = motifs.sum(0) / t
            ind_match = np.product((all_dna[i] * profile).sum(1), axis=1).argmax()
            motifs = np.append(motifs, [all_dna[i, ind_match]], axis=0)
        if motifs.sum(0).max(0).sum() > best_motifs.sum(0).max(0).sum():
            best_motifs = motifs.copy()
    return [''.join([matrix_transform[tuple(i)] for i in motif.T]) for motif in best_motifs]


with open('../data/rosalind_ba2d.txt', "r") as data:
    k, t = map(int, data.readline().split())
    dna = [line.strip() for line in data.readlines()]
result = greedy_motif_search(dna, k, t)
print('\n'.join(result))
