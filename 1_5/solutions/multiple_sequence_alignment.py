import importlib
import itertools
from Bio import pairwise2
import numpy as np

module = importlib.import_module('1_2.solutions.needleman_wunsch')
nw = module.needleman_wunsch


def get_distance_matrix(sequences, matrix, gap_penalty):
    distance_matrix = [[None for _ in range(len(sequences))] for _ in range(len(sequences))]
    for first_ind, second_ind in itertools.combinations(list(range(len(sequences))), 2):
        score = pairwise2.align.globalds(sequences[first_ind], sequences[second_ind], matrix, gap_penalty, gap_penalty,
                                         score_only=True)
        distance_matrix[first_ind][second_ind] = distance_matrix[second_ind][first_ind] = score

    return distance_matrix


def get_argmax(distance_matrix):
    distance_matrix = np.array(distance_matrix)
    ind = np.unravel_index(np.argmax(distance_matrix, axis=None), distance_matrix.shape)
    return ind


def get_pair_alignment():
    pass


def main():
    gap_penalty = -2
    match = 2
    mismatch = -1

    matrix = {}
    alphabet = ["A", "G", "T", "C"]
    for pair in itertools.product(alphabet, repeat=2):
        if pair[0] == pair[1]:
            matrix[pair] = match
        else:
            matrix[pair] = mismatch


if __name__ == '__main__':
    main()
