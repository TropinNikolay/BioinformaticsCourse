import distance
from Bio import SeqIO


def levenshtein_dist(seq_1, seq_2):
    if len(seq_1) > len(seq_2):
        seq_1, seq_2 = seq_2, seq_1

    len_1 = len(seq_1)
    len_2 = len(seq_2)
    current_row = [i for i in range(len_1 + 1)]
    for i in range(1, len_2 + 1):
        previous_row = current_row
        current_row = [0] * (len_1 + 1)
        current_row[0] = i
        for j in range(1, len_1 + 1):
            if seq_1[j - 1] != seq_2[i - 1]:
                previous_row[j - 1] += 1
            current_row[j] = min(previous_row[j] + 1, current_row[j - 1] + 1, previous_row[j - 1])

    return current_row[len_1]


if __name__ == '__main__':
    first_seq, second_sec = [seq.seq for seq in SeqIO.parse('../data/f8.fasta', 'fasta')]
    print(levenshtein_dist(first_seq, second_sec))
    print(distance.levenshtein(first_seq, second_sec))
