from Bio import SeqIO
from hamming import hamming_dist


def closest_substring(seq_1, seq_2):
    if len(seq_1) < len(seq_2):
        seq_1, seq_2 = seq_2, seq_1

    len_1 = len(seq_1)
    len_2 = len(seq_2)
    best_distance = float("+inf")
    for i in range(len_1 - len_2 + 1):
        dist = hamming_dist(seq_1[i: i + len_2], seq_2)
        if dist < best_distance:
            index = i
            best_distance = dist

    return index, seq_1[index: index + len_2], best_distance


if __name__ == '__main__':
    first_seq, second_sec = [seq.seq for seq in SeqIO.parse('data/f8.fasta', 'fasta')]
    print(closest_substring(first_seq, second_sec))
