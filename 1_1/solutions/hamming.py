from Bio import SeqIO


def hamming_dist(seq_1, seq_2):
    assert len(seq_1) == len(seq_2), "Расстояние Хэмминга не определено для строк различной длины"
    distance = 0
    for i in range(len(seq_1)):
        distance += seq_1[i] != seq_2[i]
    return distance


if __name__ == '__main__':
    first_seq, second_sec = [seq.seq for seq in SeqIO.parse('../data/gattaca.fasta', 'fasta')]
    print(hamming_dist(first_seq, second_sec))
