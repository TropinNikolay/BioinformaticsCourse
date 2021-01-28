import random
import unittest

import Bio.Align.substitution_matrices as substitution_matrices
from Bio import SeqIO, pairwise2

from needleman_wunsch import needleman_wunsch

MATRIX = substitution_matrices.load("BLOSUM62")
GAP_PENALTY = -1
SEQ_LENGTH = 50


class TestAlignmentScore(unittest.TestCase):
    def setUp(self):
        self.seq_1_first, self.seq_1_second = [seq.seq for seq in SeqIO.parse('../../1_1/data/gattaca.fasta', 'fasta')]
        self.seq_2_first, self.seq_2_second = [seq.seq for seq in SeqIO.parse('../../1_1/data/gattaca2.fasta', 'fasta')]

    @staticmethod
    def get_score(fisrt_seq, second_seq, matrix, gap_penalty):
        assert len(fisrt_seq) == len(second_seq)
        score = 0
        for i in range(len(fisrt_seq)):
            if fisrt_seq[i] == "-" or second_seq[i] == "-":
                score += gap_penalty
            else:
                score += matrix[(fisrt_seq[i], second_seq[i])]
        return score

    def handle_sequences(self, fisrt_seq, second_seq):
        alignments = pairwise2.align.globalds(fisrt_seq, second_seq, MATRIX, GAP_PENALTY, GAP_PENALTY)
        true_score = alignments[0].score
        align_1, align_2, _ = needleman_wunsch(fisrt_seq, second_seq, MATRIX, gap_penalty=GAP_PENALTY)
        score = self.get_score(align_1, align_2, MATRIX, GAP_PENALTY)
        return true_score, score

    def test_score_from_data_1(self):
        true_score, score = self.handle_sequences(self.seq_1_first, self.seq_1_second)
        self.assertEqual(score, true_score)

    def test_score_from_data_2(self):
        true_score, score = self.handle_sequences(self.seq_2_first, self.seq_2_second)
        self.assertEqual(score, true_score)

    def test_random_seq(self):
        seq_1 = ''.join(random.choice(["A", "C", "G", "T"]) for _ in range(SEQ_LENGTH))
        seq_2 = ''.join(random.choice(["A", "C", "G", "T"]) for _ in range(SEQ_LENGTH))
        true_score, score = self.handle_sequences(seq_1, seq_2)
        self.assertEqual(score, true_score)


if __name__ == "__main__":
    unittest.main()
