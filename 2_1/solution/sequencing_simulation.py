import random
import numpy as np

random.seed(42)

SEQ_LENGTH = 50000
MEAN_READ_LENGTH = 250
READ_STD = 30
LOWER_BOUNDARY = 1
UPPER_BOUNDARY = 80  # увеличил, т.к. при дефолтном значении ошибок в ридах почти нет
N_TIMES_PER_NUCLEOTIDE = 100
NUMBER_OF_READS = 50000
ALPHABET = "ATGC"


def make_single_read(seq_part):
    read_seq = ""
    quality_seq = ""
    seq_of_mistakes = ""

    for nucleotide in seq_part:
        mistake_number_of_times = np.random.randint(LOWER_BOUNDARY, UPPER_BOUNDARY + 1)
        nucleotide_pool = ALPHABET.replace(nucleotide, "")
        nucleotide_counter = {nucl: 0 for nucl in nucleotide_pool}
        nucleotide_counter[nucleotide] = N_TIMES_PER_NUCLEOTIDE - mistake_number_of_times

        for _ in range(mistake_number_of_times):
            nucleotide_readed = random.choice(nucleotide_pool)
            nucleotide_counter[nucleotide_readed] += 1

        max_count = -1
        max_nucl = None

        for nucl, count in nucleotide_counter.items():
            if count >= max_count:
                max_count = count
                max_nucl = nucl

        probability = 1 - max_count / N_TIMES_PER_NUCLEOTIDE
        quality = int(-10 * np.log10(probability))
        read_seq += max_nucl
        quality_seq += chr(quality + 33)

    assert len(seq_part) == len(read_seq)
    for i in range(len(seq_part)):
        seq_of_mistakes += "1" if read_seq[i] == seq_part[i] else "0"
    return read_seq, quality_seq, seq_of_mistakes


def make_some_reads(seq):
    with open("READS.fastq", "a") as fastq, \
            open("original_seq_with_mistake_positions.txt", "a") as data_out:
        for i in range(NUMBER_OF_READS):
            read_length = int(np.random.normal(MEAN_READ_LENGTH, READ_STD))
            start_ind = np.random.randint(0, SEQ_LENGTH - read_length)
            seq_part = seq[start_ind: start_ind + read_length]
            read, quality, mistakes = make_single_read(seq_part)
            fastq.write(f"@{i}\n{read}\n+\n{quality}\n")
            data_out.write(f"{seq_part}\n{mistakes}\n")


if __name__ == "__main__":
    sequence = ''.join(random.choices(ALPHABET, k=SEQ_LENGTH))
    make_some_reads(sequence)
