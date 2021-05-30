import random


def generate_DNA():
    alphabet = "ATGC"
    dna_length = 1000

    with open("DNA.txt", "w") as data:
        dna = ''.join(random.choices(alphabet, k=dna_length))
        data.write(dna)
    return dna


def generate_reads(dna):
    reads_length = 150
    quality = chr(126) * reads_length  # максимальное качество

    with open("READS.fastq", "a") as fastq:
        for i in range(len(dna) - reads_length + 1):
            read = dna[i: i + reads_length]
            fastq.write(f"@{i}\n{read}\n+\n{quality}\n")


def main():
    DNA = generate_DNA()
    generate_reads(DNA)


if __name__ == '__main__':
    main()
