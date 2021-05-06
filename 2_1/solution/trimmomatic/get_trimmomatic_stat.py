def main():
    number_of_deleted_seqs = 0
    number_of_unchanged_seqs = 0
    number_of_trimmed_seqs = 0
    with open("trimmomatic_log.txt", "r") as log_file, \
            open(r"..\original_seq_with_mistake_positions.txt", "r") as ground_truth:
        log = log_file.read().splitlines()
        ground_truth = ground_truth.read().splitlines()

    del ground_truth[0::2]

    for i in range(len(ground_truth)):
        id_, length_after_trim, start, end, length_of_cutoff = map(int, log[i].split())
        if length_after_trim == len(ground_truth[i]):
            number_of_unchanged_seqs += 1
        elif length_after_trim == 0:
            number_of_deleted_seqs += 1
        else:
            number_of_trimmed_seqs += 1

    print(f"Number of trimmed sequences: {number_of_trimmed_seqs}")
    print(f"Number of deleted sequences: {number_of_deleted_seqs}")
    print(f"Number of unchanged sequences: {number_of_unchanged_seqs}")


if __name__ == '__main__':
    main()
