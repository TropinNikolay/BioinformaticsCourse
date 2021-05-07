"""
Number of trimmed sequences: 26415
Number of deleted sequences: 23419
Number of unchanged sequences: 166
Percent of deleted nucleotides which were sequenced wrong: 31.39
Percent of deleted nucleotides which were sequenced properly: 89.46 (aka False positive rate)

Кажется, что trimmomatic удаляет большую часть правильно секвенированных нуклеотидов,
потому что удаляет много ридов.
"""


def main():
    number_of_deleted_seqs = 0
    number_of_unchanged_seqs = 0
    number_of_trimmed_seqs = 0
    num_all_deleted_nucl = 0  # кол-во всего удалённых нуклеотидов
    num_deleted_zero_nucl = (
        0  # кол-во всего удалённых нуклеотидов, которые были неверно секвенированы
    )
    num_deleted_one_nucl = (
        0  # кол-во всего удалённых нуклеотидов, которые были верно секвенированы
    )

    with open("trimmomatic_log.txt", "r") as log_file, open(
        r"..\original_seq_with_mistake_positions.txt", "r"
    ) as ground_truth:
        log = log_file.read().splitlines()
        ground_truth = ground_truth.read().splitlines()

    del ground_truth[0::2]

    for i in range(len(ground_truth)):
        id_, length_after_trim, start, end, length_of_cutoff = map(int, log[i].split())

        seq_length_ = len(ground_truth[i])
        right_ones_ = sum(int(x) for x in ground_truth[i][end:])
        left_ones_ = sum(int(x) for x in ground_truth[i][:start])
        all_ones_ = sum(int(x) for x in ground_truth[i])

        if length_after_trim == seq_length_:
            number_of_unchanged_seqs += 1
        elif length_after_trim == 0:
            number_of_deleted_seqs += 1
            num_all_deleted_nucl += seq_length_
            num_deleted_one_nucl += all_ones_
            num_deleted_zero_nucl += seq_length_ - all_ones_
        else:
            number_of_trimmed_seqs += 1
            num_all_deleted_nucl += seq_length_ - length_after_trim
            num_deleted_one_nucl += right_ones_ + left_ones_
            num_deleted_zero_nucl += seq_length_ - right_ones_ - left_ones_

    print(f"Number of trimmed sequences: {number_of_trimmed_seqs}")
    print(f"Number of deleted sequences: {number_of_deleted_seqs}")
    print(f"Number of unchanged sequences: {number_of_unchanged_seqs}")

    print(
        f"Percent of deleted nucleotides which were sequenced wrong: {round(num_deleted_zero_nucl / num_all_deleted_nucl * 100, 2)}"
    )
    print(
        f"Percent of deleted nucleotides which were sequenced properly: {round(num_deleted_one_nucl / num_all_deleted_nucl * 100, 2)}"
    )


if __name__ == "__main__":
    main()
