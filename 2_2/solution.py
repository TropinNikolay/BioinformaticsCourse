"""
Input format: AAA AAT ATT TTA TTT

Output: TTTAAATTA

Log файл содержит последовательные преобразования графа,
который хранится в виде матрицы смежности
"""
import pandas as pd
import numpy as np
import logging


def calculate_common_length(first_node, second_node):
    p = 1
    max_common_length = 0
    while p <= len(first_node) and p <= len(second_node):
        if first_node[-p:] == second_node[:p]:
            max_common_length = p
        p += 1
    return max_common_length if max_common_length != 0 else np.NINF


def get_distance(ind_1, ind_2, columns):
    first_node, second_node = columns[ind_1], columns[ind_2]
    if first_node == second_node:
        return len(first_node)
    return calculate_common_length(first_node, second_node)


def create_distance_matrix(columns):
    number_of_contigs = len(columns)
    matrix = [[np.NINF for _ in range(number_of_contigs)] for _ in range(number_of_contigs)]
    for i in range(number_of_contigs - 1):
        for j in range(i + 1, number_of_contigs):
            matrix[i][j] = get_distance(i, j, columns)
            matrix[j][i] = get_distance(j, i, columns)
    return matrix


def find_max(matrix):
    matrix = matrix.to_numpy()
    ind = np.unravel_index(np.argmax(matrix, axis=None), matrix.shape)
    return ind, matrix[ind]


def main():
    columns = input().split()
    rows = columns.copy()
    matrix = pd.DataFrame(np.array(create_distance_matrix(columns)), columns=columns, index=rows)

    logging.info(matrix)

    while True:
        ind, common_length = find_max(matrix)

        if common_length == np.NINF:
            break

        max_ind_row, max_ind_column = ind
        nodes = matrix.columns
        new_node = nodes[max_ind_row] + nodes[max_ind_column][int(common_length):]
        matrix.drop(index=matrix.index[[max_ind_row, max_ind_column]], columns=nodes[[max_ind_row, max_ind_column]],
                    inplace=True)
        row = []
        column = []
        for vertex in matrix.columns:
            row.append(calculate_common_length(new_node, vertex))
            column.append(calculate_common_length(vertex, new_node))
        row_to_append = pd.DataFrame([row], columns=matrix.columns, index=[new_node])
        matrix = matrix.append(row_to_append)
        column.append(np.NINF)
        matrix[new_node] = column
        logging.info(matrix)

    print("".join(matrix.columns))


if __name__ == '__main__':
    logging.basicConfig(filename='log.txt', format='%(levelname)s:\n%(message)s\n', level=logging.INFO)
    main()
