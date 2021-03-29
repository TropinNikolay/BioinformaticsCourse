def neighbor_join_matix(n, matrix):
    d_new = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                d_new[i][j] = (n - 2) * matrix[i][j] - sum(matrix[i]) - sum(matrix[j])
    return d_new


def get_min_element(n, matrix):
    min_value = matrix[0][0]
    for i in range(n):
        for j in range(i, n):
            if matrix[i][j] < min_value:
                min_value = matrix[i][j]
                position = [i, j]
    return position


def neighbor_joining(n, d, node_list=None, dct=None):
    if dct is None:
        dct = {}
    if node_list is None:
        node_list = list(range(n))
    if n == 2:
        if node_list[0] not in dct.keys():
            dct[node_list[0]] = {}
        if node_list[1] not in dct.keys():
            dct[node_list[1]] = {}
        dct[node_list[0]][node_list[1]] = d[0][1]
        dct[node_list[1]][node_list[0]] = d[0][1]
        return dct
    else:
        d_new = neighbor_join_matix(n, d)
        min_element = get_min_element(n, d_new)
        i = min_element[0]
        j = min_element[1]
        diff = (sum(d[i]) - sum(d[j])) / (n - 2)
        limb_i = (d[i][j] + diff) / 2
        limb_j = (d[i][j] - diff) / 2
        add_row = [(d[k][i] + d[k][j] - d[i][j]) / 2 for k in range(n)] + [0]
        d.append(add_row)
        for l in range(n):
            d[l].append(add_row[l])
        m = node_list[-1] + 1
        node_list.append(m)
        [pair.pop(max(i, j)) for pair in d]
        [pair.pop(min(i, j)) for pair in d]
        d.remove(d[max(i, j)])
        d.remove(d[min(i, j)])
        node_i = node_list[i]
        node_j = node_list[j]
        node_list.remove(node_i)
        node_list.remove(node_j)
        if node_i not in dct.keys():
            dct[node_i] = {}
        if node_j not in dct.keys():
            dct[node_j] = {}
        if m not in dct.keys():
            dct[m] = {}

        dct[node_i][m] = limb_i
        dct[node_j][m] = limb_j
        dct[m][node_i] = limb_i
        dct[m][node_j] = limb_j
        neighbor_joining(n - 1, d, node_list, dct)
        return dct


if __name__ == '__main__':
    with open('../data/rosalind_ba7e.txt') as data:
        n = int(data.readline())
        distance_matrix = []
        for i in range(n):
            distance_matrix.append([int(t) for t in data.readline().split()])

    tree = neighbor_joining(n, distance_matrix)
    sorted_t = {k: v for k, v in sorted(tree.items(), key=lambda item: item[0])}
    for key in sorted_t:
        for node in sorted_t[key]:
            print(f"{key}->{node}:{sorted_t[key][node]:.2f}")
