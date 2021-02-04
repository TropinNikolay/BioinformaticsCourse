import numpy as np
import itertools
import copy

deletion = -2
insertion = -2
match = 2
mismatch = -1

matrix = {}
alphabet = ["A", "G", "T", "C"]
for pair in itertools.product(alphabet, repeat=2):
    if pair[0] == pair[1]:
        matrix[pair] = match
    else:
        matrix[pair] = mismatch


def NWScore(x, y):
    rowPrev = np.zeros(len(y) + 1)
    for j, yj in enumerate(' ' + y):
        if j == 0:
            continue
        rowPrev[j] = rowPrev[j - 1] + insertion

    rowCur = np.zeros(len(y) + 1)
    lastLine = np.zeros(len(y) + 1)
    for i, xi in enumerate(' ' + x):
        if i == 0:
            continue
        rowCur[0] = rowPrev[0] + deletion
        for j, yj in enumerate(' ' + y):
            if j == 0:
                continue
            scoreSub = rowPrev[j - 1] + matrix[(xi, yj)]
            scoreDel = rowPrev[j] + deletion
            scoreIns = rowCur[j - 1] + insertion
            rowCur[j] = max([scoreIns, scoreDel, scoreSub])
        rowPrev = copy.deepcopy(rowCur)
    for j, _ in enumerate(' ' + y):
        lastLine[j] = rowCur[j]

    return lastLine


def Hirschberg(x, y, level=0):
    """
    This function calculates global alignment and
    prints out recursion tree in preorder traversal
    - every level of tree has additional indent.
    """
    print(f'{" " * level}{(x, y)}')
    z = ""
    w = ""
    if len(x) == 0:
        for _, yi in enumerate(y):
            z = z + '-'
            w = w + yi
    elif len(y) == 0:
        for _, xi in enumerate(x):
            z = z + xi
            w = w + '-'
    elif len(x) == 1:
        found = False
        for _, yi in enumerate(y):
            w = w + yi
            if yi == x and not found:
                z = z + x
                found = True
            else:
                z = z + '-'
        if not found:
            z = x + z[1:]
    elif len(y) == 1:
        found = False
        for _, xi in enumerate(x):
            z = z + xi
            if xi == y and not found:
                w = w + y
                found = True
            else:
                w = w + '-'
        if not found:
            w = y + w[1:]
    else:
        xmid = len(x) // 2
        scoreL = NWScore(x[:xmid], y)
        scoreR = NWScore(x[xmid:][::-1], y[::-1])
        revScoreR = scoreR[::-1]

        ymid = 0
        ymax = 0
        for i, _ in enumerate(scoreL):
            if scoreL[i] + revScoreR[i] > ymax:
                ymax = int(scoreL[i] + revScoreR[i])
                ymid = i
        z1, w1 = Hirschberg(x[:xmid], y[:ymid], level + 1)
        z2, w2 = Hirschberg(x[xmid:], y[ymid:], level + 1)
        z = z1 + z2
        w = w1 + w2
    return z, w


if __name__ == '__main__':
    print("Tree of recursion in preorder traversal:")
    z, w = Hirschberg("AGTACGCA", "TATGC")
    print("\nFinal alignment:")
    print(z)
    print(w)
