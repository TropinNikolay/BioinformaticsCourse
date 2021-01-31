def needleman_wunsch(a, b, matrix, gap_penalty):
    n = len(a)
    m = len(b)
    score = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i

    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + matrix[(a[j - 1], b[i - 1])]
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    alignment_a = ""
    alignment_b = ""
    i = m
    j = n
    while i > 0 and j > 0:
        score_current = score[i][j]
        score_diagonal = score[i - 1][j - 1]
        score_up = score[i][j - 1]
        score_left = score[i - 1][j]

        if score_current == score_diagonal + matrix[(a[j - 1], b[i - 1])]:
            alignment_a += a[j - 1]
            alignment_b += b[i - 1]
            i -= 1
            j -= 1
        elif score_current == score_up + gap_penalty:
            alignment_a += a[j - 1]
            alignment_b += '-'
            j -= 1
        elif score_current == score_left + gap_penalty:
            alignment_a += '-'
            alignment_b += b[i - 1]
            i -= 1

    while j > 0:
        alignment_a += a[j - 1]
        alignment_b += '-'
        j -= 1
    while i > 0:
        alignment_a += '-'
        alignment_b += b[i - 1]
        i -= 1

    alignment_a = alignment_a[::-1]
    alignment_b = alignment_b[::-1]
    return alignment_a, alignment_b, score[m][n]
