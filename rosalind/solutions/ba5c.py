def lcs(first_seq, second_seq):
    m = len(first_seq)
    n = len(second_seq)
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 or j == 0:
                dp[i][j] = 0
            elif first_seq[i - 1] == second_seq[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])
    ind = dp[m][n]
    lcs = [""] * (ind + 1)
    lcs[ind] = ""
    while m > 0 and n > 0:
        if first_seq[m - 1] == second_seq[n - 1]:
            lcs[ind - 1] = first_seq[m - 1]
            m -= 1
            n -= 1
            ind -= 1
        elif dp[m - 1][n] > dp[m][n - 1]:
            m -= 1
        else:
            n -= 1
    return lcs


with open("../data/rosalind_ba5c.txt") as data:
    first_seq = data.readline().strip()
    second_seq = data.readline().strip()
print(*lcs(first_seq, second_seq), sep="")
