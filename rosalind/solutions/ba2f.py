import random


def random_select_motifs(dna, k):
    motifs_lst = []
    for seq in dna:
        ind = random.randint(0, len(seq) - k)
        motifs_lst.append(seq[ind:ind + k])
    return motifs_lst


def profile(motif_lst, k):
    profile_lst = []
    for i in range(k):
        for j in range(len(motif_lst)):
            if j == 0:
                profile_lst.append({'A': 1, 'T': 1, 'C': 1, 'G': 1})
            profile_lst[i][motif_lst[j][i]] += 1
    return profile_lst


def motifs(profile_lst, dna):
    motifs_lst = []
    for Seq in dna:
        k = len(profile_lst)
        maxProbs = -1
        kmer = ''
        for i in range(len(Seq) - k + 1):
            Sum = 1
            for j in range(k):
                Sum *= (profile_lst[j][Seq[i + j]])
            if Sum > maxProbs:
                maxProbs = Sum
                kmer = Seq[i:i + k]
        motifs_lst.append(kmer)
    return motifs_lst


def score(motif_lst, k, t):
    profile_lst = profile(motif_lst, k)
    result = 0
    for a in range(len(profile_lst)):
        result += (4 + t - profile_lst[a][max(profile_lst[a], key=profile_lst[a].get)])
    return result


def randomized_motif_search(dna, k, t):
    motif_lst = random_select_motifs(dna, k)
    best_motifs = list(motif_lst)
    while True:
        profile_lst = profile(motif_lst, k)
        motif_lst = motifs(profile_lst, dna)
        if score(motif_lst, k, t) < score(best_motifs, k, t):
            best_motifs = list(motif_lst)
        else:
            return best_motifs


with open('../data/rosalind_ba2f.txt', "r") as data:
    k, t = map(int, data.readline().split())
    dna = [line.strip() for line in data.readlines()]

best_motifs = randomized_motif_search(dna, k, t)
for i in range(0, 1000):
    motif_lst = randomized_motif_search(dna, k, t)
    if score(motif_lst, k, t) < score(best_motifs, k, t):
        best_motifs = motif_lst
print(*best_motifs, sep="\n")
