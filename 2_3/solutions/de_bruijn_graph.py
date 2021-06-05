import networkx as nx
from Bio import SeqIO
import sys
import matplotlib.pyplot as plt


def create_graph(edges, kmer_length):
    G = nx.DiGraph()
    while edges:
        edge, coverage = edges.popitem()
        G.add_edge(edge[:kmer_length - 1], edge[-kmer_length + 1:], coverage=coverage)
    return G


def count_edge_coverage(kmer_length):
    edges = {}
    for record in SeqIO.parse("READS.fastq", "fastq"):
        read = record.seq
        for i in range(len(read) - kmer_length + 1):
            edge = read[i: i + kmer_length]
            if edge in edges:
                edges[edge] += 1
            else:
                edges[edge] = 1
    return edges


def compress_graph(graph):
    return graph


def main():
    # k = int(sys.argv[0])

    k = 50
    edges = count_edge_coverage(k)
    graph = create_graph(edges, k)
    graph = compress_graph(graph)

    # print(edges)
    # print(len(edges))
    # print(len(list(edges.keys())[0]))
    # print(graph.nodes)
    # print(len(graph.nodes))  # edges number + 1
    # print(graph.edges.data())
    # nx.draw(graph)
    # plt.show()


if __name__ == '__main__':
    main()
