import networkx as nx
from Bio import SeqIO
from copy import deepcopy
import difflib


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


def create_graph(edges, kmer_length):
    G = nx.DiGraph()
    while edges:
        edge, coverage = edges.popitem()
        G.add_edge(edge[:kmer_length - 1], edge[-kmer_length + 1:], coverage=coverage, substring=edge)
    return G


def compress_graph(G):
    nodes = deepcopy(list(G.nodes))
    for node in nodes:
        in_edges = list(G.in_edges(node, data=True))
        out_edges = list(G.out_edges(node, data=True))
        if len(in_edges) == 1 and len(out_edges) == 1:
            predecessor = in_edges[0][0]
            predecessor_edge_coverage = in_edges[0][2]["coverage"]
            predecessor_edge_substr = in_edges[0][2]["substring"]

            successor = out_edges[0][1]
            successor_edge_coverage = out_edges[0][2]["coverage"]
            successor_edge_substr = out_edges[0][2]["substring"]

            G.remove_node(node)

            coverage = (predecessor_edge_coverage * len(predecessor_edge_substr) + successor_edge_coverage * len(
                successor_edge_substr)) / (len(predecessor_edge_substr) + len(successor_edge_substr))

            match = difflib.SequenceMatcher(None, predecessor_edge_substr, successor_edge_substr)
            pred_pos, succ_pos, length = match.find_longest_match(0, len(predecessor_edge_substr), 0,
                                                                  len(successor_edge_substr))
            substring = predecessor_edge_substr + successor_edge_substr[length + 1:]

            G.add_edge(predecessor, successor, coverage=coverage, substring=substring)
    return G


def delete_bubles(G):
    pass


def delete_tails(G):
    pass


def save_graph_stats(G):
    pass


def main():
    # k = int(sys.argv[0])

    k = 50
    edges = count_edge_coverage(k)
    graph = create_graph(edges, k)
    save_graph_stats(graph)

    graph = compress_graph(graph)
    save_graph_stats(graph)


if __name__ == '__main__':
    main()
