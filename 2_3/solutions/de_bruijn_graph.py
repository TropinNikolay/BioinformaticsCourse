import networkx as nx
from Bio import SeqIO
from graphviz import Digraph
import os
import numpy as np
import matplotlib.pyplot as plt

os.environ["PATH"] += os.pathsep + r'C:\Program Files\Graphviz\bin'


def count_edge_coverage(kmer_length):
    edges = {}
    for record in SeqIO.parse("../data_input/READS.fastq", "fastq"):
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


def get_compression_nodes(G):
    nodes_to_compress = []
    for node in list(G.nodes):
        in_edges = list(G.in_edges(node, data=True))
        out_edges = list(G.out_edges(node, data=True))
        if len(in_edges) == 1 and len(out_edges) == 1:
            nodes_to_compress.append(node)
    return nodes_to_compress


def compress_graph(G):
    nodes_to_compress = get_compression_nodes(G)
    for node in nodes_to_compress:
        in_edge = list(G.in_edges(node, data=True))
        out_edge = list(G.out_edges(node, data=True))

        pred_node = in_edge[0][0]
        next_node = out_edge[0][1]

        pred_node_edge = in_edge[0][2]["substring"]
        next_node_edge = out_edge[0][2]["substring"]

        pred_node_coverage = in_edge[0][2]["coverage"]
        next_node_coverage = out_edge[0][2]["coverage"]

        coverage = (pred_node_coverage * len(pred_node_edge) + next_node_coverage * len(next_node_edge)) / (
                len(pred_node_edge) + len(next_node_edge))

        G.add_edge(pred_node, next_node, coverage=coverage, substring=pred_node_edge + next_node_edge[len(node):])

        G.remove_node(node)
    return G


def get_tail_nodes(G):
    tail_nodes = []
    for node in list(G.nodes):
        in_edges = list(G.in_edges(node))
        out_edges = list(G.out_edges(node))
        if len(in_edges) == 1 and len(out_edges) == 0:
            tail_nodes.append(node)
    return tail_nodes


def get_tail_distribution(G, tails):
    tails_distribution = []
    for node in tails:
        in_edge = list(G.in_edges(node, data=True))
        edge_length = len(in_edge[0][2]["substring"])
        edge_coverage = in_edge[0][2]["coverage"]
        tails_distribution.append((edge_length * edge_coverage, node))
    tails_distribution.sort()
    return tails_distribution


def delete_tails(G, *, threshold):
    tail_nodes = get_tail_nodes(G)
    tails_distribution = get_tail_distribution(G, tail_nodes)
    ind = np.ceil(threshold * len(tails_distribution)).astype(int)
    for _, node in tails_distribution[:ind]:
        G.remove_node(node)
    return G


def delete_bubles(G, *, buble_length):
    return G


def save_graph_stats(G, out_dir):
    with open(os.path.join(out_dir, "Stats.txt"), "w") as stat:
        stat.write(f"Number of nodes: {G.number_of_nodes()}\n")
        stat.write(f"Number of edges: {G.number_of_edges()}\n")

    coverage = []
    for from_v, to_v, data in G.edges(data=True):
        coverage.append(data["coverage"])

    plt.hist(coverage, range=[0, max(coverage) + 1])
    plt.title("Coverage distribution")
    plt.savefig(os.path.join(out_dir, "Coverage_distribution.png"))
    plt.close()

    in_degree = []
    for node, in_deg in G.in_degree:
        in_degree.append(in_deg)

    plt.hist(in_degree)
    plt.title("Degree in")
    plt.savefig(os.path.join(out_dir, "Degree_in.png"))
    plt.close()

    out_degree = []
    for node, out_deg in G.out_degree:
        out_degree.append(out_deg)

    plt.hist(out_degree)
    plt.title("Degree out")
    plt.savefig(os.path.join(out_dir, "Degree_out.png"))
    plt.close()


def vizualize(G, filename):
    dot = Digraph()
    for v in G.nodes:
        dot.node(str(v))
        for edge_info in G.out_edges(v, data=True):
            dot.edge(str(v), str(edge_info[1]), label=str(edge_info[2]["substring"]))
    dot.render(filename)


def test_assembly(G):
    assemble = "".join([str(edge[2]["substring"]) for edge in list(G.edges(data=True))])
    with open("../data_input/DNA/DNA.txt", "r") as dna:
        dna = dna.readline()
    return "Success!" if dna == assemble else "Assembled incorrectly"


def main():
    k = 50  # для работы с большой ДНК, для DNA_small k = 4
    edges = count_edge_coverage(k)

    graph = create_graph(edges, k)
    vizualize(graph, "../data_output/DNA/original/graph")
    save_graph_stats(graph, "../data_output/DNA/original")

    graph = compress_graph(graph)
    vizualize(graph, "../data_output/DNA/compressed/compressed_graph")
    save_graph_stats(graph, "../data_output/DNA/compressed")

    # print(test_assembly(graph))  # для проверки большой ДНК

    graph = delete_tails(graph, threshold=0.3)
    save_graph_stats(graph, "../data_output/DNA/without_tails")


if __name__ == '__main__':
    main()
