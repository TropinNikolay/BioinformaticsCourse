import networkx as nx
from Bio import SeqIO
from graphviz import Digraph
import os

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


def delete_bubles(G):
    pass


def delete_tails(G):
    pass


def save_graph_stats(G):
    pass


def vizualize(G, filename):
    dot = Digraph()
    for v in G.nodes:
        dot.node(str(v))
        for edge_info in G.out_edges(v, data=True):
            dot.edge(str(v), str(edge_info[1]), label=str(edge_info[2]["substring"]))
    dot.render(filename)


def test_assembly(G):
    assemble = "".join([str(edge[2]["substring"]) for edge in list(G.edges(data=True))])
    with open("../data_input/DNA.txt", "r") as dna:
        dna = dna.readline()
    return "Success!" if dna == assemble else "Assembled incorrectly"


def main():
    # k = int(sys.argv[0])

    k = 50
    edges = count_edge_coverage(k)

    graph = create_graph(edges, k)
    vizualize(graph, "../data_output/graph")
    # save_graph_stats(graph)

    graph = compress_graph(graph)
    vizualize(graph, "../data_output/compressed_graph")
    # save_graph_stats(graph)

    print(test_assembly(graph))


if __name__ == '__main__':
    main()
