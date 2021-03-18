from Bio import Phylo


def forward(node):
    pass


def backward(node):
    pass


def fitch(tree):
    forward(tree.clade)
    backward(tree.clade)
    return tree


if __name__ == '__main__':
    tree = Phylo.read("newick_tree_example.txt", format="newick")
    tree = fitch(tree)
    Phylo.draw_ascii(tree)
