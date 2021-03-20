from Bio import Phylo


def transform_set_to_string(node):
    """
    Useful function for debugging
    """
    node.name = str(node.name)
    if node.is_terminal():
        return

    up_node, down_node = node.clades
    transform_set_to_string(up_node)
    if down_node:
        transform_set_to_string(down_node)


def forward(node):
    if node.is_terminal():
        node.name = set(node.name)
        return

    up_node, down_node = node.clades
    forward(up_node)
    node.name = up_node.name
    if down_node:
        forward(down_node)
        intersect = node.name.intersection(down_node.name)
        union = node.name.union(down_node.name)
        node.name = intersect if intersect else union


def backward(root):
    root.name = root.name.pop()
    if root.is_terminal():
        return

    up_node, down_node = root.clades
    if len(up_node.name) > 1:
        up_node.name = up_node.name.intersection(root.name)

    if down_node:
        if len(down_node.name) > 1:
            down_node.name = down_node.name.intersection(root.name)

    backward(up_node)
    backward(down_node)


def fitch(tree):
    forward(tree.clade)
    # transform_set_to_string(tree.clade)
    backward(tree.clade)
    return tree


if __name__ == '__main__':
    tree = Phylo.read("newick_tree_example.txt", format="newick")
    tree = fitch(tree)
    Phylo.draw(tree)
