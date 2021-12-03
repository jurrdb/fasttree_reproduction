import newick

import main
from tree import Tree


def temp_generate_test_tree():
    # Test structure is same as example in https://en.wikipedia.org/wiki/Newick_format#Examples
    # But without the A node
    F = Tree()
    F.up_distance = 0
    F.name = 'F'
    B = Tree()
    B.up_distance = 0.2
    B.name = 'B'
    E = Tree()
    E.name = 'E'
    E.up_distance = 0.5
    F.left = B
    F.right = E
    C = Tree()
    C.up_distance = 0.3
    C.name = 'C'
    D = Tree()
    D.up_distance = 0.4
    D.name = 'D'
    E.left = C
    E.right = D
    return F


def traverse_tree_recursively(tree, newick, named_parent_nodes=True):
    parent_node_name = tree.name if named_parent_nodes else ""
    if not tree.left and not tree.right:
        return f"{tree.name}:{tree.up_distance}"
    elif tree.left and not tree.right:
        return f"(,{traverse_tree_recursively(tree.left, newick, named_parent_nodes)},){parent_node_name}:{tree.up_distance}"
    elif not tree.left and tree.right:
        return f"({traverse_tree_recursively(tree.right, newick, named_parent_nodes)}){parent_node_name}:{tree.up_distance}"
    else:
        return f"({traverse_tree_recursively(tree.left, newick, named_parent_nodes)},{traverse_tree_recursively(tree.right, newick, named_parent_nodes)}){parent_node_name}:{tree.up_distance}"


def to_newick(tree, named_parent_nodes=True):
    # # ----- DEBUG -----
    # test_tree = temp_generate_test_tree()
    # tree = test_tree
    # # Expected outcome: (B:0.2,(C:0.3,D:0.4)E:0.5)F;
    # # --- END DEBUG ---

    newick = ""
    newick = traverse_tree_recursively(tree, newick, named_parent_nodes)
    newick = f"{newick};"
    print("Created newick format:\n", newick)
    return newick


if __name__ == '__main__':
    # Debug driver function
    sequences, sequence_length = main.parse_input()
    tree = main.run(sequences, sequence_length)
    newick_string = to_newick(tree, named_parent_nodes=False)
    n_tree = newick.loads(newick_string)[0]
    print(n_tree.ascii_art())
