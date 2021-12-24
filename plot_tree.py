from tree import Tree


def calculate_branch_lengths(topology):
    pass


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


def try_calculate_branch_length(tree, side_a, side_b):
    try:
        A = tree.__getattribute__(side_a).__getattribute__(side_a).__getattribute__(side_a)
        B = tree.__getattribute__(side_a).__getattribute__(side_a).__getattribute__(side_b)
        C = tree.__getattribute__(side_a).__getattribute__(side_b)
        D = tree.__getattribute__(side_b)
        # TODO fill in formulas to calculate branch length. HIER lekker verder volgende keer jongens
        # TODO zie "Branch lengths" section in paper
    except AttributeError:
        print("calculate_branch_length: skipping leaf node (" + side_a + ")")


def calculate_branch_length(tree):
    """
             R
            / \
           P   D
          / \
          N  C
         /\
        A  B
    """
    try_calculate_branch_length(tree, 'left', 'right')
    try_calculate_branch_length(tree, 'right', 'left')
    pass


def traverse_tree_recursively(tree, newick, named_parent_nodes=True):
    parent_node_name = tree.name if named_parent_nodes else ""
    if not tree.left and not tree.right:
        # Tree is leaf
        return f"{tree.name}:{tree.up_distance}"
    elif tree.left and not tree.right:
        # Tree has a left child, no right child
        return f"(,{traverse_tree_recursively(tree.left, newick, named_parent_nodes)},){parent_node_name}:{tree.up_distance}"
    elif not tree.left and tree.right:
        return f"({traverse_tree_recursively(tree.right, newick, named_parent_nodes)}){parent_node_name}:{tree.up_distance}"
    else:
        # Tree has two children
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
