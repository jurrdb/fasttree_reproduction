from distance_measures import profile_distance_corrected as d


def traverse_tree_recursively(tree, newick, named_parent_nodes=True, sequence_len=1):
    # output the tree in newick format by recursively traversing the tree and calculate the distances of branches
    parent_node_name = tree.name if named_parent_nodes else ""
    branch_length = 0
    if not tree.left and not tree.right:
        # Tree is leaf
        if tree.parent is not None and tree.parent.parent is not None:
            # Branch length leading to leaf A:
            # d(A,BC) = ( d(A,B) + d(A,C) - d(B,C) ) / 2
            N = tree.parent
            P = tree.parent.parent
            A = tree
            B = N.left if N.left != A else N.right
            C = P.right if P.right != N else P.left
            branch_length = round((d(A, B) + d(A, C) - d(B, C)) / (2 * sequence_len), 3)
        return f"{tree.name}:{branch_length}"
    elif tree.left and not tree.right:
        # Tree has a left child, but no right child
        return f"(,{traverse_tree_recursively(tree.left, newick, named_parent_nodes, sequence_len)},){parent_node_name}:{branch_length}"
    elif not tree.left and tree.right:
        # Tree has a right child, but no left child
        return f"({traverse_tree_recursively(tree.right, newick, named_parent_nodes, sequence_len)}){parent_node_name}:{branch_length}"
    else:
        # Tree has two children
        if tree.parent is not None and tree.parent.parent is not None:
            # Branch length for internal branches:
            # d(A,BC) = ( d(A,B) + d(A,D) + d(B,C) + d(B,D) ) / 4 - ( d(A,B) + d(C,D) ) / 2
            N = tree
            P = tree.parent
            R = tree.parent.parent
            A = N.left
            B = N.right
            C = P.right if P.right != N else P.left
            D = R.right if R.right != P else R.left
            branch_length = round(
                ((d(A, C) + d(A, D) + d(B, C) + d(B, D)) / 4 - (d(A, B) + d(C, D)) / 2) / sequence_len, 3)
        return f"({traverse_tree_recursively(tree.left, newick, named_parent_nodes, sequence_len)},{traverse_tree_recursively(tree.right, newick, named_parent_nodes, sequence_len)}){parent_node_name}:{branch_length}"


def to_newick(tree, named_parent_nodes=True, sequence_len=1):
    # # ----- DEBUG -----
    # test_tree = temp_generate_test_tree()
    # tree = test_tree
    # # Expected outcome: (B:0.2,(C:0.3,D:0.4)E:0.5)F;
    # # --- END DEBUG ---

    newick = ""
    newick = traverse_tree_recursively(tree, newick, named_parent_nodes, sequence_len=sequence_len)
    newick = f"{newick};"
    return newick
