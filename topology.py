import copy
import math

import numpy as np

import distance_measures as dm
from tree import Tree

nt_map = {
    'A': 0,
    'C': 1,
    'T': 2,
    'G': 3,
}


def sequence_to_profile(sequence):
    sequence_length = len(sequence)
    profile = np.full((sequence_length, 4), 0.0)

    for i, nt in enumerate(sequence):
        index = nt_map[nt]
        profile[i][index] = 1.0

    return profile


def compute_total_profile(sequences, sequence_length):
    number_of_sequences = len(sequences)
    print('Number of sequences:', number_of_sequences)

    # Total profile is an array of dictionaries. Each index in the array represents a position in the sequence, and
    # each dictionary represents the frequency of nucleotides in that position.
    # pseudo_total = number_of_sequences + 4
    # pseudo_freq = 1.0 / pseudo_total
    total_profile = np.full((sequence_length, 4), 0.0)

    for key, sequence in sequences.items():
        # Iterate over every sequence
        for i, nt in enumerate(sequence):
            index = nt_map[nt]
            total_profile[i][index] += 1.0 / number_of_sequences

    return total_profile


def compute_total_profile_active_nodes(active_nodes):
    """
    Computes the total_profile: the average profile of all active nodes.
    :return: total_profile
    """
    sequence_length = len(active_nodes[0].profile)
    total_profile = np.full((sequence_length, 4), 0.0)
    for node in active_nodes:
        total_profile += node.profile
    total_profile /= len(active_nodes)


def sequences_to_trees(sequences, total_profile):
    """
    Converts each nucleotide sequence to a list of Tree node objects. Each Tree will initially contain only a single
    node
    :param sequences: dictionary containing nucleotide sequences
    :param total_profile: average profile of all nucleotides
    :return: list of Tree nodes
    """
    nodes = list()
    num_active_nodes = len(sequences)
    for name, s in sequences.items():
        node = Tree()
        node.name = name
        node.profile = sequence_to_profile(s)
        node.calculate_out_distance(total_profile, num_active_nodes)
        nodes.append(node)
    return nodes


def calculate_top_hits(node, node_list, target_length):
    nodes_copy = node_list.copy()
    if node in nodes_copy:
        nodes_copy.remove(node)
    sorted_nodes = sorted(nodes_copy, key=lambda node_b:
        dm.profile_distance_corrected(node, node_b) - node.out_distance - node_b.out_distance)[0:int(target_length)+1]
    node.top_hit_list = sorted_nodes
    return


def calculate_top_hits_active_nodes(active_nodes, m):
    """
    Find the 2m closest node for each node (the top hit list) according to the neighbor joining criterion
    :param active_nodes: nodes that are being considered for joining
    :param m: number of nodes in each hit list (Note: should not include safety factor of 2). By default this is sqrt(N)
    """
    target_length = min(2 * int(m), len(active_nodes))
    # initial_node = active_nodes[0]
    node_copy = copy.copy(active_nodes)
    # node_copy.remove(initial_node)
    while len(node_copy) > 0:
        initial_node = node_copy[0]
        calculate_top_hits(initial_node, active_nodes, target_length)
        node_copy.remove(initial_node)
        for b in initial_node.top_hit_list:
            if dm.profile_distance_uncorrected(initial_node, b) < 0.75*dm.profile_distance_uncorrected(initial_node, initial_node.top_hit_list[target_length]):
                calculate_top_hits(b, initial_node.top_hit_list, target_length)
                if b in node_copy:
                    node_copy.remove(b)

    return


def join_nodes(left_child, right_child, total_profile, num_active_nodes):
    # Construct new tree node
    parent = Tree()
    parent.name = left_child.name + right_child.name  # Concatenate names to create new parent name
    parent.left = left_child
    parent.right = right_child
    left_child.parent = parent
    right_child.parent = parent

    parent.profile = np.mean((left_child.profile, right_child.profile), axis=0)
    parent.calculate_up_distance()
    parent.calculate_out_distance(total_profile, num_active_nodes)

    joined_hit_list = list(set(left_child.top_hit_list + right_child.top_hit_list))

    if left_child in joined_hit_list:
        joined_hit_list.remove(left_child)
    if right_child in joined_hit_list:
        joined_hit_list.remove(right_child)

    calculate_top_hits(parent, joined_hit_list, math.sqrt(num_active_nodes))

    return parent  # should return the new parent as well as the new active_nodes


def try_interchange(tree, counter, side_a, side_b):
    try:
        # To deal with tree symmetry we use getters that can get both the left and right sided subtree nodes.
        A = tree.__getattribute__(side_a).__getattribute__(side_a).__getattribute__(side_a)
        B = tree.__getattribute__(side_a).__getattribute__(side_a).__getattribute__(side_b)
        C = tree.__getattribute__(side_a).__getattribute__(side_b)
        D = tree

        if A and B and C and D:
            # reached count?
            P = tree.__getattribute__(side_a)
            N = tree.__getattribute__(side_a).__getattribute__(side_a)
            # Calculate distances
            # d(A,B) + d(C,D) < d(A,C) + d(B,D) and d(A,B) + d(C,D) < d(A,D) + d(B,C)
            dist_ABCD = dm.profile_distance_corrected(A, B) + dm.profile_distance_corrected(C, D)
            dist_ACBD = dm.profile_distance_corrected(A, C) + dm.profile_distance_corrected(B, D)
            dist_ADBC = dm.profile_distance_corrected(A, D) + dm.profile_distance_corrected(B, C)

            if dist_ACBD < dist_ABCD and dist_ACBD < dist_ADBC:
                print("performed interchange (" + side_a + ")")
                # dist_ACBD is smallest, B and C swapped
                N.__setattr__(side_b, C)
                P.__setattr__(side_b, B)
                C.parent = N
                B.parent = P
                counter.count += 1
            elif dist_ADBC < dist_ABCD and dist_ADBC < dist_ACBD:
                print("performed interchange (" + side_a + ")")
                # dist_ADBC is smallest, D takes place of B, B takes place of C and C takes place of D
                N.__setattr__(side_b, D)
                D.parent = N
                tree.__setattr__(side_b, C)
                C.parent = tree
                tree.__getattribute__(side_a).__setattr__(side_b, B)
                B.parent = P
                counter.count += 1
            # Else dist_ABCD is smallest, do nothing
            else:
                print("performed no interchange (" + side_a + ")")
    except AttributeError:
        print("interchange_nodes: skipping leaf node (" + side_a + ")")


def interchange_nodes(tree, counter):
    """
             R
            / \
           P   D
          / \
          N  C
         /\
        A B
-----------------------
             R
            / \
           P   D
          / \
          N  B
         /\
        A C
-------------------------
             R
            / \
           P   C
          / \
          N  B
         /\
        A D
    For each non-root node N, with children A,B, sibling C, and uncle D,
     we compare the current topology AB|CD to the alternate topologies
     AC|BD and AD|BC, by using the 4 relevant profiles.
    """
    if not tree or tree.is_leaf() or counter.count >= counter.max_rounds:
        return tree

    # Algorithm Postorder(tree)
    #    1. Traverse the left subtree, i.e., call Postorder(left-subtree)
    tree.left = interchange_nodes(tree.left, counter)
    #    2. Traverse the right subtree, i.e., call Postorder(right-subtree)
    tree.right = interchange_nodes(tree.right, counter)
    #    3. Visit the root.
    # perform interchange (get relevant subtrees A B C D and calculate best topology)

    try_interchange(tree, counter, 'right', 'left')
    try_interchange(tree, counter, 'left', 'right')
    return tree
