import numpy as np

import distance_measures as dm
import tree

# u(i) = up-distance: "The average distance of the node from its children."
# Δ(i,j) = profile distance
# P(AB) = average profile
# d_u(i,j) = Distance between internal nodes: d_u(i,j) = Δ(i,j) - u(i) - u(j)

# NOT USED in Fast Tree:
# r(i) = out-distance: "average out distance of i to other active nodes"


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


def parse_input():
    with open('data/fasttree-input.aln') as f:
        lines = f.readlines()
        sequences = {}
        for i in range(0, len(lines), 2):
            sequences[lines[i].strip()] = lines[i + 1].strip()
        sequence_length = len(lines[1].strip())
        return sequences, sequence_length


def sequences_to_trees(sequences, total_profile):
    nodes = list()
    active_nodes = len(sequences)
    for name, s in sequences.items():
        node = tree.Tree()
        node.name = name
        node.profile = sequence_to_profile(s)
        node.calculate_out_distance(total_profile, active_nodes)
        nodes.append(node)
    return nodes




def interchange_nodes(tree):
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
    # Algorithm Postorder(tree)
    #    1. Traverse the left subtree, i.e., call Postorder(left-subtree)
    interchange_nodes(tree.left)
    #    2. Traverse the right subtree, i.e., call Postorder(right-subtree)
    interchange_nodes(tree.right)
    #    3. Visit the root.
    # perform interchange (get relevant subtrees A B C D and calculate best topology)
    try:
        A = tree.left.left.left
        B = tree.left.left.right
        C = tree.left.right
        D = tree.right

    except AttributeError:
        print("Error")
    if A and B and C and D:

        N = tree.left.left
        # Calculate distances
        # d(A,B) + d(C,D) < d(A,C) + d(B,D) and d(A,B) + d(C,D) < d(A,D) + d(B,C)
        dist_ABCD = dm.profile_distance_corrected(A, B) + dm.profile_distance_corrected(C, D)
        dist_ACBD = dm.profile_distance_corrected(A, C) + dm.profile_distance_corrected(B, D)
        dist_ADBC = dm.profile_distance_corrected(A, D) + dm.profile_distance_corrected(B, C)

        if dist_ACBD < dist_ABCD and dist_ACBD < dist_ADBC:
            # dist_ACBD is smallest, B and C swapped
            N.right = C
            tree.left.right = B
        elif dist_ADBC < dist_ABCD and dist_ADBC < dist_ACBD:
            # dist_ADBC is smallest, D takes place of B, B takes place of C and C takes place of D
            N.right = D
            tree.right = C
            tree.left.right = B
        # Else dist_ABCD is smallest, do nothing
    return tree


def nearest_neighbour_interchange():
    # for i in range(4 * np.log(n, 2)):
    pass


def run(sequences, sequence_length):
    total_profile = compute_total_profile(sequences, sequence_length)
    print('Total profile:')

    for nt_profile in total_profile:
        print(nt_profile)

    active_nodes = sequences_to_trees(sequences, total_profile)

    while len(active_nodes) > 1:
        # Find best candidates for join using criterion
        min = np.inf
        pair = ()
        for i, node_a in enumerate(active_nodes):
            for node_b in active_nodes[i + 1:]:
                dis = dm.profile_distance_corrected(node_a, node_b) - node_a.out_distance - node_b.out_distance
                if dis < min:
                    min = dis
                    pair = (node_a, node_b)
        # Join nodes, update values
        active_nodes.remove(pair[0])
        active_nodes.remove(pair[1])
        parent = tree.join_nodes(pair[0], pair[1], total_profile, len(active_nodes))
        active_nodes.append(parent)

    # profile0 = sequence_to_profile(sequences['>0'])
    # profile1 = sequence_to_profile(sequences['>1'])
    # d = sequence_distance_uncorrected(profile0, profile1)
    return active_nodes[0]


if __name__ == '__main__':
    sequences, sequence_length = parse_input()
    run(sequences, sequence_length)

    # print(sequence_distance_uncorrected(sequences['>0'], sequences['>1']))
