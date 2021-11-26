import numpy as np

import datastructures

# u(i) = up-distance: "The average distance of the node from its children."
# Δ(i,j) = profile distance
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
    with open('data/test-small.aln') as f:
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
        tree = datastructures.Tree()
        tree.name = name
        tree.profile = sequence_to_profile(s)
        tree.calculate_out_distance(total_profile, active_nodes)
        nodes.append(tree)
    return nodes


def run(sequences, sequence_length):
    total_profile = compute_total_profile(sequences, sequence_length)
    print('Total profile:')

    for nt_profile in total_profile:
        print(nt_profile)

    nodes = sequences_to_trees(sequences, total_profile)

    profile0 = sequence_to_profile(sequences['>0'])
    profile1 = sequence_to_profile(sequences['>1'])
    # d = sequence_distance_uncorrected(profile0, profile1)
    # print(d)


if __name__ == '__main__':
    sequences, sequence_length = parse_input()
    run(sequences, sequence_length)

    # print(sequence_distance_uncorrected(sequences['>0'], sequences['>1']))
