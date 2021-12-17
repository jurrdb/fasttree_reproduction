import copy
import math

import numpy as np

import distance_measures as dm
import helpers
from topology import sequences_to_trees, interchange_nodes, calculate_top_hits_active_nodes, compute_total_profile, \
                    compute_total_profile_active_nodes, join_nodes

# u(i) = up-distance: "The average distance of the node from its children."
# Δ(i,j) = profile distance
# P(AB) = average profile
# d_u(i,j) = Distance between internal nodes: d_u(i,j) = Δ(i,j) - u(i) - u(j)
# r(i) = out-distance: "average out distance of i to other active nodes"
#
# L = number of positions / sequence length


def parse_input():
    with open('data/fasttree-input.aln') as f:
        lines = f.readlines()
        sequences = {}
        for i in range(0, len(lines), 2):
            sequences[lines[i].strip()] = lines[i + 1].strip()
        sequence_length = len(lines[1].strip())
        return sequences, sequence_length


def run(sequences, sequence_length):
    total_profile = compute_total_profile(sequences, sequence_length)
    print('Total profile:')

    for nt_profile in total_profile:
        print(nt_profile)

    active_nodes = sequences_to_trees(sequences, total_profile)
    num_nodes = len(active_nodes)
    calculate_top_hits_active_nodes(active_nodes, math.sqrt(num_nodes))

    num_joined_nodes = 0
    active_nodes = sorted(active_nodes, key=lambda x: x.out_distance,reverse=True)
    while len(active_nodes) > 1:
        # Find best candidates for join using criterion
        min = np.inf
        pair = ()
        for i, node_a in enumerate(active_nodes):
            for node_b in node_a.top_hit_list:
                dis = dm.profile_distance_corrected(node_a, node_b) - node_a.out_distance - node_b.out_distance
                # TODO from paper: Traditional Neighbor-Joining computes all N time, and updates each out-distance
                #  after each join, which also takes O(N2) time overall
                #  To avoid this work, FastTree computes each out-distance as needed in O(La) time by
                #  using a “total profile” T which is the average of all active nodes’ profiles, as
                #  implied by: (See formula on p. 1644)
                if dis < min:
                    min = dis
                    pair = (node_a, node_b)

        # Join nodes, update values
        active_nodes.remove(pair[0])
        active_nodes.remove(pair[1])
        parent = join_nodes(pair[0], pair[1], total_profile, len(active_nodes))
        num_joined_nodes += 1
        active_nodes.append(parent)

        for node in active_nodes:
            for i, hit in enumerate(node.top_hit_list):
                if hit == pair[0] or hit == pair[1]:
                    node.top_hit_list[i] = parent

        if num_joined_nodes % 200 == 199:
            print('Recalculating total profile')
            # TODO test this
            total_profile = compute_total_profile_active_nodes(active_nodes)

    # interchange nodes postorder until log(N) + 1 rounds of interchanges
    initial_topology = copy.deepcopy(active_nodes[0])
    max_rounds = math.inf #math.log2(num_nodes) + 1
    counter = helpers.Counter(max_rounds=max_rounds)

    final_topology = interchange_nodes(active_nodes[0], counter)
    for i in range(int(math.log2(num_nodes))*4):
        final_topology = interchange_nodes(final_topology, counter)

    return initial_topology, final_topology


if __name__ == '__main__':
    sequences, sequence_length = parse_input()
    run(sequences, sequence_length)

    # print(sequence_distance_uncorrected(sequences['>0'], sequences['>1']))
