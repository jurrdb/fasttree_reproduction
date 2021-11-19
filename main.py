import numpy as np

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





def sequence_distance_uncorrected(profileA, profileB):
    # Uncorrected distance is the fraction of places that differ

    return (profileA == profileB).sum() / profileA.shape[0]


def sequence_distance_corrected(profileA, profileB):
    # Calculates the corrected (Jukes-Cantor) distance
    return -3/4 * np.log(1 - 4 / 3 * sequence_distance_corrected(profileA, profileB))


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


def run(sequences, sequence_length):
    total_profile = compute_total_profile(sequences, sequence_length)
    print('Total profile:')
    for nt_profile in total_profile:
        print(nt_profile)

    profile0 = sequence_to_profile(sequences['>0'])
    profile1 = sequence_to_profile(sequences['>1'])
    d = sequence_distance_uncorrected(profile0, profile1)
    print(d)


if __name__ == '__main__':
    sequences, sequence_length = parse_input()
    run(sequences, sequence_length)
