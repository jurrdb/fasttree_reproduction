import numpy as np


def compute_total_profile(sequences, sequence_length):
    number_of_sequences = len(sequences)
    print('Number of sequences:', number_of_sequences)

    # Total profile is an array of dictionaries. Each index in the array represents a position in the sequence, and
    # each dictionary represents the frequency of nucleotides in that position.
    total_profile = []
    for i in range(0, sequence_length):
        total_profile.append({
            'A': 0.0,
            'C': 0.0,
            'T': 0.0,
            'G': 0.0,
        })

    for key, sequence in sequences.items():
        # Iterate over every sequence
        for i, nt in enumerate(sequence):
            total_profile[i][nt] += 1.0 / number_of_sequences

    return total_profile


def parse_input():
    with open('data/test-small.aln') as f:
        lines = f.readlines()
        sequences = {}
        for i in range(0, len(lines), 2):
            sequences[lines[i].strip()] =  lines[i + 1].strip()
        sequence_length = len(lines[1].strip())
        return sequences, sequence_length


def run(sequences, sequence_length):
    total_profile = compute_total_profile(sequences, sequence_length)
    print('Total profile:')
    for nt_profile in total_profile:
        print(nt_profile)


if __name__ == '__main__':
    sequences, sequence_length = parse_input()
    run(sequences, sequence_length)
