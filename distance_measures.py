import numpy as np
import functools


def sequence_distance_uncorrected(seq_A: str, seq_B: str):
    # Uncorrected distance is the fraction of places that differ
    return sum([x == y for (x, y) in zip(seq_A, seq_B)]) / len(seq_A)


def sequence_distance_corrected(seq_A: str, seq_B: str):
    # Calculates the corrected (Jukes-Cantor) distance
    return -3 / 4 * np.log(1 - 4 / 3 * sequence_distance_uncorrected(seq_A, seq_B))


def out_distance():
    # Calculates r(i) = \sum_{k!=i} d_u(i,k) / (n-2)
    pass


def delta(Profile_A: np.ndarray, Profile_B: np.ndarray):
    # Calculate the distance between two profiles using Euclidean distance
    return np.linalg.norm(Profile_A - Profile_B, 1)


@functools.lru_cache(maxsize=None)
def profile_distance_uncorrected(A, B):
    return delta(A.profile, B.profile) - A.up_distance - B.up_distance


def profile_distance_corrected(A, B):
    # Currently we do not use the corrected transform
    return profile_distance_uncorrected(A, B)

