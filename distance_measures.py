import numpy as np
import functools

def is_leaf():
    pass

def sequence_distance_uncorrected(seq_A, seq_B):
    # Uncorrected distance is the fraction of places that differ
    return sum([x == y for (x, y) in zip(seq_A, seq_B)]) / len(seq_A)


def sequence_distance_corrected(seq_A, seq_B):
    # Calculates the corrected (Jukes-Cantor) distance
    return -3 / 4 * np.log(1 - 4 / 3 * sequence_distance_uncorrected(seq_A, seq_B))



@functools.lru_cache(maxsize=None)
def profile_distance_uncorrected(A, B):
    if A.is_leaf() and B.is_leaf():
        return sequence_distance_uncorrected(A, B)
    elif B.is_leaf():
        # d_u(A, B) = d(A.left, B) + d(A.right, B)
        return 0
    elif A.is_leaf():
        # d_u(A, B) = d(B.left, A) + d(B.right, A)
        return 0

    return profile_distance_uncorrected(A.left, A.right) / 2 + profile_distance_uncorrected(B.left, B.right) / 2


def profile_distance_corrected(profile_A, profile_B):
    pass
