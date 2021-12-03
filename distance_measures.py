import numpy as np


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


# @functools.lru_cache(maxsize=None)
def profile_distance_uncorrected(A, B):
    return delta(A.profile, B.profile) - A.up_distance - B.up_distance

    # if A.is_leaf() and B.is_leaf():
    #     # d_u(i,j) = delta (i,j) - u(i) - u(j)
    #     return delta(A.profile, B.profile) - A.calculate_up_distance() - B.calculate_up_distance()
    #     # return sequence_distance_uncorrected(A, B)
    # elif B.is_leaf():
    #     return profile_distance_uncorrected(A.left, B) + profile_distance_uncorrected(A.right, B)
    #     # d_u(A, B) = d(A.left, B) + d(A.right, B)
    # elif A.is_leaf():
    #     return profile_distance_uncorrected(B.left, A) + profile_distance_uncorrected(B.right, A)
    #     # d_u(A, B) = d(B.left, A) + d(B.right, A)
    #
    # return profile_distance_uncorrected(A.left, A.right) / 2 + profile_distance_uncorrected(B.left, B.right) / 2


def profile_distance_corrected(A, B):
    # TODO correct for the proportion of nongaps (See section 'Distance Between Profiles')
    return profile_distance_uncorrected(A, B)
