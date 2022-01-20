import functools
import numpy as np


def delta(Profile_A: np.ndarray, Profile_B: np.ndarray):
    # Calculate the distance between two profiles using Frobenius norm
    return np.linalg.norm(Profile_A - Profile_B, 2)


@functools.lru_cache(maxsize=None)
def profile_distance_uncorrected(A, B):
    return delta(A.profile, B.profile) - A.up_distance - B.up_distance


def profile_distance_corrected(A, B):
    # Currently we do not use the corrected transform
    return profile_distance_uncorrected(A, B)
