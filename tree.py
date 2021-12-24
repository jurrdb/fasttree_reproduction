import numpy as np

import distance_measures


class Tree(object):
    def __init__(self):
        self.parent = None
        self.left = None
        self.right = None
        self.name = None
        self.profile = None
        self.top_hit_list = []
        self.best_known_join = None

        self.up_distance = 0
        self.out_distance = 0

    def is_leaf(self):
        # The tree has no children, it is a leaf node
        return not self.left and not self.right

    def calculate_up_distance(self):
        if self.is_leaf():
            return 0
        self.up_distance = distance_measures.delta(self.left.profile, self.right.profile) / 2

    def calculate_out_distance(self, total_profile, num_active_nodes):
        """
        :return: out_distance = r(i) = ( n * delta(i, T) - delta(i,i) ) / (n - 2)
        """
        # "delta(i,i) is the average distance between children of i including self-comparisons"
        if self.left and self.right:
            delta_i_i = distance_measures.delta(self.left.profile, self.right.profile)
        elif self.left or self.right:
            delta_i_i = 2*self.up_distance
        else:
            delta_i_i = 0
        self.out_distance = (num_active_nodes * distance_measures.delta(self.profile, total_profile) \
                             - delta_i_i) / np.max((num_active_nodes - 2, 1))
