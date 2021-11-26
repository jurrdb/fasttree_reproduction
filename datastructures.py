import distance_measures


class Tree(object):
    def __init__(self):
        self.left = None
        self.right = None
        self.name = None
        self.profile = None
        self.up_distance = 0
        self.out_distance = 0

    def is_leaf(self):
        # The tree has no children, it is a leaf node
        return not self.left and not self.right

    def calculate_up_distance(self):
        if self.is_leaf():
            return 0
        left_dist = right_dist = 0
        if self.left:
            left_dist = self.left.calculate_up_distance()
        if self.right:
            right_dist = self.right.calculate_up_distance()
        return (left_dist + right_dist) / 2

    def calculate_out_distance(self, total_profile, active_nodes):
        self.out_distance = active_nodes * distance_measures.delta(self.profile, total_profile) \
               - distance_measures.profile_distance_corrected(self, self)
