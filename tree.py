import distance_measures


class Tree(object):
    def __init__(self):
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
        self.out_distance = (num_active_nodes * distance_measures.delta(self.profile, total_profile) \
                            - (self.up_distance * 2)) / (num_active_nodes - 2)

