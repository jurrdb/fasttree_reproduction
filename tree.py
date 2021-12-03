import numpy as np

import distance_measures


# def insert(self, data):
#     if self.data:
#         if data < self.data:
#             if self.left is None:
#                 self.left = Node(data)
#             else:
#                 self.left.insert(data)
#         else data > self.data:
#             if self.right is None:
#                 self.right = Node(data)
#             else:
#                 self.right.insert(data)
#     else:
#         self.data = data


def join_nodes(left_child, right_child, total_profile, num_active_nodes):
    # Construct new tree node
    parent = Tree()
    parent.name = left_child.name + right_child.name  # Concatenate names to create new parent name
    parent.left = left_child
    parent.right = right_child
    parent.profile = np.mean((left_child.profile, right_child.profile), axis=0)
    parent.calculate_up_distance()
    parent.calculate_out_distance(total_profile, num_active_nodes)

    return parent  # should return the new parent as well as the new active_nodes


class Tree(object):
    def __init__(self):
        self.left = None
        self.right = None
        self.name = None
        self.profile = None
        # self.top_hit_list = []
        self.up_distance = 0
        self.out_distance = 0

    def is_leaf(self):
        # The tree has no children, it is a leaf node
        return not self.left and not self.right

    def calculate_up_distance(self):
        if self.is_leaf():
            return 0
        left_dist = right_dist = 0
        self.up_distance = distance_measures.delta(self.left.profile, self.right.profile) / 2

    def calculate_out_distance(self, total_profile, num_active_nodes):
        self.out_distance = num_active_nodes * distance_measures.delta(self.profile, total_profile) \
                            - (self.up_distance * 2)
