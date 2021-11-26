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


def join_nodes(left_child, right_child, total_profile, active_nodes):
    if not left_child.parent and not right_child.parent:
        # All nodes are active, start of building the tree
        parent = Tree()
        parent.name = left_child.name + right_child.name    # Concatenate names to create new parent name
        parent.left = left_child
        parent.right = right_child
        parent.profile = np.mean((left_child.profile, right_child.profile))
        parent.up_distance = parent.calculate_up_distance()
        parent.calculate_out_distance(total_profile, active_nodes)
        left_child.parent = parent
        right_child.parent = parent
        return parent
    elif left_child.parent or right_child.parent:
        # either node is part of the tree, the other node is active node that is joined into the tree
        # Jurrian: I wrote this part but we should still discuss if nodes can reference their parents
        node_with_parent = left_child if left_child.parent else right_child
        parent = Tree()
        if node_with_parent.parent.left == node_with_parent:
            node_with_parent.parent.left = parent
        else:
            node_with_parent.parent.right = parent
        parent.name = left_child.name + right_child.name  # Concatenate names to create new parent name
        parent.left = left_child
        parent.right = right_child
        parent.profile = np.mean(left_child.profile, right_child.profile)
        parent.up_distance = parent.calculate_up_distance()
        parent.calculate_out_distance(total_profile, active_nodes)
        # TODO is it now necessary to update the parent's parent up_distance / out_distance?
        left_child.parent = parent
        right_child.parent = parent
        return parent
    else:
        # both nodes have parents, this shouldn't happen
        print('ERROR, both joining nodes have parents, this shouldn\'t happen')
        return


class Tree(object):
    def __init__(self):
        self.parent = None
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
        if self.left:
            left_dist = self.left.up_distance
        if self.right:
            right_dist = self.right.up_distance
        return (left_dist + right_dist) / 2

    def calculate_out_distance(self, total_profile, active_nodes):
        self.out_distance = active_nodes * distance_measures.delta(self.profile, total_profile) \
                            - distance_measures.profile_distance_corrected(self, self)
