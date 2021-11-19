
class Tree(object):
    def __init__(self):
        self.left = None
        self.right = None
        self.data = None

    def is_leaf(self):
        # The tree has no children, it is a leaf node
        return not self.left and not self.right
