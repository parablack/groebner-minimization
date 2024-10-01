import itertools as it
from .timing import *
from .setcover import *
import time

def simplify_not_not(node):
    if not isinstance(node, UnaryNode):
        return node
    if not isinstance(node.child, UnaryNode):
        return node
    if node.symbol == "~" and node.child.symbol == "~":
        return node.child.child
    else:
        return node

def simplify_binary_one_child(node):
    if not isinstance(node, BinaryNode):
        return node
    if len(node.children) == 1:
        return node.children[0]
    else:
        return node

def simplify_binary_same_child(node):
    if not isinstance(node, BinaryNode):
        return node
    children = []
    for child in node.children:
        if isinstance(child, BinaryNode) and child.symbol == node.symbol:
            children.extend(child.children)
        else:
            children.append(child)
    return BinaryNode(node.symbol, children)

def simplify_all(node):
    node = simplify_not_not(node)
    node = simplify_binary_one_child(node)
    node = simplify_binary_same_child(node)
    return node

def modify(node, func):
    if isinstance(node, BinaryNode):
        node.children = [modify(child, func) for child in node.children]
    if isinstance(node, UnaryNode):
        node.child = modify(node.child, func)
    return func(node)

class BinaryNode:

    def __init__(self, symbol, children):
        self.symbol = symbol
        self.children = children
    def __repr__(self):
        return f"BinaryNode(\"{self.symbol}\", {repr(self.children)})"

    def __str__(self):
        str_children = [str(child) for child in self.children]
        if len(str_children) == 0:
            if self.symbol == "|": return "0"
            if self.symbol == "&": return "1"
            else: raise ValueError("Unknown binary operator: {}".format(self.symbol))
        if len(str_children) == 1:
            return str_children[0]
        else:
            return "(" + (" " + self.symbol + " ").join(str_children) + ")"

    def copy(self):
        return BinaryNode(self.symbol, [child.copy() for child in self.children])

    def size(self):
        return 1 + sum(child.size() for child in self.children)

    def find(self, symbol):
        if self.symbol == symbol:
            yield self
        else:
            for child in self.children:
                yield from child.find(symbol)

    def leafs(self):
        for child in self.children:
            yield from child.leafs()

    def traverse(self, n, stats, depth=0):
        if len(self.children) == 0:
             return
        if stats.verbose:
            print("-" * depth, self.symbol)
        if self.symbol != "^":
            for child in self.children:
                child.traverse(n, stats, depth + 1)
        if self.symbol in ["|", "&"]:
            subsets = []
            weights = []
            for child in self.children:
                weights.append(child.size())
                formula = eval("lambda x: " + str(child) + " & 1")
                if self.symbol == "|":
                    subset = [value for value in it.product([0, 1], repeat=n) if formula(value)] # TODO n!
                else:
                    subset = [value for value in it.product([0, 1], repeat=n) if not formula(value)] # TODO n!
                subsets.append(subset)
            time_start = time.time()
            cover = set_cover(subsets, weights)
            stats.add_time_in_setcover(time.time() - time_start)
            self.children = [child for i, child in enumerate(self.children) if cover[i]]


class UnaryNode:

    def __init__(self, symbol, child):
        self.symbol = symbol
        self.child = child

    def __repr__(self):
        return f"UnaryNode(\"{self.symbol}\", {repr(self.child)})"

    def __str__(self):
        return self.symbol + str(self.child)

    def copy(self):
        return UnaryNode(self.symbol, self.child.copy())

    def size(self):
        return 1 + self.child.size()

    def find(self, symbol):
        if self.symbol == symbol:
            yield self
        else:
            yield from self.child.find(symbol)

    def leafs(self):
        yield from self.child.leafs()

    def traverse(self, n, stats, depth=0):
        global verbose
        if stats.verbose:
            print("-" * depth, self.symbol)
        self.child.traverse(n, stats, depth + 1)

class LeafNode:

    def __init__(self, value):
        self.value = value

    def __repr__(self):
        return f"LeafNode({repr(self.value)})"

    def __str__(self):
        return str(self.value)

    def copy(self):
        return LeafNode(self.value)

    def size(self):
        return 1

    def leafs(self):
        yield self

    def traverse(self, n, stats,  depth=0):
        global verbose
        if stats.verbose:
            print("-" * depth, self.value)

    def find(self, symbol):
        yield from []


def And(children):
    return BinaryNode("&", children)

def Or(children):
    return BinaryNode("|", children)

def Xor(children):
    return BinaryNode("^", children)

def Not(child):
    return UnaryNode("~", child)

def Value(value):
    return LeafNode(value)