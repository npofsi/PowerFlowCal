import os
from typing import Any
import numpy as np
from enum import Enum


class NodeType(Enum):
    PQ = 1
    PV = 2
    Slack = 3


class Node:
    def __init__(self, name, type, P=0, Q=0, V=0, Ys=0+0j, theta=0, canChangeType=True):
        self.name = name
        self.type = type
        self.canChangeType = canChangeType
        self.P = P
        self.Q = Q
        self.V = V
        self.Smin = 0
        self.Smax = 100
        self.Ymin = 0+0j
        self.Ymax = 100+100j
        self.Ys = Ys
        self.theta = theta
        self.connectedBranches = []

    def S(self):
        return self.P + self.Q * 1j

    def changeType(self, type):
        if not self.canChangeType:
            self.type = type

    def connect(self, branch):
        self.connectedBranches.append(branch)

    def __str__(self) -> str:
        print(
            f'Node {self.name}:\n\tType: {self.type}\n\tP: {self.P}\n\tQ: {self.Q}\n\tV: {self.V}\n\ttheta: {self.theta}')


class Branch:
    def __init__(self, name, type, node1, node2, Y=0+0j):
        self.name = name
        self.type = type
        self.node1 = node1
        self.node2 = node2
        self.node1.connect(self)
        self.node2.connect(self)
        self.Y = Y

    def __str__(self) -> str:
        print(f'Branch {self.name}:\n\tNode1: {self.node1.name}\n\tNode2: {self.node2.name}\n\tR: {self.R}\n\tX: {self.X}\n\tB: {self.B}')


class Profile:
    def __init__(self, path):
        self.path = path
        self.data = open(path, 'r').readlines()
        self.data = [' '.join(line.strip().split()).split(' ')
                     for line in self.data]
        # drop lines start with '*'
        self.data = [line for line in self.data if line[0][0] != '*']

    def __str__(self) -> str:
        for line in self.data:
            print(' '.join(line))


class Model:
    def __init__(self):
        self.groundNode = Node('GND', NodeType.PV, P=0, V=0)
        self.nodes = [self.groundNode]
        self.branches = []

    def addNodes(self, *anodes):
        for node in anodes:
            self.nodes.append(node)

    def addBranches(self, *abranches):
        for branch in abranches:
            for mbranch in self.branches:
                if branch.name == mbranch.name:
                    return None
            self.branches.append(branch)

    def findNodeByName(self, name):
        for node in self.nodes:
            if node.name == name:
                return node
        return None

    def findBranchByName(self, name):
        for branch in self.branches:
            if branch.name == name:
                return branch
        node = Node(name, NodeType.PQ)
        self.nodes.append(node)
        return node

    def deriveYMatrix(self):
        pass

    def compose(self, profile):
        pass