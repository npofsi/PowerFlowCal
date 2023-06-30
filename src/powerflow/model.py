import os
from typing import Any
import numpy as np
from enum import Enum


class NodeType(Enum):
    PQ = 1
    PV = 2
    Slack = 3

class BranchType(Enum):
    Impedance = 1
    Addmitance = 2
    Transformer = 3
    Power = 4


class Node:
    def __init__(self, name, type, P=0, Q=0, V=0, theta=0, canChangeType=True):
        self.name = name
        self.type = type
        self.canChangeType = canChangeType
        self.P = P
        self.Q = Q
        self.V = V
        self.theta = theta

    def changeType(self, type):
        if not self.canChangeType:
            self.type = type

    def __str__(self) -> str:
        print(
            f'Node {self.name}:\n\tType: {self.type}\n\tP: {self.P}\n\tQ: {self.Q}\n\tV: {self.V}\n\ttheta: {self.theta}')


class Branch:
    def __init__(self, name, type, node1name, node2name, R=0, X=0, G=0, B=0, k=1, P=0, Q=0):
        self.name = name
        self.type = type
        self.node1name = node1name
        self.node2name = node2name
        self.R = R
        self.X = X
        self.G = G
        self.B = B
        self.k = k
        self.P = P
        self.Q = Q

    def __str__(self) -> str:
        print(f'Branch {self.name}:\n\tNode1: {self.node1name}\n\tNode2: {self.node2name}\n\tR: {self.R}\n\tX: {self.X}\n\tB: {self.B}')


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

    def getNode(self, name):
        for node in self.nodes:
            if node.name == name:
                return node
        node = Node(name, NodeType.PQ)
        self.nodes.append(node)
        return node

    def findNodeByName(self, name):
        for node in self.nodes:
            if node.name == name:
                return node
        return None

    def findBranchByName(self, name):
        for branch in self.branches:
            if branch.name == name:
                return branch
        return None

    def deriveYMatrix(self):
        pass

    def __str__(self) -> str:
        for node in self.nodes:
            print(node)
        for branch in self.branches:
            print(branch)

    def compose(self,profile):
      for i, v in enumerate(profile.data):
          ctype = v[0]
          match ctype:
              case 'SYSBASE':
                  pass
              case 'SYSFREQ':
                  pass
              case 'THLINE':
                  node1 = self.getNode(v[2])
                  node2 = self.getNode(v[3])
                  branchLine = Branch(
                      v[1], node1.name, node2.name, R=float(v[4]), X=float(v[5]))
                  branchB1 = Branch(v[1]+'_B1', node1.name,
                                    self.groundNode.name, B=float(v[6]))
                  branchB2 = Branch(v[1]+'_B2', node2.name,
                                    self.groundNode.name, B=float(v[6]))
                  self.addBranches(branchLine, branchB1, branchB2)
                  pass
              case 'LINE':

                  # TODO: add line
                  pass
              case 'THTRFO':
                  node1 = self.getNode(v[2])
                  node2 = self.getNode(v[3])
                  branchTransformer = Branch(v[1], node1.name, node2.name, R=float(
                      v[4]), X=float(v[5]), k=float(v[6]))
                  self.addBranches(branchTransformer)
                  pass
              case 'TRFO':
                  # TODO: add transformer
                  pass
              case 'THFORB2':
                  pass
              case 'THTRPH':
                  pass
              case 'TAP':
                  pass
              case 'TAPCV':
                  pass
              case 'GENER':
                  node1 = self.getNode(v[2])
                  
                  pass
              case 'GENERCV':
                  node1 = self.getNode(v[2])
                  node1.changeType(NodeType.PV)
                  node1.P -= float(v[3])
                  node1.V = float(v[4])
                  pass
              case 'GENERDATA':
                  # TODO: add generator
                  pass
              case 'THLOAD':
                  node1 = self.getNode(v[2])
                  node1.changeType(NodeType.PQ)
                  node1.P -= float(v[3])
                  node1.Q -= float(v[4])
                  pass
              case 'LOAD2':
                  # TODO: add load
                  pass
              case 'THSLACK':
                  # TODO: add slack
                  pass
              case 'SLACKPH':
                  # TODO: add slack
                  pass
              case 'THSHUNT':
                  node1 = self.getNode(v[2])
                  branchShuntP = Branch(v[1], node1.name, self.groundNode.name, R=float(v[3]))
                  if float(v[4]) > 0:
                    branchShuntQ = Branch(v[1], node1.name, self.groundNode.name, X=float(v[4]))
                  else:
                    branchShuntQ = Branch(v[1], node1.name, self.groundNode.name, B=float(v[4]))
                  self.addBranches(branchShuntP, branchShuntQ)
                  pass
