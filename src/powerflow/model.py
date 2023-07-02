import os
from typing import Any
import numpy as np
from enum import Enum

from powerflow.component import Component

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
        self.Pmin = 0
        self.Pmax = 100
        self.Qmin = 0
        self.Qmax = 100
        self.Vmin = 0
        self.Vmax = 100
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
        return f'Node {self.name}:\n\tType: {self.type}\n\tP: {self.P}\n\tQ: {self.Q}\n\tV: {self.V}\n\ttheta: {self.theta}'


class Branch:
    def __init__(self, name, type, node1: Node, node2: Node, Y=0+0j):
        self.name = name
        self.type = type
        self.node1 = node1
        self.node2 = node2
        self.node1.connect(self)
        self.node2.connect(self)
        self.Y = Y

    def __str__(self) -> str:
        return f'Branch {self.name}:\n\tNode1: {self.node1.name}\n\tNode2: {self.node2.name}\n\tR: {self.R}\n\tX: {self.X}\n\tB: {self.B}'


class Profile:
    def __init__(self, path):
        self.path = path
        self.data = open(path, 'r').readlines()
        self.data = [' '.join(line.strip().split()).split(' ')
                     for line in self.data]
        # drop lines start with '*'
        self.data = [line for line in self.data if line[0][0] != '*']

    def __str__(self) -> str:
        result = ''
        for line in self.data:
            result += ' '.join(line)+'\n'
        return result


class Model:
    def __init__(self):
        self.nodes = []
        self.branches = []

    def compose(self, profile: Profile):
        self.profile = profile
        self.componentManager = ComponentManager()
        self.componentManager.parseProfile(self.profile)

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



class ComponentManager:
    def __init__(self):
        self.components = []

    def parseProfile(self, profile):
        for i, v in enumerate(profile.data):
            self.parse(v)

    def parse(self, strList: list[str]):
        if (len(strList) == 0):
            return
        ctype = strList[0]
        match ctype:

            case 'SYSBASE':
                pass
            case 'SYSFREQ':
                pass
            case 'THLINE':
                self.addComponent(THLINE(strList))
                pass
            case 'LINE':
                pass
            case 'THTRFO':
                self.addComponent(THTRFO(strList))
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
                self.addComponent(GENER(strList))
                pass
            case 'GENERCV':
                self.addComponent(GENERCV(strList))
                pass
            case 'GENERDATA':
                # gener=self.findComponentByName(strList[1])
                # node=gener.node1
                # node.Pmax=float(strList[2])
                # node.Pmin=float(strList[3])
                # node.Qmax=float(strList[4])
                # node.Qmin=float(strList[5])
                # node.Vmax=float(strList[6])
                # node.Vmin=float(strList[7])
                pass
            case 'THLOAD':
                self.addComponent(THLOAD(strList))
                pass
            case 'LOAD2':
                # TODO: add load
                pass
            case 'THSLACK':
                gener: GENER = self.findComponentByName(strList[1])
                node: Node = gener.node1
                node.changeType(NodeType.SLACK)
                node.theta = 0
                pass
            case 'SLACKPH':
                # TODO: add slack
                pass
            case 'THSHUNT':
                self.addComponent(THSHUNT(strList))
                pass

    def addComponent(self, component):
        self.components.append(component)

    def findComponentByName(self, name):
        for component in self.components:
            if component.name == name:
                return component
        return None


class THLINE(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.node2 = strList[3]
        self.R = float(strList[4])
        self.X = float(strList[5])
        self.nBf2 = float(strList[6])
        print(f'Parsing {self.name}...')

    def apply(self, model: Model):
        node1: Node = model.findNodeByName(self.node1)
        node2: Node = model.findNodeByName(self.node2)
        branchLine = Branch(
            self.name, node1, node2, Y=1/(self.R + self.X*1j))
        node1.Ys += self.nBf2 * 1j
        node2.Ys += self.nBf2 * 1j
        model.addBranches(branchLine)


class THTRFO(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.node2 = strList[3]
        self.R = float(strList[4])
        self.X = float(strList[5])
        self.k = float(strList[6])
        print(f'Parsing {self.name}...')

    def apply(self, model: Model):
        node1: Node = model.findNodeByName(self.node1)
        node2: Node = model.findNodeByName(self.node2)
        branchT = Branch(
            self.name, node1, node2, Y=1 / ((self.R + self.X*1j) / self.k))
        node1.Ys += (self.k - 1) / (self.k * (self.R + self.X*1j))
        node2.Ys += (1 - self.k) / (self.k**2 * (self.R + self.X*1j))
        model.addBranches(branchT)


class THSHUNT(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.G = float(strList[3])
        self.B = float(strList[4])
        print(f'Parsing {self.name}...')

    def apply(self, model: Model):
        node1 = model.findNodeByName(self.node1)
        node1.Ys += self.G + self.B * 1j


class THLOAD(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.P = float(strList[3])
        self.Q = float(strList[4])
        print(f'Parsing {self.name}...')

    def apply(self, model: Model):
        node1: Node = model.findNodeByName(self.node1)
        node1.P += self.P
        node1.Q += self.Q


class GENERCV(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.P = float(strList[3])
        self.Q = float(strList[4])
        self.V = float(strList[5])
        self.state = int(strList[6])
        print(f'Parsing {self.name}...')

    def apply(self, model: Model):
        node1: Node = model.findNodeByName(self.node1)
        node1.P += self.P
        node1.V = self.V
        node1.changeType(NodeType.PV)


class GENER(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.P = float(strList[3])
        self.Q = float(strList[4])
        self.state = int(strList[5])
        print(f'Parsing {self.name}...')

    def apply(self, model: Model):
        node1: Node = model.findNodeByName(self.node1)
        node1.P += self.P
        node1.Q += self.Q
        node1.changeType(NodeType.PQ)

