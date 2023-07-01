import numpy as np
from model import Model, Branch, Node, NodeType

class Component(object):
    name = None

    def apply(model):
        pass


class THLINE(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.node2 = strList[3]
        self.R = float(strList[4])
        self.X = float(strList[5])
        self.nBf2 = float(strList[6])

    def apply(self, model: Model):
        node1 = model.findNodebyName(self.node1)
        node2 = model.findNodeByName(self.node2)
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

    def apply(self, model):
        node1 = model.getNode(self.node1)
        node2 = model.getNode(self.node2)
        branchLine = Branch(
            self.name, node1, node2, R=self.R, X=self.X)
        node1.Ys += self.nBf2
        node2.Ys += self.nBf2
        model.addBranches(branchLine)

class THSHUNT(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.R = float(strList[4])
        self.X = float(strList[5])
        self.k = float(strList[6])

    def apply(self, model):
        node1 = model.findNodebyName(self.node1)
        node2 = model.findNodeByName(self.node2)
        branchLine = Branch(
            self.name, node1, node2, R=self.R, X=self.X)
        node1.Ys += self.nBf2
        node2.Ys += self.nBf2
        model.addBranches(branchLine)

class THLOAD(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.P = float(strList[3])
        self.Q = float(strList[4])

    def apply(self, model):
        node1 = model.getNode(self.node1)
        node2 = model.getNode(self.node2)
        branchLine = Branch(
            self.name, node1, node2, R=self.R, X=self.X)
        node1.Ys += self.nBf2
        node2.Ys += self.nBf2
        model.addBranches(branchLine)

class ComponentManager:
    def __init__(self):
        self.components = []

    def parseProfile(self, profile):
        for i, v in enumerate(profile.data):
            self.parse(v)

    def parse(self, strList):
        if(len(strList) == 0): return
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
                node1 = self.getNode(v[2])

                pass
            case 'GENERCV':
                
                pass
            case 'GENERDATA':
                # TODO: add generator
                pass
            case 'THLOAD':
                
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
                pass

    def addComponent(self, component):
        self.components.append(component)

    def findComponentByName(self, name):
        for component in self.components:
            if component.name == name:
                return component
        return None

