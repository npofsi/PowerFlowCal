import os
import numpy as np
from enum import Enum

from powerflow.component import Component
from powerflow.utils import P2C, P2Complex

global Sb
Sb = 100

#节点类型，使用枚举类型
class NodeType(Enum):
    PQ = 1
    PV = 2
    Slack = 3 #平衡节点

#节点模型
class Node:
    def __init__(self, name, type=NodeType.PQ, P=.0, Q=.0, V=1.+0j, Ys=.0+.0j, theta=0., canChangeType=True):
        self.name = name
        self.type = type #节点类型，默认是PQ节点
        self.canChangeType = canChangeType #是否可以改变节点类型

        self.P = P #节点总有功功率
        self.Pg = P #节点发电机有功功率
        self.Pd = 0. #节点负荷有功功率
        self.Q = Q #节点总无功功率
        self.Qg = Q #节点发电机无功功率
        self.Qd = 0. #节点负荷无功功率
        self.V = P2Complex(V, theta) #节点电压
        self.oV = V #原始节点电压幅值

        #节点功率和节点电压的最小值和最大值
        self.Pmin = -99999./Sb 
        self.Pmax = 99999./Sb
        self.Qmin = -99999./Sb
        self.Qmax = 99999./Sb
        self.Vmin = 0.
        self.Vmax = 100.

        #节点自导纳
        self.Ys = Ys

        #节点所连接的支路
        self.connectedBranches = []

    #V 是复数，计算相角需要一步转换
    def getTheta(self):
        return np.angle(self.V)
    #设置相角同上
    def setTheta(self, theta):
        self.V = P2Complex(self.V, theta)

    #输出节点总视在功率
    def S(self):
        return self.P + self.Q * 1j

    #更改节点类型
    def changeType(self, type):
        if self.canChangeType:
            self.type = type

    #添加节点所连接的支路
    def connect(self, branch):
        self.connectedBranches.append(branch)

    #计算节点未确定的功率，如PV和Slack节点
    def calSg(self):
        if self.type == NodeType.PQ:
            return 0 + 0 * 1j

        Flow = 0 + 0 * 1j
        for branch in self.connectedBranches:
            if branch.node1 == self:
                Flow += branch.Flow
            else:
                Flow -= branch.Flow

        self.Pg = self.Pd + self.P + Flow.real
        # if self.type == NodeType.Slack:
        #     return self.P + self.Q * 1j
        self.Qg = self.Qd + self.Q + Flow.imag
        return self.Pg + self.Qg * 1j

    #重写字符串输出
    def __str__(self) -> str:
        return f'\tNode {self.name}:\n\t\tType: {self.type}\n\t\tP: {self.P} r[{self.Pmin}~{self.Pmax}]\n\t\tQ: {self.Q} r[{self.Qmin}~{self.Qmax}]\n\t\tV: {np.abs(self.V)} ({self.V}) r[{self.Vmin}~{self.Vmax}]\n\t\ttheta: {self.getTheta()}({self.getTheta()/np.pi*180}d)'

#导纳支路模型
class Branch:
    def __init__(self, name, node1: Node, node2: Node, Y=0+0j):
        self.name = name

        #支路两端的节点
        self.node1 = node1
        self.node2 = node2
        self.node1.connect(self)
        self.node2.connect(self)

        #支路导纳
        self.Y = Y

        #支路电流，支路功率流，支路功率损耗
        self.I = 0+0j
        self.Flow = 0+0j
        self.Loss = 0+0j

        #支路额定电流
        self.Irated = 99999./Sb

    def __str__(self) -> str:
        return f'\tBranch {self.name}:\n\t\tNode1: {self.node1.name}\n\t\tNode2: {self.node2.name}\n\t\tY: {self.Y}'

#输入解析类，用于解析输入文件，解析每一行数据存入一个数组
class Profile:
    def __init__(self, path):
        self.path = path

        #数据存在data数组中，一行一个数组元素，每一行也是数组，数组里面的元素是字符串
        self.data = open(path, 'r').readlines()
        self.data = [' '.join(line.strip().split()).split(' ')
                     for line in self.data]
        # drop lines start with '*'
        self.data = [line for line in self.data if not (
            len(line[0]) == 0 or line[0][0] == '*')]

    def __str__(self) -> str:
        result = ''
        for line in self.data:
            result += ' '.join(line)+'\n'
        return result

#模型类，用于存储整个节点导纳模型
class Model:
    def __init__(self):
        #电网总损耗，
        #系统节点和支路
        self.loss = 0.+0.j
        self.nodes = []
        self.branches = []

    #生成模型
    def compose(self, profile: Profile):
        self.profile = profile
        self.componentManager = ComponentManager(self)
        #解析输入文件，分析元件，生成节点和支路
        self.componentManager.parseProfile(self.profile)

    #添加节点
    def addNodes(self, *anodes):
        for node in anodes:
            self.nodes.append(node)

    #添加支路
    def addBranches(self, *abranches):
        for branch in abranches:
            for mbranch in self.branches:
                if branch.name == mbranch.name:
                    return None
            self.branches.append(branch)

    #寻找节点，如果没有则创建
    def findNodeByName(self, name) -> Node:
        for node in self.nodes:
            if node.name == name:
                return node
        node = Node(name, NodeType.PQ)
        self.nodes.append(node)
        return node

    #寻找支路
    def findBranchByName(self, name) -> Branch:
        for branch in self.branches:
            if branch.name == name:
                return branch
        return None

    #计算节点导纳矩阵
    def deriveYMatrix(self):
        self.nodes.sort(key=lambda node: node.type.value)
        n = len(self.nodes)
        Y = np.zeros((n, n), dtype=complex)

        #遍历模型中的支路，逐渐添加支路到空的导纳矩阵中
        for branch in self.branches:
            i = self.nodes.index(branch.node1)
            j = self.nodes.index(branch.node2)
            Y[i, j] = -branch.Y
            Y[j, i] = -branch.Y
            Y[i, i] += branch.Y
            Y[j, j] += branch.Y

        #考虑节点自导纳对节点导纳矩阵的影响
        for node in self.nodes:
            i = self.nodes.index(node)
            Y[i, i] += node.Ys
        
        #打印节点导纳矩阵
        print('Y Matrix(.0f):')
        for i in range(n):
            print(f'{self.nodes[i].name}', end='\t')
            for j in range(n):
                print(f'{Y[i,j]:f}', end='\t')
            print()

        return Y

    def printTopology(self):
        print('Nodes:')
        for node in self.nodes:
            print(node)
        print('Branches:')
        for branch in self.branches:
            print(branch)

#元件管理类，用于存储元件的基本信息
class ComponentManager:
    def __init__(self, model: Model = None):
        self.components = []
        self.model = model

    #解析输入的数据，对行进行遍历
    def parseProfile(self, profile):
        for i, v in enumerate(profile.data):
            self.parse(v)

    #解析一行数据
    def parse(self, strList: list[str]):
        #跳过空行
        if (len(strList) == 0):
            return
        
        #根据第一个元素判断元件类型
        ctype = strList[0]
        match ctype:
            case 'SYSBASE':
                pass
            case 'SYSFREQ':
                pass
            #添加线缆元件
            case 'THLINE':
                self.addComponent(THLINE(strList)).apply(self.model)
                pass
            case 'LINE':
                self.addComponent(LINE(strList)).apply(self.model)
                pass
            #添加变压器元件
            case 'THTRFO':
                self.addComponent(THTRFO(strList)).apply(self.model)
                pass
            case 'TRFO':
                self.addComponent(TRFO(strList)).apply(self.model)
                pass
            case 'THFORB2':
                pass
            case 'THTRPH':
                pass
            case 'TAP':
                pass
            case 'TAPCV':
                pass
            #添加发电机元件
            case 'GENER':
                self.addComponent(GENER(strList)).apply(self.model)
                pass
            case 'GENERCV':
                self.addComponent(GENERCV(strList)).apply(self.model)
                pass
            #解析发电机额外数据
            case 'GENERDATA':
                gener = self.findComponentByName(strList[1])
                node = self.model.findNodeByName(gener.node1)
                node.Pmax = float(strList[7])/Sb
                node.Pmin = float(strList[8])/Sb
                node.Qmax = float(strList[9])/Sb
                node.Qmin = float(strList[10])/Sb
                node.Vmax = float(strList[11])
                node.Vmin = float(strList[12])
                pass
            #添加负荷元件
            case 'THLOAD':
                self.addComponent(THLOAD(strList)).apply(self.model)
                pass
            case 'LOAD2':
                self.addComponent(LOAD2(strList)).apply(self.model)
                pass
            #解析平衡节点数据，并应用于模型
            case 'THSLACK':
                print("Dealing with THSlack")
                gener: GENER = self.findComponentByName(strList[1])
                node: Node = self.model.findNodeByName(gener.node1)
                node.V = float(strList[2])
                node.setTheta(float(strList[3]))
                node.changeType(NodeType.Slack)
                node.canChangeType = False
                pass
            case 'SLACKPH':
                print("Dealing with SlackPH")
                gener: GENER = self.findComponentByName(strList[1])
                node: Node = self.model.findNodeByName(gener.node1)
                node.V = float(strList[4])
                node.setTheta(float(strList[3]))
                node.changeType(NodeType.Slack)
                node.canChangeType = False
                pass
            #添加并联功率元件
            case 'THSHUNT':
                self.addComponent(THSHUNT(strList)).apply(self.model)
                pass

    #添加元件到管理类
    def addComponent(self, component) -> Component:
        self.components.append(component)
        return component

    #根据元件名查找元件
    def findComponentByName(self, name):
        for component in self.components:
            if component.name == name:
                return component
        return None

#线路模型
class THLINE(Component):
    #解析线路数据
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.node2 = strList[3]
        self.R = float(strList[4])
        self.X = float(strList[5])
        self.nBf2 = float(strList[6])

        print(f'Parsing {self.name}...')

    #在所给的模型中创建线路对应的导纳模型
    def apply(self, model: Model):
        node1: Node = model.findNodeByName(self.node1)
        node2: Node = model.findNodeByName(self.node2)
        branchLine = Branch(
            self.name, node1, node2, Y=1/(self.R + self.X*1j))
        branchLine
        node1.Ys += -self.nBf2 * 1j
        node2.Ys += -self.nBf2 * 1j
        model.addBranches(branchLine)

#同上，支持修改工作状态
class LINE(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.node2 = strList[3]
        self.R = float(strList[4])
        self.X = float(strList[5])
        self.nBf2 = float(strList[6])
        self.Irated = float(strList[7])/Sb
        self.state = (int(strList[8])+int(strList[9]))
        print(f'Parsing {self.name}...')

    def apply(self, model: Model):
        if self.state == 0 or self.state == 1:
            return
        node1: Node = model.findNodeByName(self.node1)
        node2: Node = model.findNodeByName(self.node2)
        branchLine = Branch(
            self.name, node1, node2, Y=1/(self.R + self.X*1j))
        branchLine.Irated = self.Irated
        node1.Ys += -self.nBf2 * 1j
        node2.Ys += -self.nBf2 * 1j
        model.addBranches(branchLine)

#变压器模型
class THTRFO(Component):
    #解析变压器数据
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.node2 = strList[3]
        self.R = float(strList[4])
        self.X = float(strList[5])
        self.k = float(strList[6])/100
        print(f'Parsing {self.name}...')

    #在所给的模型中创建变压器对应的导纳模型
    def apply(self, model: Model):
        node1: Node = model.findNodeByName(self.node1)
        node2: Node = model.findNodeByName(self.node2)
        #支路导纳
        branchT = Branch(
            self.name, node1, node2, Y=1 / (self.R + self.X*1j) / self.k)
        #节点导纳
        node1.Ys += (self.k - 1) / (self.R + self.X*1j) / self.k
        node2.Ys += (1 - self.k) / (self.R + self.X*1j) / self.k**2
        #在模型中添加支路
        model.addBranches(branchT)

#同上，支持修改工作状态
class TRFO(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.node2 = strList[3]
        self.R = float(strList[4])
        self.X = float(strList[5])
        self.k = float(strList[6])/100
        self.Irated = float(strList[7])/Sb
        self.state = (int(strList[8])+int(strList[9]))
        print(f'Parsing {self.name}...')

    def apply(self, model: Model):
        if self.state == 0 or self.state == 1:
            return
        node1: Node = model.findNodeByName(self.node1)
        node2: Node = model.findNodeByName(self.node2)
        branchT = Branch(
            self.name, node1, node2, Y=1 / (self.R + self.X*1j) / self.k)
        branchT.Irated = self.Irated
        node1.Ys += (self.k - 1) / (self.R + self.X*1j) / self.k
        node2.Ys += (1 - self.k) / (self.R + self.X*1j) / self.k**2
        model.addBranches(branchT)

#并联功率元件模型
class THSHUNT(Component):
    #解析并联元件数据
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.G = float(strList[3])
        self.B = float(strList[4])
        print(f'Parsing {self.name}...')

    #在所给的模型中创建并联元件对应的导纳模型
    def apply(self, model: Model):
        #只更改节点自导纳
        node1 = model.findNodeByName(self.node1)
        node1.Ys += self.G + self.B * 1j

#负荷模型
class THLOAD(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.P = float(strList[3])
        self.Q = float(strList[4])
        print(f'Parsing {self.name}...')

    #将符合添加到节点上
    def apply(self, model: Model):
        #查找节点
        node1: Node = model.findNodeByName(self.node1)
        #修改节点功率
        node1.P -= self.P
        node1.Q -= self.Q
        node1.Pd += self.P
        node1.Qd += self.Q

#同上，支持修改工作状态
class LOAD2(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.P = float(strList[3])/Sb
        self.Q = float(strList[4])/Sb
        self.state = int(strList[10])
        print(f'Parsing {self.name}...')

    def apply(self, model: Model):
        if self.state == 0:
            return
        node1: Node = model.findNodeByName(self.node1)
        node1.P -= self.P
        node1.Q -= self.Q
        node1.Pd += self.P
        node1.Qd += self.Q

#发电机模型，支持修改工作状态
class GENERCV(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.P = float(strList[3])/Sb #标幺化
        self.Q = float(strList[4])/Sb
        self.V = float(strList[5])
        self.oV = self.V
        self.state = int(strList[6])
        print(f'Parsing {self.name}...')

    def apply(self, model: Model):
        if self.state == 0:
            return
        node1: Node = model.findNodeByName(self.node1)
        node1.P += self.P
        node1.Pg += self.P
        node1.V = self.V
        node1.oV = self.oV
        node1.changeType(NodeType.PV)

#同上，但是为PQ发电机模型，支持修改工作状态
class GENER(Component):
    def __init__(self, strList):
        self.type = strList[0]
        self.name = strList[1]
        self.node1 = strList[2]
        self.P = float(strList[3])/Sb
        self.Q = float(strList[4])/Sb
        self.state = int(strList[5])
        print(f'Parsing {self.name}...')

    def apply(self, model: Model):
        if self.state == 0:
            return
        node1: Node = model.findNodeByName(self.node1)
        node1.P += self.P
        node1.Pg += self.P
        node1.Q += self.Q
        node1.Qg += self.Q
        node1.changeType(NodeType.PQ)
