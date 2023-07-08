import numpy as np
from powerflow.model import Model, NodeType
from powerflow.utils import C2P, P2C, P2Complex

#牛顿迭代法，直角坐标法
class NewtonCartesian:
    times = 0

    #输入模型，获取节点导纳矩阵
    def __init__(self, model: Model):
        self.model = model
        self.Y = self.model.deriveYMatrix()
        self.NodeCount = len(self.model.nodes)
        self.precision = 1E-6

    #求解
    def solve(self):
        print("Solving...")
        # self.initQ()

        flag = True #标记变量，标记是否继续迭代

        #求解
        while (flag):
            flag = self.iterate()#迭代一次
            # input(f"Press Enter to continue[{self.times}]...") #调试用
            pass

        #计算节点注入功率
        self.applyPower()
        #计算支路功率
        self.calBranchesFlow()
        print(f"Nodes:")
        for node in self.model.nodes:
            print(f"{node}")

    #一次迭代
    def iterate(self):
        print(f"\nIterating {self.times}...")
        I = self.calInjectedCurrents() #计算注入电流
        Delta = self.calDelta(I) #计算ΔP,ΔQ，ΔV

        #获取Δ的最大值判断是否继续迭代
        maxDelta = np.max(np.abs(Delta)) 
        flag = maxDelta > self.precision

        if flag:
            #迭代，计算Jacobi矩阵，从而求解ΔV
            Jacob = self.calJacobMatrix(I)
            DV = self.calDeltaV(Jacob, Delta)
            #将ΔV应用于节点
            self.applyDV2Nodes(DV)

        #计数用
        self.times += 1
        return flag

    # def initQ(self):
    #     Qmax = -99999
    #     Qmin = 99999
    #     for i, node in enumerate(self.model.nodes):
    #         if node.type == NodeType.PQ:
    #             Qmax = max(Qmax, node.Q)
    #             Qmin = min(Qmin, node.Q)
    #     for i, node in enumerate(self.model.nodes):
    #         if node.type == NodeType.PV:
    #             node.Q = (Qmax + Qmin)/2

    #计算支路功率，支路电流，支路损耗
    def calBranchesFlow(self):
        for i, branch in enumerate(self.model.branches):
            branch.I = (branch.node1.V - branch.node2.V) * branch.Y #支路电流

            branch.Loss = np.abs(branch.I)**2 / branch.Y.conjugate() #支路损耗
            self.model.loss += branch.Loss #将支路损耗加到总损耗上

            branch.Flow = branch.node1.V * branch.I.conjugate() #支路功率

    #计算节点缺失的功率
    def applyPower(self):
        for i, node in enumerate(self.model.nodes):
            S = 0.+0.j
            if node.type == NodeType.PV:
                #PV节点，计算缺失的Q
                for j, node2 in enumerate(self.model.nodes):
                    S+= (node.V*node2.V)*self.Y[i][j]
                node.Q = S.imag
                node.V = P2Complex(node.oV, node.getTheta())
            if node.type == NodeType.Slack:
                #平衡节点，计算缺失的P,Q
                for j, node2 in enumerate(self.model.nodes):
                    S+= (node.V*node2.V)*self.Y[i][j]
                node.P = S.real
                node.Q = S.imag
            
            
    #计算每个节点注入电流
    def calInjectedCurrents(self):
        I = np.zeros((self.NodeCount, 1), dtype=complex)

        #遍历每个节点
        for i, node in enumerate(self.model.nodes):
            injectedCurrent = 0
            subs = self.Y[i]
            for j, sub in enumerate(subs):
                injectedCurrent += sub * self.model.nodes[j].V

            I[i] = injectedCurrent
        print(f"Injected Currents:\n{I}")
        return I

    #计算Jacobi矩阵，需要提供注入电流
    def calJacobMatrix(self, InjectedCurrents):
        I = InjectedCurrents
        MN = self.NodeCount-1 #减去平衡节点
        H = np.zeros((MN, MN), dtype=float)
        N = np.zeros((MN, MN), dtype=float)
        J = np.zeros((MN, MN), dtype=float)
        L = np.zeros((MN, MN), dtype=float)
        Jacob = np.zeros(
            (MN*2, MN*2), dtype=float)

        #遍历矩阵的每一项
        for i in range(MN):
            for j in range(MN):

                G, B = self.Y[i][j].real, self.Y[i][j].imag
                e, f = self.model.nodes[i].V.real, self.model.nodes[i].V.imag

                if(i==1 and j==3):
                    print(f"Y[{i}][{j}]:{self.Y[i][j]}")
                    print(f"G:{G}, B:{B}")
                    print(f"e:{e}, f:{f}")
                H[i][j] = -(G*e+B*f)
                N[i][j] = (B*e-G*f)
                if self.model.nodes[i].type == NodeType.PQ:
                    J[i][j] = N[i][j]
                    L[i][j] = -H[i][j]
                elif self.model.nodes[i].type == NodeType.PV:
                    J[i][j] = 0
                    L[i][j] = 0

                if i == j:
                    a, b = I[i].real, I[i].imag
                    H[i][j] += -a
                    N[i][j] += -b
                    if self.model.nodes[i].type == NodeType.PQ:
                        J[i][j] += b
                        L[i][j] += -a
                    elif self.model.nodes[i].type == NodeType.PV:
                        J[i][j] = -2*e
                        L[i][j] = -2*f

                #将四个子矩阵 交叉 填入Jacob矩阵
                Jacob[2*i][2*j] = J[i][j]
                Jacob[2*i][2*j+1] = L[i][j]
                Jacob[2*i+1][2*j] = H[i][j]
                Jacob[2*i+1][2*j+1] = N[i][j]
        with np.printoptions(linewidth=180):
            print(f"Jacob Matrix[{Jacob.shape}]:\n{Jacob}")
        return Jacob

    #计算DeltaP，DeltaQ，DeltaV^2，最后总结为一个向量，需要提供注入电流
    def calDelta(self, InjectionCurrents):
        I = InjectionCurrents
        MN = self.NodeCount-1
        Delta = np.zeros((MN*2, 1), dtype=float)

        for i in range(MN):
            
            #PQ节点,计算DeltaP，DeltaQ
            if self.model.nodes[i].type == NodeType.PQ:
                e, f = self.model.nodes[i].V.real, self.model.nodes[i].V.imag
                a, b = I[i].real, I[i].imag
                Delta[2*i] = self.model.nodes[i].Q-(f*a-e*b)
                Delta[2*i+1] = self.model.nodes[i].P-(e*a+f*b)

            #PV节点，计算DeltaP，DeltaV^2
            elif self.model.nodes[i].type == NodeType.PV:
                e, f = self.model.nodes[i].V.real, self.model.nodes[i].V.imag
                a, b = I[i].real, I[i].imag
                Delta[2*i] = np.abs(self.model.nodes[i].V)**2-(e**2+f**2)
                Delta[2*i+1] = self.model.nodes[i].P-(e*a+f*b)

        print(f"Delta[{Delta.shape}]:\n{Delta}")
        return Delta

    #计算DeltaV，需要提供Jacobi矩阵和Delta，解方程
    def calDeltaV(self, Jacob, Delta):
        DV = np.linalg.solve(Jacob, Delta) #使用 numpy 所携带的线性方程求解器求解

        print(f"DeltaV[{DV.shape}]:\n{DV}")

        return DV

    #将计算得到的DeltaV应用到节点上
    def applyDV2Nodes(self, DV):
        for i in range(self.NodeCount-1):
            e, f = self.model.nodes[i].V.real, self.model.nodes[i].V.imag
            if self.model.nodes[i].type == NodeType.PQ:
                self.model.nodes[i].V = (e-DV[2*i])+(f-DV[2*i+1])*1j
            elif self.model.nodes[i].type == NodeType.PV:
                self.model.nodes[i].V = (e-DV[2*i])+(f-DV[2*i+1])*1j

        