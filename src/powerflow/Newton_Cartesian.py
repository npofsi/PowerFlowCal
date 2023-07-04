import numpy as np
from powerflow.model import Model, NodeType
from powerflow.utils import C2P, P2C, P2Complex


class NewtonCartesian:
    times = 0

    def __init__(self, model: Model):
        self.model = model
        self.Y = self.model.deriveYMatrix()
        self.NodeCount = len(self.model.nodes)
        self.precision = 1E-4

    def solve(self):
        print("Solving...")
        flag = True
        while (flag):
            flag = self.iterate()
            # input(f"Press Enter to continue[{self.times}]...")
            pass
        self.applyPower()

    def iterate(self):
        

        print(f"\nIterating {self.times}...")
        I = self.calInjectedCurrents()
        Delta = self.calDelta(I)
        maxDelta = np.max(np.abs(Delta))
        flag = maxDelta > self.precision
        if flag:
            Jacob = self.calJacobMatrix(I)
            DV = self.calDeltaV(Jacob, Delta)
            self.applyDV2Nodes(DV)
        self.times += 1
        return flag

    def applyPower(self):
        for i, node in enumerate(self.model.nodes):
            if node.type == NodeType.PV:
                for j, node2 in enumerate(self.model.nodes):
                    if i != j:
                        self.Y[i][j] = 0
                        self.Y[j][i] = 0
        
            

    def calInjectedCurrents(self):
        I = np.zeros((self.NodeCount, 1), dtype=complex)
        for i, node in enumerate(self.model.nodes):
            injectedCurrent = 0
            subs = self.Y[i]
            for j, sub in enumerate(subs):
                injectedCurrent += sub * self.model.nodes[j].V

            I[i] = injectedCurrent
        print(f"Injected Currents:\n{I}")
        return I

    def calJacobMatrix(self, InjectedCurrents):
        I = InjectedCurrents
        MN = self.NodeCount-1
        H = np.zeros((MN, MN), dtype=float)
        N = np.zeros((MN, MN), dtype=float)
        J = np.zeros((MN, MN), dtype=float)
        L = np.zeros((MN, MN), dtype=float)
        Jacob = np.zeros(
            (MN*2, MN*2), dtype=float)

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

                Jacob[2*i][2*j] = J[i][j]
                Jacob[2*i][2*j+1] = L[i][j]
                Jacob[2*i+1][2*j] = H[i][j]
                Jacob[2*i+1][2*j+1] = N[i][j]
        with np.printoptions(linewidth=180):
            print(f"Jacob Matrix[{Jacob.shape}]:\n{Jacob}")
        return Jacob

    def calDelta(self, InjectionCurrents):
        I = InjectionCurrents
        MN = self.NodeCount-1
        Delta = np.zeros((MN*2, 1), dtype=float)

        for i in range(MN):

            if self.model.nodes[i].type == NodeType.PQ:
                e, f = self.model.nodes[i].V.real, self.model.nodes[i].V.imag
                a, b = I[i].real, I[i].imag
                Delta[2*i] = self.model.nodes[i].Q-(f*a-e*b)
                Delta[2*i+1] = self.model.nodes[i].P-(e*a+f*b)

            elif self.model.nodes[i].type == NodeType.PV:
                e, f = self.model.nodes[i].V.real, self.model.nodes[i].V.imag
                a, b = I[i].real, I[i].imag
                Delta[2*i] = np.abs(self.model.nodes[i].V)**2-(e**2+f**2)
                Delta[2*i+1] = self.model.nodes[i].P-(e*a+f*b)

        print(f"Delta[{Delta.shape}]:\n{Delta}")
        return Delta

    def calDeltaV(self, Jacob, Delta):
        DV = np.linalg.solve(Jacob, Delta)

        print(f"DeltaV[{DV.shape}]:\n{DV}")

        return DV

    def applyDV2Nodes(self, DV):
        for i in range(self.NodeCount-1):
            e, f = self.model.nodes[i].V.real, self.model.nodes[i].V.imag
            if self.model.nodes[i].type == NodeType.PQ:
                self.model.nodes[i].V = (e-DV[2*i])+(f-DV[2*i+1])*1j
                # , self.model.nodes[i].theta = C2P(
                #     e-DV[2*i], f-DV[2*i+1]*1j)
            elif self.model.nodes[i].type == NodeType.PV:
                self.model.nodes[i].V = e+(f-DV[2*i+1])*1j
                # self.model.nodes[i].V, self.model.nodes[i].theta = C2P(
                #     e, f-DV[2*i+1]*1j)

        print(f"Nodes:")
        for node in self.model.nodes:
            print(f"{node}")
