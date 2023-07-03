import numpy as np
from powerflow.model import Model


class Newton:
    times = 0

    def __init__(self, model: Model):
        self.model = model
        self.Y = self.model.deriveYMatrix()
        self.NodeCount = len(self.model.nodes)

    def solve(self):
        print("Solving...")
        precision = 1E-4

        self.calJacobMatrix()

        pass

    def iterate(self):
        self.times += 1
        print(f"Iterating {self.times}...")
        pass

    def calInjectedCurrents(self):
        I = np.zeros((self.NodeCount, 1), dtype=complex)
        for i, node in enumerate(self.model.nodes):
            injectedCurrent = 0
            subs = self.Y[i]
            for j, sub in enumerate(subs):
                injectedCurrent += sub*self.model.nodes[j].V

            I[i] = injectedCurrent
        print(f"Injected Currents:\n{I}")
        return I

    def calJacobMatrix(self):
        I = self.calInjectedCurrents()
        MN = self.NodeCount-1
        H = np.zeros((MN, MN), dtype=complex)
        N = np.zeros((MN, MN), dtype=complex)
        J = np.zeros((MN, MN), dtype=complex)
        L = np.zeros((MN, MN), dtype=complex)
        Jacob = np.zeros(
            (MN*2, MN*2), dtype=complex)
        for i in range(MN):
            for j in range(MN):
                G, B = self.Y[i][j].real, self.Y[i][j].imag
                e, f = self.model.nodes[i].V.real, self.model.nodes[j].V.imag
                H[i][j] = -(G*e+B*f)
                N[i][j] = (B*e-G*f)

                if 
                J[i][j] = N[i][j]
                L[i][j] = -H[i][j]
                if i == j:
                    a, b = I[i].real, I[i].imag
                    H[i][j] += -a
                    N[i][j] += -b
                    J[i][j] += b
                    L[i][j] += -a
                
                Jacob[2*i][2*j]+=H[i][j]
                Jacob[2*i][2*j+1]+=N[i][j]
                Jacob[2*i+1][2*j]+=J[i][j]
                Jacob[2*i+1][2*j+1]+=L[i][j]

        print(f"Jacob Matrix:\n{Jacob}")
        return Jacob
