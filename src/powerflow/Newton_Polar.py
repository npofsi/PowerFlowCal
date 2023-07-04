import numpy as np
from powerflow.model import Model, NodeType
from powerflow.utils import C2P, P2C, P2Complex


class NewtonPolar:
    def __init__(self, model: Model):
        self.model = model
        self.Y = self.model.deriveYMatrix()
        self.NodeCount = len(self.model.nodes)
        self.precision = 1E-4

    def solve(self):
        NodeData = self.genNodeData()
        NodeData = self.cal(NodeData)
        with np.printoptions(linewidth=180):
          print(f"节点\t类型\tP有功\tQ无功\t~\t~\tV电压幅值\ttheta电压相位")
          print(f"{NodeData}")
        self.applyResult(NodeData)

    def genNodeData(self):
        # 节点	类型	发电机有功	发电机无功	负荷有功	负荷无功	电压幅值	电压相位
        self.node = np.zeros((self.NodeCount, 8))
        for i in range(self.NodeCount):
            self.node[i, 0] = i+1
            self.node[i, 1] = 4 - self.model.nodes[i].type.value
            self.node[i, 2] = self.model.nodes[i].P
            self.node[i, 3] = self.model.nodes[i].Q
            self.node[i, 4] = 0  # self.model.nodes[i].P
            self.node[i, 5] = 0  # self.model.nodes[i].Q
            self.node[i, 6] = self.model.nodes[i].V
            self.node[i, 7] = self.model.nodes[i].getTheta()
        # self.node = np.flip(self.node, axis=0)
        return self.node

    def applyResult(self, node):
        # self.node = np.flip(node, axis=0)
        for i in range(self.NodeCount):
            self.model.nodes[i].P = node[i, 2]
            self.model.nodes[i].Q = node[i, 3]
            self.model.nodes[i].V = node[i, 6]
            self.model.nodes[i].setTheta(node[i, 7])
        print(f"Nodes:")
        for node in self.model.nodes:
            print(f"{node}")

    def cal(self, node):

        Y = self.model.deriveYMatrix()
        YY = Y.copy()
        ####### YY需要调整，根据PQ、PV、平衡顺序,先 PQ、再PV、再平衡，即前9个是pq，中间是4个pv，最后是一个平衡。node顺序和Y顺序都要改！！！！#######
        YR = np.real(YY)
        YI = np.imag(YY)

        n = len(self.model.nodes)

        nPQ = len([node for node in self.model.nodes if node.type == NodeType.PQ])
        err = self.precision
        snet = np.zeros((1, n), dtype=complex)  # 每个节点发电机与负荷的净注入功率
        for i in range(n):
            snet[0, i] = node[i, 2] + node[i, 3]*1j - node[i, 4] - \
                node[i, 5]*1j  # （发电机注入功率P、Q减去节点输出功率P、Q）
        # 功率的给定值
        f = 0  # 迭代标志
        nit = 1  # 当前迭代次数
        nitmax = 100
        # 开始迭代

        while (f == 0 and nit < nitmax):
            # 计算P-Q不平衡量
            str = np.zeros((1, n), dtype=complex)  # 各节点净注入功率
            for i in range(n):
                psum = 0
                qsum = 0
                for j in range(n):
                    # 节点电压方程（极坐标形式）
                    ################## 这里得改，公式如下####################
                    # Pi+=Vi*Vj*YRij*cos(thetai-thetaj) + Vi*Vj*YIij*sin(thetai-thetaj)
                    psum += node[i, 6] * node[j, 6] * YR[i, j] * np.cos(
                        node[i, 7] - node[j, 7]) + node[i, 6] * node[j, 6] * YI[i, j] * np.sin(node[i, 7] - node[j, 7])
                    # 计算每个节点传输有功
                    # Qi+=Vi*Vj*YRij*sin(thetai-thetaj) - Vi*Vj*YIij*cos(thetai-thetaj)
                    qsum += node[i, 6] * node[j, 6] * YR[i, j] * np.sin(
                        node[i, 7] - node[j, 7]) - node[i, 6] * node[j, 6] * YI[i, j] * np.cos(node[i, 7] - node[j, 7])
                    # 计算每个节点传输无功
                    str[0, i] = psum + qsum*1j

            # 记录每个节点传输的功率
            P = np.zeros((1, n))
            Q = np.zeros((1, n))
            for i in range(n):
                P[0, i] = np.real(str[0, i])  # 提取各个节点的有功功率并存入数组
                Q[0, i] = np.imag(str[0, i])  # 提取各个节点的无功功率并存入数组
                DS = snet - str
            # 得到PQ与给定的偏差

            # 判断是否满足误差要求
            DP = np.zeros((1, n-1))
            DQ = np.zeros((1, nPQ))
            # 统计前n-1个节点的有功不平衡量绝对值
            for i in range(n-1):
                DP[0, i] = np.abs(np.real(DS[0, i]))
            # 统计PQ节点的无功不平衡量绝对值
            for i in range(nPQ):
                DQ[0, i] = np.abs(np.imag(DS[0, i]))
            A = np.max(DP)
            B = np.max(DQ)
            MAX = max(A, B)
            if MAX < err:
                f = 1  # 若误差满足要求，则标志位置1，表示停止迭代
            # 判断结束，若不满足，继续迭代

            DP = np.zeros((1, n-1))
            DQ = np.zeros((1, nPQ))
            # 统计前n-1个节点的有功不平衡量
            for i in range(n-1):
                DP[0, i] = np.real(DS[0, i])
            # 统计PQ节点的无功不平衡量
            for i in range(nPQ):
                DQ[0, i] = np.imag(DS[0, i])
            # 形成PQ不平衡量
            delt_PQ = np.concatenate([DP.T, DQ.T], axis=0)

            # 形成雅可比矩阵
            H = np.zeros((n-1, n-1))
            L = np.zeros((nPQ, nPQ))
            N = np.zeros((n-1, nPQ))
            J = np.zeros((nPQ, n-1))

            # 形成矩阵H
            for i in range(n-1):
                for j in range(n-1):
                    if i != j:
                        ################## 这几个需要改，公式如下####################
                        # H[i, j]= -Vi*Vj*( YR[i, j] * np.sin(thetai-thetaj) - YI[i, j] * np.cos(thetai-thetaj) )
                        H[i, j] = -node[i, 6] * node[j, 6] * (YR[i, j] * np.sin(
                            node[i, 7] - node[j, 7]) - YI[i, j] * np.cos(node[i, 7] - node[j, 7]))
                    else:
                        #  H[i, j]=(Vi**2)*YI[i, i] + Q[0, i]
                        H[i, j] = node[i, 6]**2 * YI[i, i] + Q[0, i]
            # 形成矩阵N
            for i in range(n-1):
                for j in range(nPQ):
                    if i != j:
                      # N[i, j]= -Vi*Vj*( YR[i, j] * np.COS(thetai-thetaj) + YI[i, j] * np.sin(thetai-thetaj) )
                        N[i, j] = -node[i, 6] * node[j, 6] * (YR[i, j] * np.cos(
                            node[i, 7] - node[j, 7]) + YI[i, j] * np.sin(node[i, 7] - node[j, 7]))
                    else:
                      #  N[i, j]= -(Vi**2)*YR[i, i] - P[0, i]
                        N[i, j] = -node[i, 6]**2 * YR[i, i] - P[0, i]
            # 形成矩阵J
            for i in range(nPQ):
                for j in range(n-1):
                    if i != j:
                        #   J[i, j]=Vi*Vj*( YR[i, j] * np.COS(thetai-thetaj) + YI[i, j] * np.SIN(thetai-thetaj) )
                        J[i, j] = node[i, 6] * node[j, 6] * (YR[i, j] * np.cos(
                            node[i, 7] - node[j, 7]) + YI[i, j] * np.sin(node[i, 7] - node[j, 7]))
                    else:
                        #  J= (Vi**2)*YR[i, i] - P[0, i]
                        J[i, j] = node[i, 6]**2 * YR[i, i] - P[0, i]
            # 形成矩阵L
            for i in range(nPQ):
                for j in range(nPQ):
                    if i != j:
                        # L[i, j]= -Vi*Vj*( YR[i, j] * np.sin(thetai-thetaj) - YI[i, j] * np.cos(thetai-thetaj) )
                        L[i, j] = -node[i, 6] * node[j, 6] * (YR[i, j] * np.sin(
                            node[i, 7] - node[j, 7]) - YI[i, j] * np.cos(node[i, 7] - node[j, 7]))
                    else:
                        # L[i, j] = (Vi**2)*YI[i, i] - Q[0, i]
                        L[i, j] = node[i, 6]**2 * YI[i, i] - Q[0, i]
            JX = np.concatenate(
                [np.concatenate([H, N], axis=1), np.concatenate([J, L], axis=1)], axis=0)
            # 求电压和相角的修正值
            Ud2 = np.zeros((nPQ, nPQ))
            for i in range(nPQ):
                for j in range(nPQ):
                    if i == j:
                        Ud2[i, j] = node[i, 6]  # PQ节点的电压幅值
            delt = -np.linalg.inv(JX) @ delt_PQ
            # 分别求幅值和相位的修正量
            delt_p = delt[0:n-1, 0]  # 相位
            delt_a = Ud2 @ (delt[n-1:n+nPQ, 0].reshape(-1, 1))  # 幅值

            # 求得电压和相位的修正量

            node[0:n-1, 7] += delt_p
            node[0:nPQ, 6] += np.reshape(delt_a, -1)
            nit = nit + 1
            # 迭代结束，得到每个节点电压幅值和相位的结果

            # 开始计算发电机功率
            # 负荷功率均为给定值
            # str中记录了结果中每个节点注入的功率
            for i in range(n):
                if node[i, 1] == 2:  # PV节点，需求解注入的无功
                    node[i, 3] = np.imag(str[0, i]) + node[i, 5]
                elif node[i, 1] == 1:  # 平衡节点，需求解注入的有功无功
                    node[i, 2] = np.real(str[0, i]) + node[i, 4]
                    node[i, 3] = np.imag(str[0, i]) + node[i, 5]

        return node