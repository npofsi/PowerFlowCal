import numpy as np
import os
import sys
from powerflow.Newton_Cartesian import NewtonCartesian
from powerflow.Newton_Polar import NewtonPolar
from powerflow.model import Model, Profile, NodeType
from powerflow.report import Report
root = os.getcwd()

#主程序
if __name__ == '__main__':
    original_stdout = sys.stdout
    # sys.stdout = open('output.txt', 'w')

    # 读取数据文件
    profile = Profile(os.path.join(root, "tests\\IEEE-14.th"))
    print(profile)

    # 生成一个空的模型
    model = Model()

    # 从数据文件中读取信息，创建元件模型，并依赖元件模型，生成节点导纳模型
    model.compose(profile)

    print(f"Node Count:{len(model.nodes)}, Branch Count:{len(model.branches)}")

    # 从节点导纳模型生成节点导纳矩阵
    model.deriveYMatrix()  # nxn matrix

    #输出节点导纳模型信息
    model.printTopology()

    #创建计算模型使用的对象，此处使用的是极坐标系下的牛顿拉夫逊法，也可以使用直角坐标系下的牛顿拉夫逊法，即使用NewtonCartesian类
    cal = NewtonPolar(model)

    #计算求解模型
    cal.solve()

    #创建输出结果报告对象
    report = Report(model)

    #输出结果报告
    #
    #输出节点信息表
    report.listNodes()
    #输出支路信息表
    report.listBranches()
    #输出网络总损耗
    report.showTotalLoss()
