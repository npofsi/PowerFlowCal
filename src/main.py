import numpy as np
import os
import sys
from powerflow.Newton_Cartesian import NewtonCartesian
from powerflow.Newton_Polar import NewtonPolar
from powerflow.model import Model, Profile, NodeType
root = os.getcwd()


if __name__ == '__main__':
    original_stdout = sys.stdout
    sys.stdout = open('output.txt', 'w')
    profile = Profile(os.path.join(root, "tests\\IEEE-30.th"))
    print(profile)
    model = Model()
    model.compose(profile)

    print(f"Node Count:{len(model.nodes)}, Branch Count:{len(model.branches)}")

    model.deriveYMatrix()  # nxn matrix

    model.printTopology()

    cal = NewtonCartesian(model)

    cal.solve()
