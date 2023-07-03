import numpy as np
import os
from powerflow.cal import Newton
from powerflow.model import Model, Profile, NodeType
root = os.getcwd()


if __name__ == '__main__':


    profile = Profile(os.path.join(root, "tests\\IEEE-14.th"))
    print(profile)
    model = Model()
    model.compose(profile)

    print(f"Node Count:{len(model.nodes)}, Branch Count:{len(model.branches)}")

    model.deriveYMatrix()  # nxn matrix

    model.printTopology()

    cal = Newton(model)

    cal.solve()
