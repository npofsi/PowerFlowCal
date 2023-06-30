import os
from powerflow.model import Profile, NodeType
root = os.getcwd()


if __name__ == '__main__':
    print(Profile(os.path.join(root,"tests\\IEEE-14.th")))

    