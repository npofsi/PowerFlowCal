import sys, os
import numpy as np
from powerflow.model import Model, NodeType
from powerflow.utils import C2P, P2C, P2Complex

class Report:
    def __init__(self, model: Model) -> None:
        self.model = model
        
    def listNodes(self):
        table=[]
        print(f"Nodes\ttype\t\tP\tQ\tV\t\ttheta\tSg\t\tPmax\tPmin\tQmax\tQmin\tVmax\tVmin")
        for i, node in enumerate(self.model.nodes):
            table.append('\t'.join([str(value) for value in (
                node.name,
                str(node.type),
                '%.2f'%node.P,
                '%.2f'%node.Q,
                '%.2f`%.2f'%(np.abs(node.V),np.angle(node.V)),
                '%.2f'%node.getTheta(),
                '%.2f+%.2fj'%(node.calSg().real,node.calSg().imag),
                node.Pmax,
                node.Pmin,
                node.Qmax,
                node.Qmin,
                node.Vmax,
                node.Vmin
            )]))
           
        print('\n'.join(table))
        print()
            

    def listBranches(self):
        table=[]
        print(f"Branchesfrom\tto\tY\t\tFlow\t\tLoss\t\tI\t\tIrated")
        for i, branch in enumerate(self.model.branches):
            table.append('\t'.join([str(value) for value in (
                branch.name,
                branch.node1.name,
                branch.node2.name,
                '%.2f`%.2f'%(np.abs(branch.Y),np.angle(branch.Y)),
                '%.2f`%.2f'%(np.abs(branch.Flow),np.angle(branch.Flow)),
                '%.3f`%.2f'%(np.abs(branch.Loss),np.angle(branch.Loss)),
                '%.2f`%.2f'%(np.abs(branch.I),np.angle(branch.I)),
                '%.2f'%branch.Irated
            )]))
        print('\n'.join(table))
        print()

    def showTotalLoss(self):
        print(f"Total Loss: {self.model.loss}[{'%.2f`%.2f'%(np.abs(self.model.loss),np.angle(self.model.loss))}]")