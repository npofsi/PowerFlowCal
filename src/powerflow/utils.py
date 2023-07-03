import numpy as np


def C2P(V, theta):
    return V * np.exp(theta * 1j)

def P2C(S):
    return np.abs(S), np.angle(S)