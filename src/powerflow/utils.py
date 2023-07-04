import numpy as np


def P2C(V, theta):
    return V * np.exp(theta * 1j).real, V * np.exp(theta * 1j).imag

def C2P(Real, Imag):
    return np.abs(Real + Imag * 1j), np.angle(Real + Imag * 1j)

def P2Complex(V, theta):
    return V * np.exp(theta * 1j)

def Complex2P(Complex):
    return np.abs(Complex), np.angle(Complex)