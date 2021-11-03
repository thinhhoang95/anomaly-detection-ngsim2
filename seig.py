import numpy as np
from numpy import linalg as npla

def eig(A):
    eigenValues, eigenVectors = npla.eig(A)
    idx = np.argsort(eigenValues)[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    return (eigenValues, eigenVectors)