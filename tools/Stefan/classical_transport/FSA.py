from numpy import sum
from numba import jit

def FSA_2(X,B):
    # assumes BOOZER coordinates
    return sum(X/B**2)/sum(1/B**2)

def FSA(X,B):
    # assumes BOOZER coordinates
    # sums the two last dimensions
    return sum(X/B**2,axis=(-1,-2))/sum(1/B**2)

@jit
def jit_FSA(X,B):
    # assumes BOOZER coordinates
    return sum(X/B**2)/sum(1/B**2)


if __name__=="__main__":
    import numpy as np
    import timeit

    Nx = 91
    Ny = 99
    x=np.linspace(0,1,Nx)
    y=np.linspace(0,1,Ny)

    Y,X = np.meshgrid(y,x)
    B = 2 + np.sin(X)
    N=2
    A = np.zeros((N,Nx,Ny))
    A[:] = np.cos(Y+X)

    def test_FSA(X,B):
        def _FSA():
            return FSA(X,B)
        return _FSA

    def test_jit_FSA(X,B):
        def _jit_FSA():
            return jit_FSA(X,B)
        return _jit_FSA

    t=timeit.Timer(test_FSA(A,B))
    print min(t.repeat(100,1))
    
    t=timeit.Timer(test_FSA(A[0],B))
    print N*min(t.repeat(100,1))

    tjit=timeit.Timer(test_jit_FSA(A[0],B))
    print N*min(tjit.repeat(100,1))
