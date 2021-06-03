from __future__ import division, print_function
import numpy as np
import f90nml
import os

from filename2ref import filename2ref

def get_database(equilibrium_filename):
    fn=os.path.abspath(__file__).rsplit("/",1)[0] + "/equilibria.namelist"
    eqnml = f90nml.read(fn)["w7x"]
    eqtype = eqnml[filename2ref(equilibrium_filename)]
    return eqtype + ".dat"

def get_resolution(database,cols):
    mincol = min(cols)
    maxcol = max(cols)

    try:
        i = database.index('/')
    except ValueError:
        database = os.path.abspath(__file__).rsplit("/",1)[0] + "/" + database

    with open(database,'r') as f:
        ls = f.readlines()
    mincols = [float(x) for x in ls[1][:-1].split(',')]
    Nthetas = [int(float(x)) for x in ls[2][:-1].split(',')]
    Nzetas = [int(float(x)) for x in ls[3][:-1].split(',')]
    Nxis = [int(float(x)) for x in ls[4][:-1].split(',')]
    maxcols = [float(x) for x in ls[5][:-1].split(',')]
    Nxs = [int(float(x)) for x in ls[6].split(',')]

    Ntheta = np.interp(mincol,mincols,Nthetas)
    Nzeta = np.interp(mincol,mincols,Nzetas)
    Nxi = np.interp(mincol,mincols,Nxis)
    Nx = np.interp(maxcol,maxcols,Nxs)
    
    Ntheta = np.round(Ntheta)
    Nzeta = np.round(Nzeta)
    Nxi = np.round(Nxi)
    Nx = np.round(Nx)
    

    if Ntheta % 2 == 0:
        Ntheta = Ntheta + 1
    if Nzeta % 2 == 0:
        Nzeta = Nzeta + 1
    
    return (int(Ntheta),int(Nzeta),int(Nxi),int(Nx))


if __name__ == "__main__":
    Ntheta,Nzeta,Nxi,Nx = get_resolution("EIM.dat",[0.3,1.0])
    print(Ntheta)
    print(Nzeta)
    print(Nxi)
    print(Nx)
