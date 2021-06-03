#!/usr/bin/env python

from sfincs_simulation import Sfincs_simulation, Sfincs_input
from filename2ref import filename2ref
from get_resolution import get_database, get_resolution

import f90nml

def set_resolution(dirname):
    s = Sfincs_simulation(dirname)
    database  = get_database(s.input.equilibrium_name)
    Ntheta,Nzeta,Nxi,Nx = get_resolution(database,s.collisionality)
    print("simul '{}' has geometry '{}'".format(dirname,database))
    print("Ntheta: {}->{}".format(s.Ntheta,Ntheta))
    print("Nzeta: {}->{}".format(s.Nzeta,Nzeta))
    print("Nxi: {}->{}".format(s.Nxi,Nxi))
    print("Nx: {}->{}".format(s.Nx,Nx))
    s.input.changevar("resolutionparameters","ntheta",Ntheta)
    s.input.changevar("resolutionparameters","nzeta",Nzeta)
    s.input.changevar("resolutionparameters","nxi",Nxi)
    s.input.changevar("resolutionparameters","nx",Nx)
    
if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    if argc > 1:
        dirname = argv[1]
        set_resolution(dirname)
    else:
        print("usage: ./set_resolution.py <sfincs simulation dirname>")
 
