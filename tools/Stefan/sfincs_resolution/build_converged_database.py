#!/usr/bin/env python
from __future__ import division,print_function

from filename2ref import filename2ref
from sfincs_simulation import Sfincs_simulation

import matplotlib.pyplot as plt
import f90nml
import numpy as np
import os

converged_simuldirs = []
csfn = os.path.abspath(__file__).rsplit("/",1)[0] + "/converged_simulations.dat"
with open(csfn,'r') as f:
    ls = f.readlines()
    for l in ls:
        if len(l[:-1])>0:
            converged_simuldirs.append(l[:-1])

print(converged_simuldirs)

simuls = [Sfincs_simulation(d) for d in converged_simuldirs]

cols = [s.collisionality for s in simuls]
Nthetas = [s.Ntheta for s in simuls]
Nzetas = [s.Nzeta for s in simuls]
Nxs = [s.Nx for s in simuls]
Nxis = [s.Nxi for s in simuls]
eqnames = [s.input.equilibrium_name for s in simuls]

eqnmlfn = os.path.abspath(__file__).rsplit("/",1)[0] + "/equilibria.namelist"
eqnml = f90nml.read(eqnmlfn)["w7x"]
eqtypes = [eqnml[filename2ref(n)] for n in eqnames]

# for each eqtype, create a database
types = set(eqtypes)
coltol = 1e-6 # collisionalities within this distance are treated as being the same
use = "max" # use the highest resolution when there is a collision in collisionality


for t in types:
    print(t)
    _cols = [x for x,e in zip(cols,eqtypes) if e==t]
    _colmins = [np.min(c) for c in _cols]
    _colmaxs = [np.max(c) for c in _cols]
    
    _Nthetas = [x for x,e in zip(Nthetas,eqtypes) if e==t]
    _Nzetas = [x for x,e in zip(Nzetas,eqtypes) if e==t]
    _Nxs = [x for x,e in zip(Nxs,eqtypes) if e==t]
    _Nxis = [x for x,e in zip(Nxis,eqtypes) if e==t]
    tosort = zip(_cols,_Nthetas,_Nzetas,_Nxs,_Nxis)
    print(tosort)
    sortedarray = sorted(zip(_colmins,_Nthetas,_Nzetas,_Nxis))
    _colmins = [x[0] for x in sortedarray]
    _Nthetas = [x[1] for x in sortedarray]
    _Nzetas = [x[2] for x in sortedarray]
    _Nxis = [x[3] for x in sortedarray]
    sortedarray = sorted(zip(_colmaxs,_Nxs))
    _colmaxs = [x[0] for x in sortedarray]
    _Nxs = [x[1] for x in sortedarray]
    Nsimuls = len(_colmins)

    newcolmins = [_colmins[0]]
    newNthetas = [_Nthetas[0]]
    newNzetas = [_Nzetas[0]]
    newNxis = [_Nxis[0]]
    newcolmaxs = [_colmaxs[0]]
    newNxs = [_Nxs[0]]
    
    for i in range(1,Nsimuls):
        if np.fabs(_colmins[i-1] - _colmins[i]) < coltol:
            newNthetas[-1] = eval(use + "([_Nthetas[i-1],_Nthetas[i]])")
            newNzetas[-1] = eval(use + "([_Nzetas[i-1],_Nzetas[i]])")
            newNxis[-1] = eval(use + "([_Nxis[i-1],_Nxis[i]])")
            # do not update collisionality to prevent infinitely lumping together
            # points finely spaced in collisionality
        else:
            newNthetas.append(_Nthetas[i])
            newNzetas.append(_Nzetas[i])
            newNxis.append(_Nxis[i])
            newcolmins.append(_colmins[i])
            
        if np.fabs(_colmaxs[i-1] - _colmaxs[i]) < coltol:
            newNxs[-1] = eval(use + "([_Nxs[i-1],_Nxs[i]])")
        else:
            newNxs.append(_Nxs[i])
            newcolmaxs.append(_colmaxs[i])
            
    
    with open(t + ".dat",'w') as f:
        f.write("# rows: col(min), Ntheta, Nzeta, Nxi, col(max), Nx\n")
        f.write(",".join(["{:.6e}".format(x) for x in newcolmins]) + "\n")
        f.write(",".join(["{:.6e}".format(x) for x in newNthetas]) + "\n")
        f.write(",".join(["{:.6e}".format(x) for x in newNzetas]) + "\n")
        f.write(",".join(["{:.6e}".format(x) for x in newNxis]) + "\n")
        f.write(",".join(["{:.6e}".format(x) for x in newcolmaxs]) + "\n")
        f.write(",".join(["{:.6e}".format(x) for x in newNxs]))
        

