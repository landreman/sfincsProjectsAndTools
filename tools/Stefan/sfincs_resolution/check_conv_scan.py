#!/usr/bin/env python

from __future__ import division, print_function

import numpy as np
import os
from sfincs_simulation import Sfincs_simulation, Sfincs_input


def check_scandir(dirname):
    tol = 1e-6
    convTol = 1e-1
    to_quit = False
    subdirs = [os.path.join(dirname, f) for f in os.listdir(dirname) if os.path.isdir(os.path.join(dirname, f))]
    output =  ""
    simuls = []
    baseSimul = None
    for d in subdirs:
        s = Sfincs_simulation(d,load_geometry=False)
        simul = {}
        simul["Ntheta"] = s.Ntheta
        simul["Nzeta"] = s.Nzeta
        simul["Nxi"] = s.Nxi
        simul["Nx"] = s.Nx
        simul["d"] = d
        try:
            simul["Gamma"] = s.GammaHat
            simul["FSABFlow"] = s.FSABFlow
            simul["Q"] = s.QHat
            
        except KeyError:
            output = output + d + ": No output\n"
            if d.rsplit("/",1)[-1] == "baseCase":
                to_quit = True
        else:
            if d.rsplit("/",1)[-1] == "baseCase":
                baseSimul = simul
        simuls.append(simul)
    if to_quit:
        # if baseCase does not exist, we can't do any further checks
        return output

    for s in simuls:
        dX = None
        problems = False
        for k in ["Ntheta","Nzeta","Nxi","Nx"]:
            D = s[k]/baseSimul[k] - 1.0
            if np.fabs(D) > tol:
                dX = D
        if dX is not None:
            try:
                dGamma = s["Gamma"]/baseSimul["Gamma"] - 1.0
                dF = s["FSABFlow"]/baseSimul["FSABFlow"] - 1.0
                dQ = s["Q"]/baseSimul["Q"] - 1.0
            except KeyError:
                pass
            dGdX = dGamma/dX
            dFdX = dF/dX
            dQdX = dQ/dX
            if (np.fabs(dGdX) > convTol).any():
                output = output + s["d"] + " : Gamma :  " + str(dGdX) + " ; "
                problems = True
            if (np.fabs(dFdX) > convTol).any():
                output = output + s["d"] + " : FSABFlow :  " + str(dFdX) + " ; "
                problems = True
            if (np.fabs(dQdX) > convTol).any():
                output = output + s["d"] + " : Q :  " + str(dQdX) + " ; "
                problems = True
            if problems:
                output = output + "\n"
    return output
if __name__=="__main__":
    import sys
    argv  = sys.argv
    argc = len(argv)
    if argc>1:
        dirname = argv[1]
    else:
        print("usage:")
        print("./check_conv_scan.py <dirname of scanType1 dir>")
    print(check_scandir(dirname), end = '')
