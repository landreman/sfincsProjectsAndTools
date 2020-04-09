from __future__ import division, print_function

from sfincs_simulation import Sfincs_simulation

from scipy.interpolate import interp1d
import numpy as np

import os



def interpolate_resolution(reff,interpolators):
    resolution = []
    for i,interpolator in enumerate(interpolators):
        N = int(np.ceil(interpolator(reff)))
        if (i == 0) or (i == 1):
            if N%2 == 0:
                N = N + 1
        resolution.append(N)
    return resolution

def read_Er_interpolators(filename):
    reff = []
    Er = []
    interval = []
    with open(filename) as f:
        for l in f.readlines():
            ls = l.strip()
            if len(ls) == 0:
                continue
            if ls[0] == '#':
                continue
            _reff, _Er, _interval = ls.split()
            reff.append(float(_reff))
            Er.append(float(_Er))
            interval.append(int(_interval))
    if len(Er) == 0:
        return [], [], []
    Nintervals = max(interval)+1

    reffs = [[] for i in range(Nintervals)]
    Ers = [[] for i in range(Nintervals)]
    all_reffs = [[] for i in range(Nintervals)]
    for i in range(Nintervals):
        for _reff, _Er, _interval in zip(reff,Er,interval):
            if _interval == i:
                print(_reff)
                all_reffs[i].append(_reff)
                if not np.isnan(_Er):
                    Ers[i].append(_Er)
                    reffs[i].append(_reff)
    interpolators = []
    minreff = []
    maxreff = []
    for all_reff,reff,Er in zip(all_reffs,reffs,Ers):
        interpolators.append(interp1d(reff,Er,fill_value="extrapolate"))
        minreff.append(min(all_reff))
        maxreff.append(max(all_reff))
        
    return interpolators, minreff, maxreff 


def read_resolution_interpolators(filename = "resolution.dat"):
    reff = []
    Ntheta = []
    Nzeta = []
    Nxi = []
    Nx = []
    NL = []    
    with open(filename) as f:
        for l in f.readlines():
            ls = l.strip()
            if len(ls) == 0:
                continue
            if ls[0] == '#':
                continue
            _reff, _Ntheta, _Nzeta, _Nxi, _Nx, _NL = ls.split()
            reff.append(float(_reff))
            Ntheta.append(int(_Ntheta))
            Nzeta.append(int(_Nzeta))
            Nxi.append(int(_Nxi))
            Nx.append(int(_Nx))
            NL.append(int(_NL))
    return interp1d(reff,Ntheta,fill_value="extrapolate"), interp1d(reff,Nzeta,fill_value="extrapolate"), interp1d(reff,Nxi,fill_value="extrapolate"), interp1d(reff,Nx,fill_value="extrapolate"), interp1d(reff,NL,fill_value="extrapolate")
         
            
def update_simulation(simul,resolution_interpolators,i_interpolator, e_interpolator, Ermargin = 0.2, NErs = 3):
    dPhiHatdrhat = simul.input.get_value_from_input_or_defaults("physicsParameters","dPhiHatdrHat")
    if dPhiHatdrhat>0:
        root = "i"
        Er_interpolator = i_interpolator
    else:
        root = "e"
        Er_interpolator = e_interpolator
    reff = simul.input.get_value_from_input_or_defaults("geometryParameters","rHat_wish")
    Ntheta, Nzeta, Nxi, Nx, NL = interpolate_resolution(reff,resolution_interpolators)
    reffmin = Er_interpolator[1]
    reffmax = Er_interpolator[2]
    Er_interpolator = Er_interpolator[0]
    for i,(rmin, rmax) in enumerate(zip(reffmin,reffmax)):
        if (reff >= rmin) and (reff <= rmax):
            Er = Er_interpolator[i](reff)
            break
    tol = Ermargin
    if Er > tol:
        Ermin = Er * (1.0 - Ermargin)
        Ermax = Er * (1.0 + Ermargin)
    elif Er < -tol:
        Ermin = Er * (1.0 + Ermargin)
        Ermax = Er * (1.0 - Ermargin)
    else:
        # we are very close to zero, do not use relative
        Ermin = -Ermargin
        Ermax = Ermargin
    print(Ermax)
    print(Ermin)
    simul.input.changessvar("dPhiHatdrHatMin",-Ermax)
    simul.input.changessvar("dPhiHatdrHatMax",-Ermin)
    simul.input.changessvar("NErs",NErs)
    simul.input.changessvar("scanType",2)
    
# main starts here    
this_dir = "."
subdirs = [name for name in os.listdir(this_dir) if os.path.isdir(os.path.join(this_dir, name))]

resolution_interpolators = read_resolution_interpolators()
i_interpolator = read_Er_interpolators("i-root.dat")
e_interpolator = read_Er_interpolators("e-root.dat")

for subdir in subdirs:
    print(subdir)
    
    simul = Sfincs_simulation(subdir,load_geometry=False)
    update_simulation(simul,resolution_interpolators,i_interpolator, e_interpolator)
