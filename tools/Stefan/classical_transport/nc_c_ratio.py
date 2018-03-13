#!/usr/bin/env python
import sys
from neoclassical_flux import neoclassical_flux
from classical_flux import classical_flux

def nc_c_ratio(inputFilename,max_m,maxabs_n,symmetry,signcorr,BBar=1,RBar=1):
    nc = neoclassical_flux(inputFilename)
    c = classical_flux(inputFilename,max_m,maxabs_n,symmetry,signcorr,BBar,RBar)
    nc_over_c = nc/c
    return nc_over_c

if __name__=="__main__":
    argv = sys.argv
    if len(argv) == 2:
        namelist_filename = sys.argv[1]
        signcorr = 1
    elif len(argv) == 3:
        namelist_filename = sys.argv[1]
        signcorr = int(sys.argv[2])
    else:
        raise Warning("Too many input arguments to " + argv[0])
        namelist_filename = sys.argv[1]
        signcorr = int(sys.argv[2])
    
    max_m = float("inf")
    maxabs_n = float("inf")
    symmetry = "StelSym"
    nc_over_c = nc_c_ratio(namelist_filename,max_m,maxabs_n,symmetry,signcorr)
    print str(nc_over_c)[1:-1]
