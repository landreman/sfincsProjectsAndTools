import numpy as np
import sys,os
import h5py


def read_neoclassical_flux(radialCoordinate,sfincsOutputFilename):   
    if radialCoordinate == 0:
        # wish coordinate is psiHat
        fluxname = "particleFlux_vd_psiHat"
    elif radialCoordinate == 1:
        # coordinate is psiN
        fluxname = "particleFlux_vd_psiN"
    elif radialCoordinate == 2:
        # coordinate is rHat
        fluxname = "particleFlux_vd_rHat"
    elif radialCoordinate == 3:
        # coordinate is rN
        fluxname = "particleFlux_vd_rN"

    try:
        sfincsOutput = h5py.File(sfincsOutputFilename,'r')
    except IOError:
        raise IOError("SFINCS output file '" + sfincsOutputFilename +"' not found.")

    try:
        particleFlux = sfincsOutput["/"+fluxname][:,-1] #get last iteration for all species
    except:
        pass
    return particleFlux


if __name__=="__main__":
    sfincsOutputFilename = sys.argv[1]
    inputRadialCoordinate = 3
    particleFlux = read_neoclassical_flux(inputRadialCoordinate,sfincsOutputFilename)
    print particleFlux
