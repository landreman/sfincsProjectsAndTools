import numpy as np
import sys,os
import h5py


def read_neoclassical_flux(radialCoordinate,sfincsOutputFilename):   
    if radialCoordinate == 0:
        # wish coordinate is psiHat
        particlefluxname = "particleFlux_vd_psiHat"
        heatfluxname = "heatFlux_vd_psiHat"
    elif radialCoordinate == 1:
        # coordinate is psiN
        particlefluxname = "particleFlux_vd_psiN"
        heatfluxname = "heatFlux_vd_psiN"
    elif radialCoordinate == 2:
        # coordinate is rHat
        particlefluxname = "particleFlux_vd_rHat"
        heatfluxname = "heatFlux_vd_rHat"
    elif radialCoordinate == 3:
        # coordinate is rN
        particlefluxname = "particleFlux_vd_rN"
        heatfluxname = "heatFlux_vd_rN"

    try:
        sfincsOutput = h5py.File(sfincsOutputFilename,'r')
    except IOError:
        raise IOError("SFINCS output file '" + sfincsOutputFilename +"' not found.")

    try:
        particleFlux = sfincsOutput["/"+particlefluxname][:,-1] #get last iteration for all species
        heatFlux = sfincsOutput["/"+heatfluxname][:,-1] #get last iteration for all species
    except:
        pass
    return particleFlux, heatFlux


if __name__=="__main__":
    sfincsOutputFilename = sys.argv[1]
    inputRadialCoordinate = 3
    particleFlux = read_neoclassical_flux(inputRadialCoordinate,sfincsOutputFilename)
    print particleFlux
