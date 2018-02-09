#!/usr/bin/env python
import sys,os
from read_namelist import readVariable
from read_neoclassical_flux import read_neoclassical_flux

def neoclassical_flux(inputFilename):

    inputRadialCoordinate = readVariable("inputRadialCoordinate",inputFilename,"int",False)
    if inputRadialCoordinate is None:
        inputRadialCoordinate = 3
    outputFilename = readVariable("outputFilename",inputFilename,"string",False)
    if outputFilename is None:
        outputFilename = "sfincsOutput.h5"
    else:
        if outputFilename[0] == '"' and outputFilename[-1] == '"':
            outputFilename = outputFilename[1:-1]

    outputFilename = os.path.dirname(os.path.abspath(inputFilename)) +"/" + outputFilename

    particleFlux = read_neoclassical_flux(0,outputFilename)
    return particleFlux

if __name__=="__main__":
    inputFilename = sys.argv[1]
    particleFlux = neoclassical_flux(inputFilename)
    print str(particleFlux)[1:-1]
