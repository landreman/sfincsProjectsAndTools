#!/usr/bin/env python 

import numpy 
import numpy as np 
import os, sys, inspect, math
import warnings  
import subprocess
from scipy import interpolate

print ("This is "+ inspect.getfile(inspect.currentframe()))

sfincsHome = os.environ.get('SFINCS_HOME') 
sfincsProjectsAndToolsHome = os.environ.get('SFINCS_PROJECTS_AND_TOOLS_HOME')

##INPUTS##

outputFilename = 'W7X_180919.055_XICS.dat'
inputFiles = ['rho.txt', 'normAr16Flux.txt', 'normAr16FluxErrors.txt', 'D.txt', 'Derror.txt', 'v.txt', 'verror.txt']
outputLabels = ['r/a', 'GammaAr16+/nAr16+', 'ErrorGammaAr16+/nAr16+', 'DAr16+', 'ErrorDAr16+', 'VAr16+', 'ErrorVAr16+']
TransformFactors = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

MinFloat = pow(10, -sys.float_info.dig) 

##############################
##########END INPUTS##########
##############################

OutputData = np.array([]).transpose()
counter = -1

for inputFile in inputFiles:
    counter += 1
    
    inputData = np.genfromtxt(inputFile, dtype=None, comments="#", skip_header=0)
    print(str(inputData))
    if OutputData.size == 0:
        OutputData = np.array([inputData]).transpose()
    else :
        OutputData = np.concatenate((OutputData, TransformFactors[counter]*np.array([inputData]).transpose()), axis=1)

print("Saving data to " + outputFilename)
np.savetxt(outputFilename, OutputData, delimiter='\t', newline='\n', header='\t'.join(outputLabels))

sys.exit(0)

