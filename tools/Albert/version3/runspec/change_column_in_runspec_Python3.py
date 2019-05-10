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

externalRadiusFile = 'rho.txt'
#externalDataFile = 'nArHe_180919055.txt'
externalDataFile = 'nArHe_180919046.txt'

#runspecFilename = 'shot180919055_profAL_TileTerunspec21.dat'
runspecFilename = 'shot180919046_profAL_TileTerunspec21.dat'

outputSuffix = "_real_nz"

radiusColumn = 0

columnsToChange = [9, 15] #Density and density gradient
TransformFactor = 10**(-6)

minorRadius = 0.51092

MinFloat = pow(10, -sys.float_info.dig) 

##############################
##########END INPUTS##########
##############################

inputRadius = np.genfromtxt(externalRadiusFile, dtype=None, comments="!", skip_header=0)
inputData = np.genfromtxt(externalDataFile, dtype=None, comments="!", skip_header=0)

print(inputRadius)
print("Input radii:")
print(inputRadius * minorRadius)
print("Input data:")
print(inputData)

interpFunc = interpolate.interp1d(inputRadius * minorRadius, inputData, fill_value="extrapolate")

inputRunspec = np.genfromtxt(runspecFilename, dtype='float', comments="!", skip_header=0)

with open(runspecFilename) as f:
    dataHeader = f.readline().strip()
print("Header:")
print(dataHeader)

print("Data in file to change:")
print(inputRunspec.T)

print("Columns to change:")
print(columnsToChange[0])
print(columnsToChange[1])

print("Radii:")
print(inputRunspec[:, radiusColumn])

#for columnTC in columnsToChange:
print("Old density column")
print(inputRunspec[:, columnsToChange[0]])
inputRunspec[:, columnsToChange[0]] = TransformFactor * interpFunc(inputRunspec[:, radiusColumn])
print("New density column")
print(inputRunspec[:, columnsToChange[0]])

print("Old density derivative column")
print(inputRunspec[:, columnsToChange[1]])
inputRunspec[:, columnsToChange[1]] = TransformFactor * np.gradient(interpFunc(inputRunspec[:, radiusColumn]), inputRunspec[:, radiusColumn])
print("New density derivative column")
print(inputRunspec[:, columnsToChange[1]])

print("New data in file to change:")
print(inputRunspec.T)


outputFilename = (runspecFilename.split('.', maxsplit=1))[0] + outputSuffix + '.' + (runspecFilename.split('.', maxsplit=1))[1]

print("Saving data to " + outputFilename)
np.savetxt(outputFilename, inputRunspec, delimiter=' ', newline='\n', header=dataHeader, comments='')

sys.exit(0)
