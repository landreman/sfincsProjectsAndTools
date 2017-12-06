#!/usr/bin/env python
# -*- coding: utf-8 -*-

outputFilename = "sfincsOutput.h5"

import matplotlib
#import matplotlib.pyplot as plt
import h5py
import numpy
import numpy as np
import inspect, math, os
import pickle
import sys
from subprocess import call

print ("This is "+ inspect.getfile(inspect.currentframe()))

import matplotlib.pyplot as plt

#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True

numRuns = 0
scanVariableValues = []
outputs = []

##INPUTS##
e = 1.6021766208*10**(-19) #Electron charge
mbar = 1.672621777*10**(-27) #Proton mass
nbar = 10.0**20 # in 10^20
Tbar = 1000*e 
Rbar = 1.0
Bbar = 1.0
vbar = np.sqrt(2.0 * Tbar / mbar)

filename = 'sfincsOutput.h5'


#xVariableName = "nHats_1"


MinFloat = pow(10, -sys.float_info.dig)
try:
    f = h5py.File(filename,'r') 
except:
    print ("Cannot open " + filename)
    print ("Must be an HDF5 file in current folder")
    sys.exit()
Zs = f["Zs"][()]
mHats = f["mHats"][()]
THats = f["THats"][()]
nHats = f["nHats"][()]
B0OverBBar = f["B0OverBBar"][()]
GHat = f["GHat"][()]
IHat = f["IHat"][()]
iota = f["iota"][()]
nu_n = f["nu_n"][()]
FSABFlow = f["FSABFlow"][()]
FSABjHat = f["FSABjHat"][()]
includePhi1 = f["includePhi1"][()]
integerToRepresentTrue = f["integerToRepresentTrue"][()]
if includePhi1 == integerToRepresentTrue:
    particleFlux_vd_rHat = f["particleFlux_vd_rHat"][()]
else:
    particleFlux_vd_rHat = f["particleFlux_vm_rHat"][()]
f.close()  

##nu_n = nuBar*Rbar/vbar
##nuBar = 4*sqrt(2*pi)*nbar*e^4*ln(Lambda) / (3*sqrt(mbar)*Tbar^(3/2)) in Gaussian units

DensityToNuFactor = np.absolute((GHat + iota*IHat)*nu_n*Zs**4 / (B0OverBBar*THats**2))

print ("*********************************") 
print ("vbar: " + str(vbar)) 
print ("*********************************") 
print ("Zs: " + str(Zs))
print ("mHats: " + str(mHats))
print ("THats: " + str(THats))
print ("nHats: " + str(nHats))
print ("B0OverBBar: " + str(B0OverBBar))
print ("GHat: " + str(GHat))
print ("IHat: " + str(IHat))
print ("iota: " + str(iota))
print ("nu_n: " + str(nu_n))
print ("NuPrime_aa: " + str(DensityToNuFactor*nHats))
print ("NuPrime_a = SUM_b (NuPrime_ab): " + str(DensityToNuFactor*nHats * sum(Zs**2 * nHats) / (Zs**2 * nHats))) 
print ("*********************************")
print ("particleFlux_vd_rHat: " + str(particleFlux_vd_rHat[:,-1]))
print ("particleFlux_vd_rHat * vbar / nHat: " + str(particleFlux_vd_rHat[:,-1] * vbar / nHats))
