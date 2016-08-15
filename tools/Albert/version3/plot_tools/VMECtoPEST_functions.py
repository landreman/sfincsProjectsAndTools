#!/usr/bin/env python

#This python file contains functions for transforming from VMEC to PEST coordinates

import matplotlib
import matplotlib.pyplot as plt
import h5py
import numpy as np
import numpy
import scipy.linalg
import os, sys, inspect
import warnings
import matplotlib.ticker as ticker
import math
from numpy import pi
from numpy import nan
from scipy import interpolate


##THIS FILE IS NOT IMPLEMENTED YET##

def interp2_cyclic(u, v, F, uout, vout, N) :

#% u must be 2*pi periodic.
#% v must be 2*pi/N periodic.

    vDirection = np.sign(v[0,-1] - v[0,0])
    uDirection = np.sign(u[-1,0] - u[0,0])

    #print v
    #print v.shape
    #print v[:, 0] + vDirection*2*pi/N
    #print ((v[:, 0] + vDirection*2*pi/N).shape)[0]
    #print (v[:, 0].reshape((v[:, 0].shape)[0], 1) + vDirection*2*pi/N)
    #print (v[:, 0].reshape((v[:, 0].shape)[0], 1) + vDirection*2*pi/N).shape

    vbig = np.concatenate((v, (v[:, 0].reshape((v[:, 0].shape)[0], 1) + vDirection*2*pi/N)), axis=1)    
    vbig = np.concatenate((vbig, vbig[0,:].reshape((vbig[0,:].shape)[0], 1).T))

    ubig = np.concatenate((u, u[0,:].reshape((u[0,:].shape)[0], 1).T + uDirection*2*pi))
    ubig = np.concatenate((ubig, ubig[:,0].reshape((ubig[:,0].shape)[0], 1)), axis=1)

    Fbig = np.concatenate((F, F[:,0].reshape((F[:,0].shape)[0], 1)), axis=1)
    Fbig = np.concatenate((Fbig, Fbig[0,:].reshape((Fbig[0,:].shape)[0], 1).T))

#    return interp2(vbig, ubig, Fbig, np.mod(vout,2*pi/N), np.mod(uout,2*pi))
    #print "vbig"
    #print vbig
    #print "ubig"
    #print ubig
    #print "Fbig"
    #print Fbig

    fInt = interpolate.interp2d(vbig, ubig, Fbig)

    #print np.mod(vout,2*pi/N).shape
    #print np.mod(uout,2*pi).shape
    #print fInt(np.mod(vout,2*pi/N).flatten(), np.mod(uout,2*pi).flatten()).shape

    #print "test"

    fToReturn = nan*numpy.ones(np.mod(uout,2*pi).shape)
    for ii in range(0, (fToReturn.shape)[0]) :
        for jj in range(0, (fToReturn.shape)[1]) :
            fToReturn[ii,jj] = fInt(np.mod(vout,2*pi/N)[ii,jj], np.mod(uout,2*pi)[ii,jj])

    #print fToReturn
    #print fToReturn.shape

    #return fInt(np.mod(vout,2*pi/N), np.mod(uout,2*pi))
    return fToReturn





#####################################################################

def griddatacyclic(unoneq, vnoneq, Fnoneq, N):
#% u must be 2*pi periodic.
#% v must be 2*pi/N periodic.
#% (unoneq, vnoneq) is the non-equidistant grid on which F is given.
#% (u,v) is the equidistant grid generated with ndgrid which to interpolate to.


    Nu = Fnoneq.shape[0]
    Nv = Fnoneq.shape[1]

    Dv = 2*pi/Nv/N
    Du = 2*pi/Nu
    vvec = numpy.arange(0, Nv)*Dv
    uvec = (numpy.arange(0, Nu)).conj().transpose()*Du
    u, v = numpy.meshgrid(uvec,vvec)
    u = u.transpose() ##MAYBE REMOVE THIS
    v = v.transpose() ##MAYBE REMOVE THIS

    Mleftadd  = numpy.maximum(0, numpy.ceil(numpy.amax(vnoneq[:,0])/(2*pi/N)))
    Mrightadd = numpy.maximum(0, numpy.ceil(1 - numpy.amin(vnoneq[:,-1])/(2*pi/N)))
    Mv = Mleftadd + 1 + Mrightadd
    Mv = Mv.astype(int)

    Mbottomadd  = numpy.maximum(0, numpy.ceil(numpy.amax(unoneq[0,:])/(2*pi)))
    Mtopadd     = numpy.maximum(0, numpy.ceil(1 - numpy.amin(unoneq[0,:])/(2*pi)))
    Mu = Mtopadd + 1 + Mbottomadd
    Mu = Mu.astype(int)
    
    unoneqBig = nan*numpy.zeros((Nu*Mu,Nv*Mv))
    vnoneqBig = nan*numpy.zeros((Nu*Mu,Nv*Mv))
    FnoneqBig = nan*numpy.zeros((Nu*Mu,Nv*Mv))
    
    for mui in range(1, Mu + 1):
        for mvi in range(1, Mv + 1):
            vadd = (mui - (Mleftadd + 1))*2*pi/N
            uadd = (mvi - (Mbottomadd + 1))*2*pi
            unoneqBig[(mui-1)*Nu:mui*Nu][:, (mvi-1)*Nv:mvi*Nv] = unoneq + uadd
            vnoneqBig[(mui-1)*Nu:mui*Nu][:, (mvi-1)*Nv:mvi*Nv] = vnoneq + vadd
            FnoneqBig[(mui-1)*Nu:mui*Nu][:, (mvi-1)*Nv:mvi*Nv] = Fnoneq

    #F = griddata(unoneqBig,vnoneqBig,FnoneqBig,u,v)
    #unoneqBig = unoneqBig.flatten()
    #vnoneqBig = vnoneqBig.flatten()
    #FnoneqBig = FnoneqBig.flatten()
    #test = interpolate.interp2d(unoneqBig,vnoneqBig,FnoneqBig)
    #print "u"
    #print u
    #print u.shape
    #print u.flatten()
    #print u.flatten().shape
    #print uvec
    #print "v"
    #print v
    #print v.shape
    #print v.flatten()
    #print v.flatten().shape
    #print vvec


    
    #F = (interpolate.interp2d(unoneqBig,vnoneqBig,FnoneqBig))(u,v)
    F = (interpolate.interp2d(unoneqBig,vnoneqBig,FnoneqBig))(uvec,vvec)
    F = F.transpose() ##MAYBE REMOVE THIS
    #print "F"
    #print F
    #print F.shape

    return (F, u, v, Nu, Nv, Du, Dv)
