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


def interp2_cyclic(u, v, F, uout, vout, N) :

#% Returns F(uout, vout) by interpolating F(u, v).
#% (u, v) is the equidistant input grid and (uout, vout) 
#% the output grid (which can be non-equidistant).
#% u and uout must be 2*pi periodic.
#% v and vout must be 2*pi/N periodic.
#% The size of u, v and F must be equal.
#% The size of uout and vout must be equal.

    ##Check that input is OK
    assert u.shape == v.shape and u.shape == F.shape
    assert uout.shape == vout.shape
    assert N > 0

    vDirection = np.sign(v[0,-1] - v[0,0])
    uDirection = np.sign(u[-1,0] - u[0,0])

    vbig = np.concatenate((v, (v[:, 0].reshape((v[:, 0].shape)[0], 1) + vDirection*2*pi/N)), axis=1)    
    vbig = np.concatenate((vbig, vbig[0,:].reshape((vbig[0,:].shape)[0], 1).T))

    ubig = np.concatenate((u, u[0,:].reshape((u[0,:].shape)[0], 1).T + uDirection*2*pi))
    ubig = np.concatenate((ubig, ubig[:,0].reshape((ubig[:,0].shape)[0], 1)), axis=1)


    Fbig = np.concatenate((F, F[:,0].reshape((F[:,0].shape)[0], 1)), axis=1)
    Fbig = np.concatenate((Fbig, Fbig[0,:].reshape((Fbig[0,:].shape)[0], 1).T))

    fInt = interpolate.interp2d(ubig[:,0], vbig[0,:], Fbig.transpose())

    ##TO CHECK THAT ALL ARRAY VALUES ARE FILLED WE START WITH NaN:s
    fToReturn = nan*numpy.ones(np.mod(uout,2*pi).shape)

    for ii in range(0, (fToReturn.shape)[0]) :
        for jj in range(0, (fToReturn.shape)[1]) :
            fToReturn[ii,jj] = fInt(np.mod(uout,2*pi)[ii,jj], np.mod(vout,2*pi/N)[ii,jj])

    return fToReturn



#####################################################################

def griddatacyclic(unoneq, vnoneq, Fnoneq, N):

#% Takes a function Fnoneq(unoneq, vnoneq) on a non-equidistant grid and returns
#% it interpolated to an equidistant grid of the same dimension.
#% unoneq must be 2*pi periodic. u is also 2*pi periodic.
#% vnoneq must be 2*pi/N periodic. v is also 2*pi/N periodic.
#% (unoneq, vnoneq) is the non-equidistant grid on which F is given.
#% (u, v) is the equidistant grid generated with meshgrid which to interpolate to.


    assert unoneq.shape == vnoneq.shape and unoneq.shape == Fnoneq.shape
    assert N > 0

    Nu = Fnoneq.shape[0]
    Nv = Fnoneq.shape[1]

    Dv = 2*pi/Nv/N
    Du = 2*pi/Nu

    Mleftadd  = numpy.maximum(0, numpy.ceil(numpy.amax(vnoneq[0,:])/(2*pi/N)))
    Mrightadd = numpy.maximum(0, numpy.ceil(1 - numpy.amin(vnoneq[-1,:])/(2*pi/N)))
    Mv = Mleftadd + 1 + Mrightadd
    Mv = Mv.astype(int)

    Mbottomadd  = numpy.maximum(0, numpy.ceil(numpy.amax(unoneq[:,0])/(2*pi)))
    Mtopadd     = numpy.maximum(0, numpy.ceil(1 - numpy.amin(unoneq[:,0])/(2*pi)))
    Mu = Mtopadd + 1 + Mbottomadd
    Mu = Mu.astype(int)

    #TO CHECK THAT ALL ARRAY VALUES ARE FILLED WE START WITH NaN:s
    unoneqBig = nan*numpy.ones((Nu*Mu,Nv*Mv))
    vnoneqBig = nan*numpy.ones((Nu*Mu,Nv*Mv))
    FnoneqBig = nan*numpy.ones((Nu*Mu,Nv*Mv))
    
    for mui in range(1, Mu + 1):
        for mvi in range(1, Mv + 1):
            vadd = (mui - (Mleftadd + 1))*2*pi/N
            uadd = (mvi - (Mbottomadd + 1))*2*pi
            unoneqBig[(mui-1)*Nu:mui*Nu][:, (mvi-1)*Nv:mvi*Nv] = unoneq + uadd
            vnoneqBig[(mui-1)*Nu:mui*Nu][:, (mvi-1)*Nv:mvi*Nv] = vnoneq + vadd
            FnoneqBig[(mui-1)*Nu:mui*Nu][:, (mvi-1)*Nv:mvi*Nv] = Fnoneq

    vvec = numpy.arange(0, Nv)*Dv
    uvec = (numpy.arange(0, Nu))*Du
    u, v = numpy.meshgrid(uvec,vvec)

    #NOTE! THE scipy FUNCTION interp2d SEEMS TO BE PERFORMING BAD FOR SOME CASES
    #USE griddata INSTEAD

    unoneqBig_flatten = unoneqBig.flatten()
    vnoneqBig_flatten = vnoneqBig.flatten()
    FnoneqBig_flatten = FnoneqBig.flatten()

    assert len(unoneqBig_flatten) == len(vnoneqBig_flatten) and len(unoneqBig_flatten) == len(FnoneqBig_flatten)


    #F = (interpolate.interp2d(unoneqBig_flatten,vnoneqBig_flatten,FnoneqBig_flatten))(uvec,vvec)
    F = interpolate.griddata(np.array((unoneqBig_flatten, vnoneqBig_flatten)).T, FnoneqBig_flatten, (u, v))

    F = F.transpose() ##MAYBE REMOVE THIS
    u = u.transpose() ##MAYBE REMOVE THIS
    v = v.transpose() ##MAYBE REMOVE THIS


    return (F, u, v)
