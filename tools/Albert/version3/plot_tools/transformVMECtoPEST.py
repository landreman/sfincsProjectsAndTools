#!/usr/bin/env python

import matplotlib
import matplotlib.pyplot as plt
import h5py
import numpy as np
import scipy.linalg
import os, sys, inspect
import warnings
import matplotlib.ticker as ticker
from numpy import pi
from numpy import nan
import math
from netCDF4 import Dataset
from shutil import copyfile, copy, copy2
from VMECtoPEST_functions import interp2_cyclic, griddatacyclic

#############
def transformVMECtoPEST(woutin, s_wish, Vmec_theta, Vmec_zeta, handedness):
##INPUT PARAMETERS##
# woutin (string) is the file name of the netcdf file of the equilibrium
# s_wish (float in [0, 1]) is the wanted flux surface in normalized flux
# Vmec_theta (array of floats) is the grid points in theta
# Vmec_zeta (array of floats) is the grid points in zeta
# handedness = 1  => right-handed output structs
# handedness =-1  => left-handed output structs
####################

   copyfile = '.'.join((woutin.split('.'))[:-1]) + '_COPY.'+ (woutin.split('.'))[-1]

   ##OPEN NETCDF FILE##
   copy2(woutin, copyfile)

   wout = Dataset(copyfile, "r+", format="NETCDF4")

   signchange = handedness*np.float_(wout["signgs"][...]) # is handedness*(-1), because vmec is left handed

   Geom_Nperiods = np.float_(wout["nfp"][...])    #!< number of field periods

   if len(wout.get_variables_by_attributes(name='iasym')) < 1:
      wout.createVariable('iasym', 'i1')
      wout['iasym'][...] = (len(wout.get_variables_by_attributes(name='bmns')) > 0)

   Geom_StelSym  = not(wout['iasym'][...]) #!<  defines stellarator symmetry for iasym=0, otherwise =1

   Geom_fullgrid_s     = wout['phi'][...]/(wout['phi'][...])[wout['ns'][...] - 1]    #full grid

   skip = 1 #this is how many elements are skipped at low radii when going to half grid

   Geom_s = (np.array(Geom_fullgrid_s[skip-1:-1]) + np.array(Geom_fullgrid_s[skip:]))/2 #half grid

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% specific for this radius
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   rindh = np.argmin(np.absolute(Geom_s-s_wish))
   skrindh = skip + rindh #shorthand notation.

#now use a RH coordinate system (u,w,s) with w=-v
   if not(Geom_StelSym):
      print "Non-stelllarator symmmetric case not implemented!"
      raise

#   Vmec_handedness = handedness

   Vmec_vmecu, Vmec_vmecw = np.meshgrid(Vmec_theta, Vmec_zeta) ##Could be meshgrid() instead of mgrid()
   Vmec_vmecu = Vmec_vmecu.transpose()
   Vmec_vmecw = Vmec_vmecw.transpose()
   
   Vmec_llist_m = np.float_(wout['xm'][...])
   Vmec_llist_n = signchange*np.float_(wout['xn'][...])/Geom_Nperiods
   Vmec_llist_data = (wout['lmns'][...])[skrindh,:]

   assert len(Vmec_llist_data) == len(Vmec_llist_m) and len(Vmec_llist_data) == len(Vmec_llist_n)

   Vmec_l = np.zeros(np.shape(Vmec_vmecu))

   for i in range(0, len(Vmec_llist_data)):
      Vmec_l += Vmec_llist_data[i] * np.sin(Vmec_llist_m[i]*Vmec_vmecu - Vmec_llist_n[i]*Geom_Nperiods*Vmec_vmecw) 



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Make the Pest coordinates
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Vmec_Dpthetavmecu = Vmec_l
   Vmec_ptheta = Vmec_Dpthetavmecu + Vmec_vmecu #This is the "transformation"
   Vmec_pzeta = Vmec_vmecw

   Pest_pzeta = Vmec_vmecw
   Pest_ptheta = Vmec_vmecu

   Pest_Dpthetavmecu = griddatacyclic(Vmec_ptheta, Vmec_pzeta, Vmec_Dpthetavmecu, Geom_Nperiods)[0]
   Pest_vmecu = Pest_ptheta - Pest_Dpthetavmecu #This is the "transformation"
   Pest_vmecw = Pest_pzeta
   
   wout.close()
   os.remove(copyfile)

   return (Pest_vmecu, Pest_vmecw, Vmec_vmecu, Vmec_vmecw, Geom_Nperiods)



