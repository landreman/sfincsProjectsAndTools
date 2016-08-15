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
#import subprocess
from VMECtoPEST_functions import interp2_cyclic, griddatacyclic

#############



def transformVMECtoPEST(woutin,s_wish,Nu,Nw,handedness):

# woutin is the wout file name or just the netcdf variables from the wout file
# (Too old matlab versions do not have the necessary netcdf routines.)
#
# s_wish is the wanted flux surface
#
# (Nu, Nw) is the (poloidal, toroidal) resolution
#
# handedness = 1  => right-handed output structs
# handedness =-1  => left-handed output structs

   ##SECTION COMMENTED BY AM##
#!   if not(isstruct(woutin)):#if not already loaded, assume woutin is a string with the file name
#!      wout=struct()
#!      if isstr(woutin):
#!         tmp=ncinfo(woutin)
#!         for vi=1:length(tmp.Variables):
#!            wout=setfield(wout,tmp.Variables(vi).Name, ncread(woutin,tmp.Variables(vi).Name))
#!      else:
#!         error('Unknown input! woutin must be a file name or the netcdf variables from the wout file!')
   ###########################
         
   #    if handedness=1 (right handed) let w = -v!
   # else if handedness=-1 (left handed) let w = v!

   inputfile = "/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/COPY_wout_lhd2.nc"
   copyfile = "/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/wout_lhd2_COPY.nc"

   ##OPEN NETCDF FILE BY AM##
   #copy2("/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/COPY_wout_lhd2.nc", "/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/wout_lhd2_COPY.nc")
   #copy2("/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/Input.nb", "/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/Input_COPY.nb")
   #subprocess.call(["cp", "/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/COPY_wout_lhd2.nc", "/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/wout_lhd2_COPY.nc"], shell=True)
   #os.system ("cp %s %s" % ("/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/COPY_wout_lhd2.nc", "/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/wout_lhd2_COPY.nc"))
   #subprocess.call(["cp", "/afs/ipp/home/a/almo/Phi1/LHD/lhd2_B_III/Input/COPY_wout_lhd2.nc", "/afs/ipp/home/a/almo/Phi1/LHD/lhd2_B_III/Input/wout_lhd2_COPY.nc"])

   copy2(inputfile, copyfile)

   wout = Dataset("/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/wout_lhd2_COPY.nc", "r+", format="NETCDF4")
   print wout["ntor"][...]
   print wout["ntor"].getValue()
   print wout["ntor"]
   #print wout["qwerty"]
   #print wout["qwerty"].getValue()
   #print wout.variables
   print wout.get_variables_by_attributes(name='northward_sea_water_velocity')
   print wout.get_variables_by_attributes(name='ntor')
   print len(wout.get_variables_by_attributes(name='northward_sea_water_velocity'))
   print len(wout.get_variables_by_attributes(name='ntor'))
   print len(wout.get_variables_by_attributes(name='northward_sea_water_velocity')) < 1
   print len(wout.get_variables_by_attributes(name='ntor')) < 1
   #wout.close()
   #sys.exit()
   ##########################
   #wout[""].getValue()


   signchange = handedness*np.float_(wout["signgs"][...]) # is handedness*(-1), because vmec is left handed

   Geom_headertext_input_extension = wout["input_extension"][...]#!<  suffix of the vmec-input file: input.<input_extension>

   Geom_m0b      = np.float_(wout["mpol"][...])    #!< number of poloidal fourier modes
   Geom_n0b      = np.float_(wout["ntor"][...])    #!< upper bound toroidal fourier modes: -ntor <= n <= ntor
   Geom_nsurf    = np.float_(wout["ns"][...])    #!< number of radial surfaces
   Geom_Nperiods = np.float_(wout["nfp"][...])    #!< number of field periods
   
   #if not(isfield(wout,'iasym')):
   #   wout.iasym = isfield(wout,'bmns')

   if len(wout.get_variables_by_attributes(name='iasym')) < 1:
      print "TEST"
      wout.createVariable('iasym', 'i1')
      wout['iasym'][...] = (len(wout.get_variables_by_attributes(name='bmns')) > 0)

   print wout['iasym'][...]

   print Geom_Nperiods



   Geom_handedness = handedness
   Geom_StelSym  = not(wout['iasym'][...]) #!<  defines stellarator symmetry for iasym=0, otherwise =
   #Geom_torfluxtot = wout.phi(wout.ns)*signchange  #!<  total toroidal flux within the boundary (s=1)
   Geom_torfluxtot = (wout['phi'][...])[wout['ns'][...] - 1]*signchange
   Geom_psi_a = Geom_torfluxtot/2/pi

   Geom_minorradiusVMEC = wout['Aminor_p'][...]  #!<  minor plasma radius
   Geom_majorradiusLastbcR00 = nan #Calculate this below (not necessary)
   Geom_minorradiusW7AS = nan #Calculate this below (not necessary)
   Geom_majorradiusVMEC =  wout['Rmajor_p'][...]  #!<  major plasma radius

   Geom_fullgrid_s     = wout['phi'][...]/(wout['phi'][...])[wout['ns'][...] - 1]    #full grid
   Geom_fullgrid_rnorm = np.sqrt(Geom_fullgrid_s)     #full grid

   skip = 1 #this is how many elements are skipped at low radii when going to half grid

   Geom_s = (np.array(Geom_fullgrid_s[skip-1:-2]) + np.array(Geom_fullgrid_s[skip:-1]))/2 #half grid
   Geom_rnorm = np.sqrt(Geom_s) #half grid

#Geom.dVdsoverNper=dVdsoverNper*signchange;
   Geom_Bphi  = wout['bvco'][...]*signchange#direction switch sign
   Geom_Btheta = wout['buco'][...]#*signchange

   Geom_Bphi = Geom_Bphi[skip:-1]
   Geom_Btheta = Geom_Btheta[skip:-1]

   Geom_fullgrid_iota = wout['iotaf'][...]*signchange
   Geom_iota = wout['iotas'][...]*signchange #half mesh
   Geom_iota = Geom_iota[skip:-1]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% specific for this radius
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   #[dummy,rindh]=min(np.absolute(Geom_s-s_wish))
   rindh = np.argmin(np.absolute(Geom_s-s_wish))
   #print rindh
   dummy = np.amin(np.absolute(Geom_s-s_wish))
   #print dummy
   s = Geom_s[rindh]

   Pest_handedness = handedness
   Pest_s = Geom_s[rindh]
   Pest_rnorm = Geom_rnorm[rindh]
   Pest_iota = Geom_iota[rindh]
   Pest_Bphi = Geom_Bphi[rindh]
   Pest_Btheta = Geom_Btheta[rindh]

   G = Geom_Bphi[rindh]   #=wout.bsubvmnc(n0ind,m0ind,rind)
   I = Geom_Btheta[rindh] #=wout.bsubumnc(n0ind,m0ind,rind)
   iota = Geom_iota[rindh]

   skrindh = skip + rindh #shorthand notation.
   rindf_plus = skrindh
   rindf_minus = skrindh - 1



#now use a RH coordinate system (u,w,s) with w=-v
   if not(Geom_StelSym):
      #error('Non-stelllarator symmmetric case not implemented!')
      print "Non-stelllarator symmmetric case not implemented!"
      raise

   Vmec_handedness = handedness
   Vmec_Nvmecu = Nu
   Vmec_Nvmecw = Nw

   ##Make VMEC grid (moved from below)
   Vmec_Dvmecw = 2*pi/Vmec_Nvmecw/Geom_Nperiods
   Vmec_Dvmecu = 2*pi/Vmec_Nvmecu
   Vmec_vmecuvec = np.arange(0,Vmec_Nvmecu)*Vmec_Dvmecu
   Vmec_vmecwvec = (np.arange(0,Vmec_Nvmecw)).conj().transpose()*Vmec_Dvmecw
   Vmec_vmecu, Vmec_vmecw = np.meshgrid(Vmec_vmecuvec, Vmec_vmecwvec) ##Could be meshgrid() instead of mgrid()
   ##

   print Vmec_Dvmecw
   print Vmec_Dvmecu
   print Vmec_vmecuvec
   print len(Vmec_vmecuvec)
   print ""
   print Vmec_vmecwvec
   print len(Vmec_vmecwvec)
   print ""
   print Vmec_vmecu
   print Vmec_vmecu.shape
   print ""
   print Vmec_vmecw
   print Vmec_vmecw.shape
   print ""
   
   Vmec_Rmnlist_m = np.float_(wout['xm'][...])
   Vmec_Rmnlist_n = signchange*np.float_(wout['xn'][...])/Geom_Nperiods
   Vmec_Rmnlist_cosparity = np.ones(np.shape(Vmec_Rmnlist_m))
   Vmec_Rmnlist_data = ((wout['rmnc'][...])[rindf_minus,:] + (wout['rmnc'][...])[rindf_plus,:])/2

   print Vmec_Rmnlist_m
   print len(Vmec_Rmnlist_m)
   print ""
   print Vmec_Rmnlist_n
   print len(Vmec_Rmnlist_n)
   print ""
   print Vmec_Rmnlist_cosparity
   print len(Vmec_Rmnlist_cosparity)
   print ""
   print Vmec_Rmnlist_data
   print len(Vmec_Rmnlist_data)
   print ""

#   Vmec_Rmn = mnmat(Vmec_Rmnlist,Vmec_Nvmecu,Vmec_Nvmecw,'forceSize')               
#   Vmec_R = ifftmn(Vmec_Rmn,Geom_Nperiods,Vmec_Nvmecu,Vmec_Nvmecw,'forceSize')       

   Vmec_R = np.zeros(np.shape(Vmec_vmecu))
   print Vmec_R
   print Vmec_R.shape
   print ""

   for i in range(0, len(Vmec_Rmnlist_data)):
      Vmec_R += Vmec_Rmnlist_data[i] * np.cos(Vmec_Rmnlist_m[i]*Vmec_vmecu - Vmec_Rmnlist_n[i]*Geom_Nperiods*Vmec_vmecw)

   print Vmec_R
   print Vmec_R.shape
   print ""
   
   Vmec_Zmnlist_m = np.float_(wout['xm'][...])
   Vmec_Zmnlist_n = signchange*np.float_(wout['xn'][...])/Geom_Nperiods
   Vmec_Zmnlist_cosparity = 0*np.ones(np.shape(Vmec_Zmnlist_m))
   Vmec_Zmnlist_data = ((wout['zmns'][...])[rindf_minus,:] + (wout['zmns'][...])[rindf_plus,:])/2
#   Vmec_Zmn = mnmat(Vmec_Zmnlist, Vmec_Nvmecu, Vmec_Nvmecw, 'forceSize')
#   Vmec_Z = ifftmn(Vmec_Zmn, Geom_Nperiods, Vmec_Nvmecu, Vmec_Nvmecw,'forceSize')
   Vmec_Z = np.zeros(np.shape(Vmec_vmecu))

   for i in range(0, len(Vmec_Zmnlist_data)):
      Vmec_Z += Vmec_Zmnlist_data[i] * np.sin(Vmec_Zmnlist_m[i]*Vmec_vmecu - Vmec_Zmnlist_n[i]*Geom_Nperiods*Vmec_vmecw)
   
   Vmec_Bmnlist_m = np.float_(wout['xm_nyq'][...])
   Vmec_Bmnlist_n = signchange*np.float_(wout['xn_nyq'][...])/Geom_Nperiods
   Vmec_Bmnlist_cosparity = np.ones(np.shape(Vmec_Bmnlist_m))
   Vmec_Bmnlist_data = (wout['bmnc'][...])[skrindh,:]
#   Vmec_Bmn = mnmat(Vmec_Bmnlist, Vmec_Nvmecu, Vmec_Nvmecw, 'forceSize')
#   Vmec_B = ifftmn(Vmec_Bmn, Geom_Nperiods, Vmec_Nvmecu, Vmec_Nvmecw, 'forceSize')
   Vmec_B = np.zeros(np.shape(Vmec_vmecu))

   for i in range(0, len(Vmec_Bmnlist_data)):
      Vmec_B += Vmec_Bmnlist_data[i] * np.cos(Vmec_Bmnlist_m[i]*Vmec_vmecu - Vmec_Bmnlist_n[i]*Geom_Nperiods*Vmec_vmecw)
   
   Vmec_Jmnlist_m = np.float_(wout['xm_nyq'][...])
   Vmec_Jmnlist_n = signchange*np.float_(wout['xn_nyq'][...])/Geom_Nperiods
   Vmec_Jmnlist_cosparity = np.ones(np.shape(Vmec_Jmnlist_m))
   Vmec_Jmnlist_data = (wout['gmnc'][...])[skrindh,:] * signchange
#   Vmec_Jmn = mnmat(Vmec_Jmnlist, Vmec_Nvmecu, Vmec_Nvmecw, 'forceSize')
#   Vmec_J = ifftmn(Vmec_Jmn, Geom_Nperiods, Vmec_Nvmecu, Vmec_Nvmecw, 'forceSize')
   Vmec_J = np.zeros(np.shape(Vmec_vmecu))

   for i in range(0, len(Vmec_Jmnlist_data)):
      Vmec_J += Vmec_Jmnlist_data[i] * np.cos(Vmec_Jmnlist_m[i]*Vmec_vmecu - Vmec_Jmnlist_n[i]*Geom_Nperiods*Vmec_vmecw)   

#   Vmec_B_ulist_m = np.float_(wout.xm_nyq)
#   Vmec_B_ulist_n = signchange*np.float_(wout.xn_nyq)/Geom_Nperiods
#   Vmec_B_ulist_cosparity = np.ones(np.shape(Vmec_B_ulist_m))
#   Vmec_B_ulist_data = wout.bsubumnc[:,skrindh]
#   Vmec_B_umn = mnmat(Vmec_B_ulist, Vmec_Nvmecu, Vmec_Nvmecw, 'forceSize')
#   Vmec_B_umntilde = remove00(Vmec_B_umn)                                        
   
#   Vmec_B_wlist_m = np.float_(wout.xm_nyq)
#   Vmec_B_wlist_n = signchange*np.float_(wout.xn_nyq)/Geom_Nperiods
#   Vmec_B_wlist_cosparity = np.ones(np.shape(Vmec_B_wlist_m))
#   Vmec_B_wlist_data = wout.bsubvmnc[:,skrindh] * signchange
#   Vmec_B_wmn = mnmat(Vmec_B_wlist, Vmec_Nvmecu, Vmec_Nvmecw, 'forceSize')
#   Vmec_B_wmntilde = remove00(Vmec.B_wmn)
   
   Vmec_llist_m = np.float_(wout['xm'][...])
   Vmec_llist_n = signchange*np.float_(wout['xn'][...])/Geom_Nperiods
   Vmec_llist_cosparity = 0*np.ones(np.shape(Vmec_llist_m))
   Vmec_llist_data = (wout['lmns'][...])[skrindh,:]
#   Vmec_lmn = mnmat(Vmec_llist, Vmec_Nvmecu, Vmec_Nvmecw, 'forceSize')
#   Vmec_l = ifftmn(Vmec_lmn, Geom_Nperiods, Vmec_Nvmecu, Vmec_Nvmecw, 'forceSize')

   Vmec_l = np.zeros(np.shape(Vmec_vmecu))

   for i in range(0, len(Vmec_llist_data)):
      Vmec_l += Vmec_llist_data[i] * np.cos(Vmec_llist_m[i]*Vmec_vmecu - Vmec_llist_n[i]*Geom_Nperiods*Vmec_vmecw) 



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Make the Pest coordinates
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##   MOVED SECTION UP ##
#   Vmec_Dvmecw = 2*pi/Vmec_Nvmecw/Geom_Nperiods
#   Vmec_Dvmecu = 2*pi/Vmec_Nvmecu
#   Vmec_vmecuvec = np.arange(0,Vmec_Nvmecu)*Vmec_Dvmecu
#   Vmec_vmecwvec = (np.arange(0,Vmec_Nvmecw)).conj().transpose()*Vmec_Dvmecw
#   Vmec_vmecu, Vmec_vmecw = np.mgrid(Vmec_vmecuvec, Vmec_vmecwvec) ##Could be meshgrid() instead of mgrid()
########################

   #Vmec_Dpthetavmecu = ifftmn(Vmec_lmn, Geom_Nperiods, Vmec_Nvmecu, Vmec_Nvmecw, 'forceSize')
   Vmec_Dpthetavmecu = Vmec_l
   Vmec_ptheta = Vmec_Dpthetavmecu + Vmec_vmecu #This is the "transformation"
   Vmec_pzeta = Vmec_vmecw

   Pest_Nptheta = Vmec_Nvmecu
   Pest_Npzeta = Vmec_Nvmecw
   Pest_Dptheta = Vmec_Dvmecu
   Pest_Dpzeta = Vmec_Dvmecw
   Pest_pzeta = Vmec_vmecw
   Pest_ptheta = Vmec_vmecu

   ##ONLY NEED TO KEEP STUFF RELATED TO Vmec_l TO DO THE TRANSFORMATION

   #Pest_Dpthetavmecu = griddatacyclic(Vmec_ptheta, Vmec_pzeta, Vmec_Dpthetavmecu, Geom_Nperiods)
   Pest_Dpthetavmecu = griddatacyclic(Vmec_ptheta, Vmec_pzeta, Vmec_Dpthetavmecu, Geom_Nperiods)[0]
   Pest_vmecu = Pest_ptheta - Pest_Dpthetavmecu #This is the "transformation"
   Pest_vmecw = Pest_pzeta
   
   Pest_B = interp2_cyclic(Vmec_vmecu, Vmec_vmecw, Vmec_B, Pest_vmecu, Pest_vmecw, Geom_Nperiods)
   Pest_R = interp2_cyclic(Vmec_vmecu, Vmec_vmecw, Vmec_R, Pest_vmecu, Pest_vmecw, Geom_Nperiods)
   Pest_Z = interp2_cyclic(Vmec_vmecu, Vmec_vmecw, Vmec_Z, Pest_vmecu, Pest_vmecw, Geom_Nperiods)
#The routine interp2_cyclic can also be used to interpolate other quantities from
#vmec to pest or vice versa.


   #Pest_mnmat_B = fftmn(Pest_B) ##CHECK THIS FUNCTION
   #Pest_mnmat_R = fftmn(Pest_R)
   #Pest_mnmat_Z = fftmn(Pest_Z)
   
   Pest_R00 = np.mean(np.mean(Pest_R))
   Pest_B00 = np.mean(np.mean(Pest_B))
   


  
##   if 0:
##      figure(1)
##      #surf(Vmec.u,Vmec.w,Vmec.Dzetaw);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
##      #surf(Vmec.u,Vmec.w,Vmec.zeta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
## #!TEMPORARY COMMENTED      surf(Vmec.u,Vmec.w,Vmec.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
##
##      figure(2)
##      #surf(Vmec.u,Vmec.w,Vmec.Dthetau);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
## #!TEMPORARY COMMENTED      surf(Vmec.u,Vmec.w,Vmec.ptheta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
##      #surf(Vmec.u,Vmec.w,Vmec.u);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
##      
##      #fig(5)
##     #surf(Booz_u,Booz_w,Pest.Dzetaw);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
##      #surf(Booz_u,Booz_w,Pest.zeta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
##      #surf(Booz_u,Booz_w,Pest.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
##
##
##      figure(3)
## #!TEMPORARY COMMENTED      surf(Vmec.ptheta,Vmec.pzeta,Vmec.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
##      figure(7)
## #!TEMPORARY COMMENTED      surf(Pest.ptheta,Pest.pzeta,Pest.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])

      
   
   wout.close()
   os.remove(copyfile)
   #wout.close()
   sys.exit()


   #return (Pest, Vmec)
   #I AM HERE
   return (Pest_vmecu, Pest_vmecw, )


transformVMECtoPEST("/hydra/u/almo/Phi1/LHD/lhd2_B_III/Input/wout_lhd2.nc", 0.5, 10, 10, -1)
