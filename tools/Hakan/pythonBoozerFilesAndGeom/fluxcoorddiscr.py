#!/usr/bin/env python
from __future__ import division
import numpy as np
import sys
from scipy import interpolate

import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import geomlib
import mnFourierlib

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This class, fluxcoorddiscr, contains for one flux surface uniformly discretised flux coordinates
# (Boozer, Hamada or Pest) and physical and geometric quantities related to them as well as their
# relation to the other coordinates.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constructor function 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Typical usage:
#
# Booz   = fluxcoorddiscr(Geom,rind,Npol,Ntor,name='Boozer')
# Hamada = fluxcoorddiscr(Geom,rind,Npol,Ntor,name='Hamada')
# Pest   = fluxcoorddiscr(Geom,rind,Npol,Ntor,name='Pest')
#
# Inputs:
# Geom: An equilibrium loaded with the constructor to the bcgeom or vmecgeom classes.
#
# rind: Index of the chosen flux surface
#
# Npol and Ntor: The resolution of the spatial grid. Npol and Ntor must both be odd!
# Typical values are 101 - 151. Larger values make the calculation of Hamada coordinates heavy.
# Npol and Ntor should be chosen to resolve the necessary toroidal
# and poloidal mode numbers of the magnetic field.
#
# name: one of ('Boozer','Hamada','Pest'). The output is a uniform discretisation on these coordinates.
#
# Output:
# The outputs are the same in all three discretisations. 
# Names of coordinates (right-handed):
#
# Boozer:      psi,  theta,  zeta
# Hamada:      psi,  vthet,  phi   (Note: psi is the radial coordinate, not V)
# Pest:        psi,  ptheta, pzeta
# Jacobians are given as
# Jacob_psi_vthet_phi
# Jacob_psi_theta_zeta
# 
# Some differences between the coordinates are also given. They are denoted with a D
# followed by the names of the two coordinates. For instance, Dzetaphi is the difference
# between zeta and phi at the discretisation points. 
# 
# Some metric elements are also given. The notation is that, e.g.,
# g_thetazeta is g_{\theta\zeta} = \mathbf{e}_\theta \cdot \mathbf{e}_\zeta
# gpsipsi     is g^{\psi\psi} = \nabla\psi\cdot\nabla\psi
# 
# The quantity u satisfies
# (B dot grad) u = 2 iota |B|^-3 B x nabla psi dot nabla |B|
# 
# G and I are defined by the expression for the magnetic field
# \mathbf{B} = I\nabla\theta + G\nabla\zeta + B_\psi(\theta,\zeta)\nabla\psi
# \mathbf{B} = I\nabla\vthet + G\nabla\phi + \nabla H(\psi,\vthet,\phi)
# 
# where Bpsitilde is B_\psi up to an additive constant. 
# (The m=0,n=0 comp. of Booz.Bpsitilde is set to zero.)
# 
# Other quantities:
# h = 1/B^2
# FSAB2       = <B^2>   (FSA stands for flux surface average)
# FSAu2B2     = <u^2 B^2>
# FSAgpsipsi  = <g^{\psi\psi}>
# FSAg_phiphi = <g_{\phi\phi}>
# (R,Z,cylphi) or (R,-cylphi,Z): Coordinates of discretisation points in
#                                             the "major radius cylinder"
# (X,Y,Z) : Cartesian coordinates of discretisation points
# For the full list, use the function disp()
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# interp2FromOtherCoords() function 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Interpolate a quantity given on a uniform grid in another set of coordinates to
# this set.
#
# Typical usage (from Booz to Pest):
# Pest.interpFromOtherCoords(quant_in_boozcoord,'Booz')
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# disp() function 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Shows values of scalars and descriptions and datatypes for all elements
#
# Typical usage: Booz.disp()
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot() function 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Makes a 2d plot over the discretisation variables of an input quantity
# Typical usage:
# fig,plt=Booz.plot('gpsipsi',cmap='jet')                           shows a quantity from Booz
# fig,plt=Booz.plot(nparrayquantity,title='The quantity',cmap=None) shows some input quantity
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot3d() function 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Make a 3d plot of the flux surface with the input quantity shown as the color
# Typical usage:
# fig,plt=Booz.plot3d('B',title='Magnetic field',torstride=2,polstride=2,cmap=None)
# fig,plt=Booz.plot3d(nparrayquantity,title='The quantity',torstride=1,polstride=1,cmap=None)


class fluxcoorddiscr:

  @staticmethod
  def interp2_cyclic(uin,vin,F,ueval,veval,N):#general function to interpolate 2d cyclic on regular grids
      # u must be a 2*pi periodic column vector.
      # v must be a 2*pi/N periodic row vector.

      if uin.ndim==2 and vin.ndim==2:
            v=vin[0,:]
            u=uin[:,0]
      elif uin.ndim==1 and vin.ndim==1:
            v=vin
            u=uin
      else:
            sys.exit('Not foreseen case!')
            
      vDirection=np.sign(v[1]-v[0])
      uDirection=np.sign(u[1]-u[0])
      
      vbig=np.append(v,np.array([v[0]+vDirection*2*np.pi/N]),axis=0)
      ubig=np.append(u,np.array([u[0]+uDirection*2*np.pi]),axis=0)
      Fbig=np.append(F,np.array([F[:,0]]).T,axis=1)
      Fbig=np.append(Fbig,np.array([Fbig[0,:]]),axis=0)
      vmod=np.mod(veval,2*np.pi/N)
      umod=np.mod(ueval,2*np.pi)

      ip=interpolate.RectBivariateSpline(ubig,vbig,Fbig)
      shp=umod.shape

      return ip(umod.flatten(),vmod.flatten(),grid=False).reshape(shp)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Constructor function 
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  def __init__(self,Geom,rind,Npol,Ntor,name='Boozer',u_zeroout_Deltaiota=-1.0,extended=False,makeSmith=False):

      #A function for interpolating on a nonuniform grid
      def griddatacyclic(unoneq,vnoneq,Fnoneq,Nperiods):
          # u must be 2*pi periodic.
          # v must be 2*pi/N periodic.
          # (unoneq, vnoneq) is the non-equidistant grid on which F is given.
          # (u,v) is the equidistant grid generated with ndgrid which to interpolate to.
          Nu=Fnoneq.shape[0]
          Nv=Fnoneq.shape[1]
          Dv=2*np.pi/Nv/Nperiods
          Du=2*np.pi/Nu
          u, v = np.mgrid[0.0:2.0*np.pi-Du:1j*Nu,0.0:2.0*np.pi/Nperiods-Dv:1j*Nv]

          Nleftadd=np.int(np.ceil(max(vnoneq[:,0])/(2*np.pi/Nperiods)))
          Nrightadd=np.int(np.ceil(1-min(vnoneq[:,-1])/(2*np.pi/Nperiods)))
              
          tmp=np.zeros((Nu*3,Nv*(Nleftadd+1+Nrightadd)))
          vnoneqBig      =np.nan*np.zeros((Nu*3,Nv*(Nleftadd+1+Nrightadd)))
          unoneqBig      =np.nan*np.zeros((Nu*3,Nv*(Nleftadd+1+Nrightadd)))
          FnoneqBig =np.nan*np.zeros((Nu*3,Nv*(Nleftadd+1+Nrightadd)))

          for addi in range(Nleftadd+1+Nrightadd):
              addci=addi-Nleftadd #is =0 for the original
              toaddtor=2*np.pi/Nperiods*addci
              vnoneqBig[   0:1*Nu, addi*Nv:(addi+1)*Nv]=vnoneq+toaddtor
              vnoneqBig[  Nu:2*Nu, addi*Nv:(addi+1)*Nv]=vnoneq+toaddtor
              vnoneqBig[2*Nu:3*Nu, addi*Nv:(addi+1)*Nv]=vnoneq+toaddtor
              unoneqBig[   0:1*Nu, addi*Nv:(addi+1)*Nv]=unoneq-2*np.pi
              unoneqBig[  Nu:2*Nu, addi*Nv:(addi+1)*Nv]=unoneq
              unoneqBig[2*Nu:3*Nu, addi*Nv:(addi+1)*Nv]=unoneq+2*np.pi
              FnoneqBig[   0:1*Nu, addi*Nv:(addi+1)*Nv]=Fnoneq
              FnoneqBig[  Nu:2*Nu, addi*Nv:(addi+1)*Nv]=Fnoneq
              FnoneqBig[2*Nu:3*Nu, addi*Nv:(addi+1)*Nv]=Fnoneq

          F=interpolate.griddata(np.concatenate((np.expand_dims(vnoneqBig.flatten(),axis=0).T,
                                                 np.expand_dims(unoneqBig.flatten(),axis=0).T),axis=1),
                                                FnoneqBig.flatten(),(v,u))
          return F

      ##############################################################
      #The Constructor starts here:
      ##############################################################   
      
      self.name=name
      if not((name=='Boozer') or (name=='Hamada') or (name=='Pest') or (name=='Smith')):
        sys.exit('The coordinate name must be one of Boozer, Hamada, Pest or Smith!')

      Ntheta=Npol
      Nzeta=Ntor
      
      mu0=1.256637061e-06

      if extended:
        ext_B = [None]*3
        ext_R = [None]*3
        ext_Z = [None]*3
        ext_X = [None]*3
        ext_Y = [None]*3
        ext_Dzetacylphi = [None]*3
        ext_Dzetaphi = [None]*3
        if rind==0:
          ext_rindstart=0
          ext_relrind=0
        elif rind==Geom.nsurf-1:
          ext_rindstart=Geom.nsurf-3
          ext_relrind=2
        else:
          ext_rindstart=rind-1
          ext_relrind=1
          
      if isinstance(Geom,geomlib.bcgeom):
          if rind > Geom.nsurf-1:
              sys.exit('Error in input to fluxcoorddiscr.py: The radius index rind = '+
                       str(rind)+' is greater than max rind = '+str(Geom.nsurf-1))
          Nperiods=Geom.Nperiods
          I=Geom.Btheta[rind]
          G=Geom.Bphi[rind]
          iota=Geom.iota[rind]
          B00=Geom.B00[rind]

          NHarmonics=Geom.nmodes[rind]
          mu0dpdpsi=Geom.dpds[rind]/Geom.psi_a*mu0
          psi_a=Geom.psi_a

          #Boozer discretisation
          Dtheta=2*np.pi/Ntheta
          Dzeta=2*np.pi/Nzeta/Geom.Nperiods
          theta, zeta = np.mgrid[0.0:2.0*np.pi-Dtheta:1j*Ntheta,0.0:2.0*np.pi/Geom.Nperiods-Dzeta:1j*Nzeta]

          Bmn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='B')
          Rmn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='R')
          Zmn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='Z')
          Dzetacylphimn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='Dphi')
          Dzetacylphi=Dzetacylphimn.ifft()
          type(Bmn)
          B=Bmn.ifft()
          R=Rmn.ifft()
          Z=Zmn.ifft()
          if extended:
            ext_s   =Geom.s[ext_rindstart:ext_rindstart+3]
            ext_G   =Geom.Bphi[ext_rindstart:ext_rindstart+3]
            ext_I   =Geom.Btheta[ext_rindstart:ext_rindstart+3]
            ext_iota=Geom.iota[ext_rindstart:ext_rindstart+3]
            for ri in range(3):
              ext_B[ri]=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=ext_rindstart+ri,quantity='B').ifft()
              ext_R[ri]=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=ext_rindstart+ri,quantity='R').ifft()
              ext_Z[ri]=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=ext_rindstart+ri,quantity='Z').ifft()
              ext_Dzetacylphi[ri]=(
                      mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=ext_rindstart+ri,quantity='Dphi').ifft())
              loc_cylphi=zeta-ext_Dzetacylphi[ri] #cylphi is minus the geometrical toroidal angle.
              loc_geomang=-loc_cylphi         #(R,Z,cylphi) and (R,geomang,Z) are right handed systems. 
              ext_X[ri]=ext_R[ri]*np.cos(loc_geomang)
              ext_Y[ri]=ext_R[ri]*np.sin(loc_geomang)

            
      elif isinstance(Geom,geomlib.vmecgeom):
          #VMEC LH coords: (u,v,s). Let w=-v => RH coords (u,w,s)
          #PEST coords pzeta=w=cylphi, ptheta=u+lambda. Here lambda is stored in the VMEC file.

          signchange=float(Geom.signgs) #is -1, because vmec is left handed
          skip=Geom.skip #=1,this is how many elements are skipped at low radii when going to half grid
          if rind > Geom.ns-1-skip: #Note: max(rind) on the half grid is Geom.ns-2
              sys.exit('Error in input to fluxcoorddiscr.py: The radius index rind = '+
                       str(rind)+' is greater than max rind = '+str(Geom.nsurf-1-skip))
          fullgrid_s      = Geom.phi/Geom.phi[Geom.ns-1] #full grid
          halfgrid_s      = (fullgrid_s[1:]+fullgrid_s[0:-1])/2
          halfgrid_iota   = Geom.iotas[skip:]*signchange
          halfgrid_Bphi   = Geom.bvco[skip:]*signchange
          halfgrid_Btheta = Geom.buco[skip:]
          halfgrid_dpds   = np.diff(Geom.presf[skip-1:])/np.diff(fullgrid_s[skip-1:])
          psi_a=Geom.phi[Geom.ns-1]*signchange/2.0/np.pi
          
          Nperiods = Geom.nfp
          iota=halfgrid_iota[rind]
          G=halfgrid_Bphi[rind]
          I=halfgrid_Btheta[rind]
          mu0dpdpsi=halfgrid_dpds[rind]/psi_a*mu0

          uw_Bmn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='B')
          uw_Rmn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='R')
          uw_Zmn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='Z')
          uw_B_umntilde=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='B_u').remove00()
          uw_B_wmntilde=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='B_w').remove00()
          uw_lmn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='lambda')
          #print uw_lmn.ifft()
          #alphamn1=uw_B_umntilde.invgrad(uw_B_wmntilde,method=1)
          alphamn2=uw_B_umntilde.invgrad(uw_B_wmntilde,method=2)
          #if B_umn and B_wmn are good, then alphamn1 should equal alphamn2 (any problem here would be due to vmec)
          
          uw_pmn=(alphamn2-I*uw_lmn)*(1.0/(G+iota*I))
          uw_Dthetau=(uw_lmn+iota*uw_pmn).ifft()
          uw_Dzetaw =uw_pmn.ifft()
          uw_Dpthetau=uw_lmn.ifft()
      
          #u,w discretisation
          Dtheta=2.0*np.pi/Ntheta
          Dzeta=2.0*np.pi/Nzeta/Nperiods
          uniform_pol,uniform_tor = np.mgrid[0.0:2.0*np.pi-Dtheta:1j*Ntheta,0.0:2.0*np.pi/Nperiods-Dzeta:1j*Nzeta]
          uw_u = uniform_pol
          uw_w = uniform_tor
          uw_theta=uw_Dthetau+uw_u
          uw_zeta =uw_Dzetaw +uw_w
          #uw_pzeta=uw_w
          #uw_ptheta=uw_u+uw_Dpthetau
          #print 'uw_Dthetau='+str(uw_Dthetau)
          #print 'uw_Dzetaw='+str(uw_Dthetau)
          #Pest discretisation
          #Pest_ptheta=uniform_pol
          #Pest_pzeta =uniform_tor
          theta=uniform_pol
          zeta =uniform_tor

          Booz_Dthetau =griddatacyclic(uw_theta,uw_zeta,uw_Dthetau,Nperiods)
          Booz_Dzetaw  =griddatacyclic(uw_theta,uw_zeta,uw_Dzetaw, Nperiods)
          Dzetacylphi=Booz_Dzetaw
          #Pest_Dpthetau=griddatacyclic(uw_ptheta,uw_pzeta,uw_Dpthetau,self.Nperiods)
    
          Booz_u=theta-Booz_Dthetau        
          Booz_w=zeta-Booz_Dzetaw
    
          #Pest_u=Pest_ptheta-Pest_Dpthetau
          #Pest_w=Pest_pzeta
    
          uw_B=uw_Bmn.ifft()
          uw_R=uw_Rmn.ifft()
          uw_Z=uw_Zmn.ifft()
          
          B=fluxcoorddiscr.interp2_cyclic(uw_u,uw_w,uw_B,Booz_u,Booz_w,Nperiods)
          R=fluxcoorddiscr.interp2_cyclic(uw_u,uw_w,uw_R,Booz_u,Booz_w,Nperiods)
          Z=fluxcoorddiscr.interp2_cyclic(uw_u,uw_w,uw_Z,Booz_u,Booz_w,Nperiods)
    
          Bmn=mnFourierlib.mnmat(B,Nperiods=Nperiods)
          Rmn=mnFourierlib.mnmat(R,Nperiods=Nperiods)
          Zmn=mnFourierlib.mnmat(Z,Nperiods=Nperiods)
          Dzetacylphimn=mnFourierlib.mnmat(Dzetacylphi,Nperiods=Nperiods)

          B00=Bmn.get00()

          if extended:
            ext_s   =halfgrid_s[ext_rindstart:ext_rindstart+3]
            ext_iota=halfgrid_iota[ext_rindstart:ext_rindstart+3]
            ext_G   =halfgrid_Bphi[ext_rindstart:ext_rindstart+3]
            ext_I   =halfgrid_Btheta[ext_rindstart:ext_rindstart+3]
            for ri in range(3):
              locrind=ext_rindstart+ri
              loc_iota=halfgrid_iota[locrind]
              loc_G=halfgrid_Bphi[locrind]
              loc_I=halfgrid_Btheta[locrind]
              loc_mu0dpdpsi=halfgrid_dpds[locrind]/psi_a*mu0

              loc_uw_Bmn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=locrind,quantity='B')
              loc_uw_Rmn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=locrind,quantity='R')
              loc_uw_Zmn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=locrind,quantity='Z')
              loc_uw_B_umntilde=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=locrind,quantity='B_u').remove00()
              loc_uw_B_wmntilde=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=locrind,quantity='B_w').remove00()
              loc_uw_lmn=mnFourierlib.mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=locrind,quantity='lambda')
              #print uw_lmn.ifft()
              #alphamn1=uw_B_umntilde.invgrad(uw_B_wmntilde,method=1)
              loc_alphamn2=loc_uw_B_umntilde.invgrad(loc_uw_B_wmntilde,method=2)
              #if B_umn and B_wmn are good, then alphamn1 should equal alphamn2 (any problem here would be due to vmec)

              loc_uw_pmn=(loc_alphamn2-loc_I*loc_uw_lmn)*(1.0/(loc_G+loc_iota*loc_I))
              loc_uw_Dthetau=(loc_uw_lmn+loc_iota*loc_uw_pmn).ifft()
              loc_uw_Dzetaw =loc_uw_pmn.ifft()
              loc_uw_Dpthetau=loc_uw_lmn.ifft()
            
              loc_uw_theta=loc_uw_Dthetau+uw_u
              loc_uw_zeta =loc_uw_Dzetaw +uw_w

              loc_Booz_Dthetau =griddatacyclic(loc_uw_theta,loc_uw_zeta,loc_uw_Dthetau,Nperiods)
              loc_Booz_Dzetaw  =griddatacyclic(loc_uw_theta,loc_uw_zeta,loc_uw_Dzetaw, Nperiods)
              ext_Dzetacylphi[ri]  =loc_Booz_Dzetaw
                  
              loc_Booz_u=theta-loc_Booz_Dthetau        
              loc_Booz_w=zeta-loc_Booz_Dzetaw
        
              loc_uw_B=loc_uw_Bmn.ifft()
              loc_uw_R=loc_uw_Rmn.ifft()
              loc_uw_Z=loc_uw_Zmn.ifft()
          
              ext_B[ri]=interp2_cyclic(uw_u,uw_w,loc_uw_B,loc_Booz_u,loc_Booz_w,Nperiods)
              ext_R[ri]=interp2_cyclic(uw_u,uw_w,loc_uw_R,loc_Booz_u,loc_Booz_w,Nperiods)
              ext_Z[ri]=interp2_cyclic(uw_u,uw_w,loc_uw_Z,loc_Booz_u,loc_Booz_w,Nperiods)
    
              #ext_Bmn[ri]=mnFourierlib.mnmat(ext_B[ri],Nperiods=Nperiods)
              #ext_Rmn[ri]=mnFourierlib.mnmat(ext_R[ri],Nperiods=Nperiods)
              #ext_Zmn[ri]=mnFourierlib.mnmat(ext_Z[ri],Nperiods=Nperiods)
              #ext_Dzetacylphimn[ri]=mnFourierlib.mnmat(ext_Dzetacylphi[ri],Nperiods=Nperiods)
              loc_cylphi=zeta-ext_Dzetacylphi[ri] #cylphi is minus the geometrical toroidal angle.
              loc_geomang=-loc_cylphi         #(R,Z,cylphi) and (R,geomang,Z) are right handed systems. 
              ext_X[ri]=ext_R[ri]*np.cos(loc_geomang)
              ext_Y[ri]=ext_R[ri]*np.sin(loc_geomang)
              
          #end if extended
      #end if bcgeom or vmecgeom
      
      self.Nperiods=Nperiods
      self.iota=iota
      self.G=G
      self.I=I
      self.mu0dpdpsi=mu0dpdpsi

      dBdthetamn,dBdzetamn=Bmn.grad()
      dRdthetamn,dRdzetamn=Rmn.grad()
      dZdthetamn,dZdzetamn=Zmn.grad()
      dDzetacylphidthetamn,dDzetacylphidzetamn=Dzetacylphimn.grad()

      d2Rdtheta2mn,    d2Rdthetadzetamn=dRdthetamn.grad()
      d2Rdthetadzetamn,d2Rdzeta2mn     =dRdzetamn.grad()
      d2Zdtheta2mn,    d2Zdthetadzetamn=dZdthetamn.grad()
      d2Zdthetadzetamn,d2Zdzeta2mn     =dZdzetamn.grad()
      d2Dzetacylphidtheta2mn,    d2Dzetacylphidthetadzetamn=dDzetacylphidthetamn.grad()
      d2Dzetacylphidthetadzetamn,d2Dzetacylphidzeta2mn     =dDzetacylphidzetamn.grad()

      dBdtheta=dBdthetamn.ifft()
      dBdzeta =dBdzetamn.ifft()
      dRdtheta=dRdthetamn.ifft()
      dRdzeta =dRdzetamn.ifft()
      dZdtheta=dZdthetamn.ifft()
      dZdzeta =dZdzetamn.ifft()
 
      dDzetacylphidtheta=dDzetacylphidthetamn.ifft()
      dDzetacylphidzeta =dDzetacylphidzetamn.ifft()

      d2Rdtheta2              = d2Rdtheta2mn.ifft()
      d2Rdthetadzeta          = d2Rdthetadzetamn.ifft()
      d2Rdzeta2               = d2Rdzeta2mn.ifft()
      d2Zdtheta2              = d2Zdtheta2mn.ifft()
      d2Zdthetadzeta          = d2Zdthetadzetamn.ifft()
      d2Zdzeta2               = d2Zdzeta2mn.ifft()
      d2Dzetacylphidtheta2    = d2Dzetacylphidtheta2mn.ifft()
      d2Dzetacylphidthetadzeta= d2Dzetacylphidthetadzetamn.ifft()
      d2Dzetacylphidzeta2     = d2Dzetacylphidzeta2mn.ifft()

      cylphi=zeta-Dzetacylphi #cylphi is minus the geometrical toroidal angle.
      geomang=-cylphi         #(R,Z,cylphi) and (R,geomang,Z) are right handed systems. 

      dgeomangdtheta      = dDzetacylphidtheta
      dgeomangdzeta       = dDzetacylphidzeta - 1
      d2geomangdtheta2    = d2Dzetacylphidtheta2
      d2geomangdzeta2     = d2Dzetacylphidzeta2
      d2geomangdthetadzeta= d2Dzetacylphidthetadzeta
      X=R*np.cos(geomang)
      Y=R*np.sin(geomang)

      dXdtheta=dRdtheta*np.cos(geomang)-R*dgeomangdtheta*np.sin(geomang)
      dXdzeta =dRdzeta *np.cos(geomang)-R*dgeomangdzeta *np.sin(geomang) 
      dYdtheta=dRdtheta*np.sin(geomang)+R*dgeomangdtheta*np.cos(geomang)
      dYdzeta =dRdzeta *np.sin(geomang)+R*dgeomangdzeta *np.cos(geomang)

      d2Xdtheta2=(d2Rdtheta2*np.cos(geomang) 
        -2*dRdtheta*dgeomangdtheta*np.sin(geomang)
        -R*d2geomangdtheta2*np.sin(geomang) 
        -R*dgeomangdtheta**2*np.cos(geomang))
      d2Xdthetadzeta=(d2Rdthetadzeta*np.cos(geomang) 
        -(dRdtheta*dgeomangdzeta+dRdzeta*dgeomangdtheta)*np.sin(geomang) 
        -R*d2geomangdthetadzeta*np.sin(geomang) 
        -R*dgeomangdtheta*dgeomangdzeta*np.cos(geomang))
      d2Xdzeta2=(d2Rdzeta2*np.cos(geomang) 
        -2*dRdzeta*dgeomangdzeta*np.sin(geomang)
        -R*d2geomangdzeta2*np.sin(geomang) 
        -R*dgeomangdzeta**2*np.cos(geomang))

      d2Ydtheta2=(d2Rdtheta2*np.sin(geomang) 
        +2*dRdtheta*dgeomangdtheta*np.cos(geomang)
        +R*d2geomangdtheta2*np.cos(geomang) 
        -R*dgeomangdtheta**2*np.sin(geomang))
      d2Ydthetadzeta=(d2Rdthetadzeta*np.sin(geomang) 
        +(dRdtheta*dgeomangdzeta+dRdzeta*dgeomangdtheta)*np.cos(geomang) 
        +R*d2geomangdthetadzeta*np.cos(geomang) 
        -R*dgeomangdtheta*dgeomangdzeta*np.sin(geomang))
      d2Ydzeta2=(d2Rdzeta2*np.sin(geomang) 
        +2*dRdzeta*dgeomangdzeta*np.cos(geomang)
        +R*d2geomangdzeta2*np.cos(geomang) 
        -R*dgeomangdzeta**2*np.sin(geomang))

      g_thetatheta=dXdtheta**2+dYdtheta**2+dZdtheta**2
      g_zetazeta  =dXdzeta**2 +dYdzeta**2 +dZdzeta**2
      g_thetazeta =dXdtheta*dXdzeta+dYdtheta*dYdzeta+dZdtheta*dZdzeta

      gradpsi_X=B**2/(G+iota*I)*(dYdtheta*dZdzeta-dZdtheta*dYdzeta)
      gradpsi_Y=B**2/(G+iota*I)*(dZdtheta*dXdzeta-dXdtheta*dZdzeta)
      gradpsi_Z=B**2/(G+iota*I)*(dXdtheta*dYdzeta-dYdtheta*dXdzeta)
      gpsipsi=gradpsi_X**2+gradpsi_Y**2+gradpsi_Z**2

      Booz_XYZ_gradpsi=np.zeros((Ntheta,Nzeta,3))
      Booz_XYZ_gradpsi[:,:,0]=gradpsi_X
      Booz_XYZ_gradpsi[:,:,1]=gradpsi_Y
      Booz_XYZ_gradpsi[:,:,2]=gradpsi_Z

      Booz_XYZ_e_theta=np.zeros((Ntheta,Nzeta,3))
      Booz_XYZ_e_theta[:,:,0]=dXdtheta
      Booz_XYZ_e_theta[:,:,1]=dYdtheta
      Booz_XYZ_e_theta[:,:,2]=dZdtheta
      
      Booz_XYZ_e_zeta=np.zeros((Ntheta,Nzeta,3))
      Booz_XYZ_e_zeta[:,:,0]=dXdzeta
      Booz_XYZ_e_zeta[:,:,1]=dYdzeta
      Booz_XYZ_e_zeta[:,:,2]=dZdzeta
      
      Booz_XYZ_d2rdtheta2=np.zeros((Ntheta,Nzeta,3))
      Booz_XYZ_d2rdtheta2[:,:,0]=d2Xdtheta2
      Booz_XYZ_d2rdtheta2[:,:,1]=d2Ydtheta2
      Booz_XYZ_d2rdtheta2[:,:,2]=d2Zdtheta2
      
      Booz_XYZ_d2rdzeta2=np.zeros((Ntheta,Nzeta,3))
      Booz_XYZ_d2rdzeta2[:,:,0]=d2Xdzeta2
      Booz_XYZ_d2rdzeta2[:,:,1]=d2Ydzeta2
      Booz_XYZ_d2rdzeta2[:,:,2]=d2Zdzeta2
      
      Booz_XYZ_d2rdthetadzeta=np.zeros((Ntheta,Nzeta,3))
      Booz_XYZ_d2rdthetadzeta[:,:,0]=d2Xdthetadzeta
      Booz_XYZ_d2rdthetadzeta[:,:,1]=d2Ydthetadzeta
      Booz_XYZ_d2rdthetadzeta[:,:,2]=d2Zdthetadzeta

      Booz_CX=(d2Xdzeta2+2*iota*d2Xdthetadzeta+iota**2*d2Xdtheta2)*(B**2/(G+iota*I))**2;
      Booz_CY=(d2Ydzeta2+2*iota*d2Ydthetadzeta+iota**2*d2Ydtheta2)*(B**2/(G+iota*I))**2;
      Booz_CZ=(d2Zdzeta2+2*iota*d2Zdthetadzeta+iota**2*d2Zdtheta2)*(B**2/(G+iota*I))**2;
      
      Booz_BdotgradabsB=B**2/(G+iota*I)*(iota*dBdtheta+dBdzeta);
      Booz_BxgradpsidotgradabsB=B**2/(G+iota*I)*(G*dBdtheta-I*dBdzeta);

      #Direction of B
      Booz_XYZ_B=np.zeros((Ntheta,Nzeta,3))
      Booz_XYZ_B[:,:,0]=(Booz_XYZ_e_zeta[:,:,0]+iota*Booz_XYZ_e_theta[:,:,0])*B**2/(G+iota*I)
      Booz_XYZ_B[:,:,1]=(Booz_XYZ_e_zeta[:,:,1]+iota*Booz_XYZ_e_theta[:,:,1])*B**2/(G+iota*I)
      Booz_XYZ_B[:,:,2]=(Booz_XYZ_e_zeta[:,:,2]+iota*Booz_XYZ_e_theta[:,:,2])*B**2/(G+iota*I)

      Booz_BX=Booz_XYZ_B[:,:,0]
      Booz_BY=Booz_XYZ_B[:,:,1]
      Booz_BZ=Booz_XYZ_B[:,:,2]
      Booz_BR=Booz_BX*np.cos(geomang)+Booz_BY*np.sin(geomang)
      Booz_Bgeomang=-Booz_BX*np.sin(geomang)+Booz_BY*np.cos(geomang)

      Booz_XYZ_curv=np.zeros((Ntheta,Nzeta,3))
      Booz_XYZ_curv[:,:,0]=(Booz_CX+Booz_BX/B*Booz_BdotgradabsB)/B**2
      Booz_XYZ_curv[:,:,1]=(Booz_CY+Booz_BY/B*Booz_BdotgradabsB)/B**2
      Booz_XYZ_curv[:,:,2]=(Booz_CZ+Booz_BZ/B*Booz_BdotgradabsB)/B**2

      Booz_curv_normal=np.sum(Booz_XYZ_gradpsi*Booz_XYZ_curv,axis=2)/np.sqrt(gpsipsi)
      Booz_curv_geodes1=Booz_BxgradpsidotgradabsB/np.sqrt(gpsipsi)/B**2
      Booz_curv_geodes2=np.sum(np.cross(Booz_XYZ_B,Booz_XYZ_gradpsi)*Booz_XYZ_curv,axis=2)/np.sqrt(gpsipsi)/B
      Booz_gradpsidotgradB=Booz_curv_normal*np.sqrt(gpsipsi)*B-mu0dpdpsi*gpsipsi
    
      #% The following is for equilibria close to axi-symmetry, where we want to calculate 
      #% the deviation from axi-symmetry in terms of dBinsurf and dBperp.
      #% Approximating that dr/dzeta is roughly in the geometrical toroidal direction I 
      #% assume that poloidal components of B in the "closest" axisymmetric equilibrium can
      #% be obtained by averaging BR and BZ over the Boozer toroidal angle. Then we project
      #% the deviation from axi-symmetry on the two directions perpedicular to the toroidal direction.
      BRmean=np.mean(Booz_BR, axis=1, keepdims=True) #mean over the Boozer toroidal angle
      BZmean=np.mean(Booz_BZ, axis=1, keepdims=True) #mean over the Boozer toroidal angle

      bRmean=BRmean/np.sqrt(BRmean**2+BZmean**2)  #unit vector tangent to flux surface
      bZmean=BZmean/np.sqrt(BRmean**2+BZmean**2)  #in the poloidal plane
      
      dBR=Booz_BR-BRmean
      dBZ=Booz_BZ-BZmean

      Booz_dBinsurf=dBR*bRmean+dBZ*bZmean #dot product
      Booz_dBperp=dBR*bZmean-dBZ*bRmean   #cross product 

      #% ---------------------------------------------------------------------------------------
      #% Calculate parallel current u from harmonics of 1/B**2. 
      #% \nabla_\parallel u = (2/B^4) \nabla B \times \vector{B} \cdot \iota \nabla \psi 
      #% ---------------------------------------------------------------------------------------
      u = np.zeros((Ntheta,Nzeta))
      #Dzetaphi        =np.zeros(Ntheta,Nzeta) #difference between zeta (tor Boozer coord)
                                               #and phi (tor Hamada coord) 
      h=1/(B**2)
      VPrimeHat=np.sum(h)*4*np.pi**2/(Nzeta*Ntheta) #Note: VPrime=VPrimeHat*(G+iota*I)
      FSAB2=4*np.pi**2/VPrimeHat
      #h00=VPrimeHat/(4*pi**2);

      umn=mnFourierlib.mnmat(h,Nperiods=Nperiods).calcu(G,I,iota,zeroout_Deltaiota=u_zeroout_Deltaiota) 
      u=umn.ifft()

      Dzetaphimn=mnFourierlib.mnmat(1-h*FSAB2,Nperiods=Nperiods).invJacBdotgrad(iota)
      Dzetaphi=Dzetaphimn.ifft()
      
      dDzetaphidthetamn,dDzetaphidzetamn=Dzetaphimn.grad()
      dDzetaphidtheta=dDzetaphidthetamn.ifft()
      dDzetaphidzeta=dDzetaphidzetamn.ifft()
      
      Booz_B_psi_tilde=-mu0dpdpsi/FSAB2*(G+iota*I)*Dzetaphi
      Booz_dBpsidtheta=-mu0dpdpsi/iota*(u-iota*I*(B**-2 - 1/FSAB2))
      Booz_dBpsidzeta = mu0dpdpsi/iota*(u+G*(B**-2 - 1/FSAB2))
      Booz_Jacob_psi_theta_zeta = h*(G+iota*I)

      Jacob_psi_vthet_phi  = (G+iota*I)/FSAB2
      
      XtozmXzot=dXdtheta*dDzetaphidzeta-dXdzeta*dDzetaphidtheta
      YtozmYzot=dYdtheta*dDzetaphidzeta-dYdzeta*dDzetaphidtheta
      ZtozmZzot=dZdtheta*dDzetaphidzeta-dZdzeta*dDzetaphidtheta

      dXdphi=B**2/FSAB2*(dXdzeta+iota*XtozmXzot)
      dYdphi=B**2/FSAB2*(dYdzeta+iota*YtozmYzot)
      dZdphi=B**2/FSAB2*(dZdzeta+iota*ZtozmZzot)

      dXdvthet=B**2/FSAB2*(dXdtheta-XtozmXzot)
      dYdvthet=B**2/FSAB2*(dYdtheta-YtozmYzot)
      dZdvthet=B**2/FSAB2*(dZdtheta-ZtozmZzot)

      g_phiphi     = dXdphi**2 + dYdphi**2 + dZdphi**2
      g_vthetvthet = dXdvthet**2 + dYdvthet**2 + dZdvthet**2
      g_vthetphi   = dXdphi*dXdvthet + dYdphi*dYdvthet + dZdphi*dZdvthet

      Booz_Dzetaphi=Dzetaphi
      Booz_Dzetapzeta=Dzetacylphi
      
      Booz_phi   =zeta-Dzetaphi
      Booz_vthet =theta-iota*Dzetaphi
      Booz_pzeta =zeta-Dzetacylphi       #pzeta=cylphi
      Booz_ptheta=theta-iota*Dzetacylphi 

      dDzetapzetadthetamn,dDzetapzetadzetamn=mnFourierlib.mnmat(Booz_Dzetapzeta,Nperiods=Nperiods).grad()
      Booz_Jacob_psi_ptheta_pzeta=Booz_Jacob_psi_theta_zeta*(1+iota*dDzetapzetadthetamn.ifft()+dDzetapzetadzetamn.ifft())
      
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      #% Calculate Smith coordinates (spol,stor) 
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if name=='Smith':
        makeSmith=True
      if makeSmith:
        zetjumptolerance=np.pi/Nperiods*0.1
        Dzetaeta=np.zeros_like(Dzetaphi)
        zetv0=np.linspace(-1.0,1.0,num=Nzeta,endpoint=False)*np.pi/Nperiods
        thev=iota*zetv0
        zetguess=zetv0[np.argmax(Bmn.evalpoint(thev,zetv0))]
        #(th0,zetguess) is the starting point for the max search along B
        th0s=np.linspace(0.0,2.0*np.pi,num=max(Ntheta,111),endpoint=False)
        ths=np.nan*np.zeros_like(th0s)
        zes=np.nan*np.zeros_like(th0s)
        absd2Bdz2_alongB=np.nan*np.zeros_like(th0s)
        th0step=th0s[2]-th0s[1]
        shortzetv0=np.linspace(-1.0,1.0,num=11,endpoint=False)*min(2.0*Dzeta,np.pi/Nperiods*0.1)
        for th0i in range(len(th0s)):
            print('th0i='+str(th0i))
            shortzetv=shortzetv0+zetguess
            shortthev=th0s[th0i]+iota*shortzetv
            dBdzetv=iota*dBdthetamn.evalpoint(shortthev,shortzetv)+dBdzetamn.evalpoint(shortthev,shortzetv)
            Li=np.where(np.diff(np.sign(dBdzetv))<0.0)[0][0]
            zes[th0i]=shortzetv[Li]+(shortzetv[Li+1]-shortzetv[Li])*dBdzetv[Li]/(dBdzetv[Li]-dBdzetv[Li+1])
            ths[th0i]=th0s[th0i]+iota*zes[th0i]
            absd2Bdz2_alongB[th0i]=-(dBdzetv[Li+1]-dBdzetv[Li])/(shortzetv[Li+1]-shortzetv[Li])
            if th0i<len(th0s)-1:
              #find next guess
              zes_extrap=zes[th0i]
              #if th0i>0:
              #    zes_extrap=2.0*zes[th0i]-zes[th0i-1]
              zetv=zetv0+zes_extrap
              thev=th0s[th0i+1]+iota*zetv
              zetguess=zetv[np.argmax(Bmn.evalpoint(thev,zetv))]              
              
              if abs(zetguess-zes[th0i])>zetjumptolerance:
                  makeSmith=False
                  plt.plot(zetv0,th0s[th0i]+iota*zetv0,'r--',
                           zes,ths,'b-x',
                           zes_extrap,th0s[th0i+1]+iota*zes_extrap,'b+')
                  plt.show()
                  print('Smith not possible! Decrease th0step or try fewer mn modes in the equilibrium!')
        #end for th0i
      #end if makeSmith

      if not(makeSmith):
          if name=='Smith':
            sys.exit('Could not generate Smith coordinates.')
            
      if makeSmith:
        Bmaxs=Bmn.evalpoint(ths,zes)
        possibleshifts=np.array([-3,-2,1,0,1,2,3])*2.0*np.pi/Nperiods
        thezshift=possibleshifts[np.argmin(np.abs(zes[-1]-zes[0]-possibleshifts))]
        zes_app=np.append(zes,[zes[0]+thezshift])
        #ths_app=np.append(ths,[ths[0]])
        th0s_app=np.append(th0s,[2*np.pi])
        func2=interpolate.interp1d(th0s_app,zes_app)
        Smith_spol=theta #use a uniform grid
        Smith_stor=zeta  #use a uniform grid
        SmithDzetastor=func2(Smith_spol)
        SmithDzetastormn=mnFourierlib.mnmat(SmithDzetastor,Nperiods=Nperiods)
        isStelSym=not(np.any(np.abs(Bmn.s[np.where(np.logical_not(np.isnan(Bmn.s)))])>0.0))
        if isStelSym: #if Stellarator symmetric
          SmithDzetastormn.c[np.where(np.logical_not(np.isnan(Bmn.c)))]=0
          SmithDzetastor=SmithDzetastormn.ifft()
        SmithdDzetastordspolmn,dummy=SmithDzetastormn.grad()
        SmithdDzetastordspol=SmithdDzetastordspolmn.ifft()
        Smith_zeta=Smith_stor+SmithDzetastor
        Smith_theta=Smith_spol+iota*SmithDzetastor

        BoozDzetastor=griddatacyclic(Smith_theta,Smith_zeta,SmithDzetastor,self.Nperiods)
        Booz_dDzetastordspol=griddatacyclic(Smith_theta,Smith_zeta,SmithdDzetastordspol,self.Nperiods)
        Booz_spol=theta-iota*BoozDzetastor
        Booz_stor=zeta-BoozDzetastor
        Booz_Jacob_psi_spol_stor=Booz_Jacob_psi_theta_zeta*(1+iota*Booz_dDzetastordspol)
        dXdspol=(1+iota*Booz_dDzetastordspol)*dXdtheta+Booz_dDzetastordspol*dXdzeta
        dYdspol=(1+iota*Booz_dDzetastordspol)*dYdtheta+Booz_dDzetastordspol*dYdzeta
        dZdspol=(1+iota*Booz_dDzetastordspol)*dZdtheta+Booz_dDzetastordspol*dZdzeta
        dXdstor=dXdzeta
        dYdstor=dYdzeta
        dZdstor=dZdzeta
        g_storstor = dXdstor**2 + dYdstor**2 + dZdstor**2
        g_spolspol = dXdspol**2 + dYdspol**2 + dZdspol**2
        g_spolstor = dXdstor*dXdspol + dYdstor*dYdspol + dZdstor*dZdspol
        Smith_dBdl=sqrt(dXdstor**2+dYdstor**2+dZdstor**2)
      #end if makeSmith
      
      ###########################################################
      # Nonlocal quantities are calculated if extended==True
      ###########################################################
      if extended:
        for ri in range(3):
          loc_h=1/(ext_B[ri]**2)
          loc_VPrimeHat=np.sum(loc_h)*4*np.pi**2/(Nzeta*Ntheta) #Note: VPrime=VPrimeHat*(G+iota*I)
          loc_FSAB2=4*np.pi**2/loc_VPrimeHat
          ext_Dzetaphi[ri]=mnFourierlib.mnmat(1-loc_h*loc_FSAB2,Nperiods=Nperiods).invJacBdotgrad(ext_iota[ri]).ifft()
                      
        self.ext_s=ext_s[ext_relrind]
        self.ext_psi=self.ext_s*Geom.psi_a
        self.ext_psi_a=Geom.psi_a
        #prepare for interpolating derivatives with resp to s evaluated at ext_s[ext_relrind]
        #f=c2*(s-s0)^2+c1*(s-s0)+c0
        #c2=d^-1*( f0*(s1-s2)+f1*(s2-s0)-f2*(s1-s0) )
        #c1=d^-1*( f0*(s2^2-s1^2)-f1*(s2-s0)^2+f2*(s1-s0)^2 )
        #d=(s1-s0)*(s2-s0)*(s1-s2)
        #df/ds(s)=2*c2*(s-s0)+c1
        ds10=ext_s[1]-ext_s[0]
        ds20=ext_s[2]-ext_s[0]
        ds12=ext_s[1]-ext_s[2]
        invdet=1/(ds10*ds20*ds12)
        dss0=self.ext_s-ext_s[0]
        dds_coef=invdet*np.array([ds12*(2*dss0-ds20-ds10),ds20*(2*dss0-ds20),ds10*(-2*dss0+ds10)])

        #Now perform the interpolations
        diotads=ext_iota[0]*dds_coef[0]+ext_iota[1]*dds_coef[1]+ext_iota[2]*dds_coef[2]
        dXds_appr=ext_X[0]*dds_coef[0]+ext_X[1]*dds_coef[1]+ext_X[2]*dds_coef[2]
        dYds_appr=ext_Y[0]*dds_coef[0]+ext_Y[1]*dds_coef[1]+ext_Y[2]*dds_coef[2]
        dZds_appr=ext_Z[0]*dds_coef[0]+ext_Z[1]*dds_coef[1]+ext_Z[2]*dds_coef[2]
        dDzetacylphids_appr=ext_Dzetacylphi[0]*dds_coef[0]+ext_Dzetacylphi[1]*dds_coef[1]+ext_Dzetacylphi[2]*dds_coef[2]
        dDzetaphids_appr=ext_Dzetaphi[0]*dds_coef[0]+ext_Dzetaphi[1]*dds_coef[1]+ext_Dzetaphi[2]*dds_coef[2]
        diotadpsi=diotads/self.ext_psi_a
        dDzetaphidpsiB_appr=dDzetaphids_appr/self.ext_psi_a
        # There is an uncertainty in the numerical s derivative. I treat it the following way:
        # First, assume that the above method gives the right direction of the vector e_psi but perhaps
        # not the right length. Fix the length by requiring that 
        # e_s dot e_theta x e_zeta = (G+iota*I)/B^2*psi_a
        # e_s_appr*corrfac dot e_theta x e_zeta = (G+iota*I)/B^2*psi_a
        # corrfac = (G+iota*I)/B^2*psi_a/(e_s_appr dot e_theta x e_zeta)
        tripprod_appr=( dXds_appr*(dYdtheta*dZdzeta-dYdzeta*dZdtheta)
                        -dYds_appr*(dXdtheta*dZdzeta-dXdzeta*dZdtheta)
                        +dZds_appr*(dXdtheta*dYdzeta-dXdzeta*dYdtheta))
        corrfac=Booz_Jacob_psi_theta_zeta*self.ext_psi_a/tripprod_appr
        dXdpsiB_appr2=dXds_appr*corrfac/self.ext_psi_a
        dYdpsiB_appr2=dYds_appr*corrfac/self.ext_psi_a
        dZdpsiB_appr2=dZds_appr*corrfac/self.ext_psi_a

        # e_psi has to have the form
        # e_psi=grad(psi)/gpsipsi + B_psi/B^2\vek{B} + kappa/B^2 \vek{B} x grad(psi)
        # with B_psi=B_psi00-(G+iota I)/<B^2> mu0 p' Dzetaphi
        # We first determine kappa and B_psi00=mean(B dot e_psi)=mean(J^-1(e_zeta+iota e_theta) dot e_psi)
        B_psi00=np.mean(((dXdzeta+iota*dXdtheta)*dXdpsiB_appr2+
                        (dYdzeta+iota*dYdtheta)*dYdpsiB_appr2+
                        (dZdzeta+iota*dZdtheta)*dZdpsiB_appr2)*B**2/(G+iota*I))
        B_psi=B_psi00-mu0dpdpsi*(G+iota*I)/FSAB2*Dzetaphi
        #and kappa=e_psi dot B x grad(psi)/gpsipsi=(G e_theta-I e_zeta) dot e_psi/(J gpsipsi)
        kappa= ((G*dXdtheta-I*dXdzeta)*dXdpsiB_appr2+
                (G*dYdtheta-I*dYdzeta)*dYdpsiB_appr2+
                (G*dZdtheta-I*dZdzeta)*dZdpsiB_appr2)/(G+iota*I)*B**2/gpsipsi
        g_psiBpsiB=1/gpsipsi+(B_psi**2+kappa*gpsipsi)/B**2
        g_psitheta=(B_psi*(g_thetazeta+iota*g_thetatheta)+kappa*(G*g_thetatheta-I*g_thetazeta))/(G+iota*I)
        g_zetapsi =(B_psi*(g_zetazeta +iota*g_thetazeta) +kappa*(G*g_thetazeta -I*g_zetazeta ))/(G+iota*I)

        #    g_psiBpsiB_appr2=dXdpsiB_appr2**2+dYdpsiB_appr2**2+dZdpsiB_appr2**2
        #    g_psitheta_appr2=dXdpsiB_appr2*dXdthetar+dYdpsiB_appr2*dYdtheta+dZdpsiB_appr2*dZdtheta
        #    g_zetapsi_appr2 =dXdpsiB_appr2*dXdzeta  +dYdpsiB_appr2*dYdzeta +dZdpsiB_appr2*dZdzeta

        #Other way to calculate the Jacobian
        #J^2=( g_psipsi*(g_thetatheta*g_zetazeta-g_thetazeta**2)
        #     -g_psitheta*(g_psitheta*g_zetazeta-g_psizeta*g_thetazeta)
        #     +g_psizeta*(g_psitheta*g_thetazeta-g_thetazeta*g_psizeta) )

        gthetatheta=(g_zetazeta*g_psiBpsiB-g_zetapsi**2)/Booz_Jacob_psi_theta_zeta**2
        gthetazeta =(g_zetapsi*g_psitheta-g_psiBpsiB*g_thetazeta)/Booz_Jacob_psi_theta_zeta**2
        gzetazeta  =(g_psiBpsiB*g_thetatheta-g_psitheta**2)/Booz_Jacob_psi_theta_zeta**2
        gpsitheta  =(g_thetazeta*g_zetapsi-g_psitheta*g_zetazeta)/Booz_Jacob_psi_theta_zeta**2
        gzetapsi   =(g_psitheta*g_thetazeta-g_zetapsi*g_thetatheta)/Booz_Jacob_psi_theta_zeta**2

        #Now turn to the Hamada metric tensor. See notes_BoozerToHamada.pdf for details
        #First the coordinate transformation matrix:
        dthetadpsiH=(iota*dDzetaphidpsiB_appr*B**2/FSAB2
              +Dzetaphi*diotadpsi*(1+iota*B**2/FSAB2*dDzetaphidtheta))
        dzetadpsiH=B**2/FSAB2*(dDzetaphidpsiB_appr+Dzetaphi*diotadpsi*dDzetaphidtheta)
        dthetadvthet=B**2/FSAB2*(1-dDzetaphidzeta)
        dzetadvthet =B**2/FSAB2*dDzetaphidtheta
        dthetadphi  =B**2/FSAB2*iota*dDzetaphidzeta
        dzetadphi   =B**2/FSAB2*(1-iota*dDzetaphidtheta)

        #and now the metric
        g_psiHpsiH=(g_psiBpsiB+dthetadpsiH**2*g_thetatheta+dzetadpsiH**2*g_zetazeta+
                    2*dthetadpsiH*g_psitheta+2*dzetadpsiH*g_zetapsi+2*dthetadpsiH*dzetadpsiH*g_thetazeta)
        g_psivthet=(dthetadvthet*g_psitheta+dzetadvthet*g_zetapsi+
                    dthetadpsiH*dthetadvthet*g_thetatheta+dzetadpsiH*dzetadvthet*g_zetazeta+
                    (dthetadpsiH*dzetadvthet+dzetadpsiH*dthetadvthet)*g_thetazeta)
        g_phipsi  =(dthetadphi*g_psitheta+dzetadphi*g_zetapsi+
                    dthetadpsiH*dthetadphi*g_thetatheta+dzetadpsiH*dzetadphi*g_zetazeta+
                    (dthetadpsiH*dzetadphi+dzetadpsiH*dthetadphi)*g_thetazeta)
        #dXdpsiH=dXdpsiB+dthetadpsiH_appr*dXdtheta+dzetadpsiH_appr*dXdzeta
        #dYdpsiH=dYdpsiB+dthetadpsiH_appr*dYdtheta+dzetadpsiH_appr*dYdzeta
        #dZdpsiH=dZdpsiB+dthetadpsiH_appr*dZdtheta+dzetadpsiH_appr*dZdzeta
        #g_psiHpsiH=dXdpsiH**2+dYdpsiH**2+dZdpsiH**2
        #g_psivthet=dXdpsiH*dXdvthet+dYdpsiH*dYdvthet+dZdpsiH*dZdvthet
        #g_phipsi  =dXdpsiH*dXdphi +dYdpsiH*dYdphi +dZdpsiH*dZdphi
        
        gvthetvthet=(g_phiphi*g_psiHpsiH-g_phipsi**2)/Jacob_psi_vthet_phi**2
        gvthetphi  =(g_phipsi*g_psivthet-g_psiHpsiH*g_vthetphi)/Jacob_psi_vthet_phi**2
        gphiphi    =(g_psiHpsiH*g_vthetvthet-g_psivthet**2)/Jacob_psi_vthet_phi**2
        gpsivthet  =(g_vthetphi*g_phipsi-g_psivthet*g_phiphi)/Jacob_psi_vthet_phi**2
        gphipsi    =(g_psivthet*g_vthetphi-g_phipsi*g_vthetvthet)/Jacob_psi_vthet_phi**2

        if np.any(np.abs(corrfac-1.0)>0.05):
          print('Probably too few surfaces in the equilibrium to get a good '+
                   'approximation for e_psi! Turning of extended!')
          extended=False
      #end if extended
      
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      #% Store variables in the self struct
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      self.FSAB2=FSAB2
      self.FSAu2B2=sum(sum(u**2*B**2*h))/sum(sum(h))
      self.FSAg_phiphi=sum(sum(g_phiphi*h))/sum(sum(h))
      self.FSAgpsipsi=sum(sum(gpsipsi*h))/sum(sum(h))
      self.FSAsqrtgpsipsi=sum(sum(np.sqrt(gpsipsi)*h))/sum(sum(h))
      self.Jacob_psi_vthet_phi=Jacob_psi_vthet_phi #This is a scalar
      self.dpsidrGraz=self.FSAsqrtgpsipsi*np.sign(psi_a)
      self.dsdrGraz  =self.FSAsqrtgpsipsi/abs(psi_a)
      self.r_eff=np.sqrt(self.FSAg_phiphi*self.FSAgpsipsi)/abs(self.G)
      self.r_eff_appr1=np.sqrt(sum(sum(gpsipsi*h/B**2))/sum(sum(h)))
      if (self.FSAg_phiphi-self.FSAu2B2-self.G**2/FSAB2)>=0:
          self.r_eff_appr2=np.sqrt((self.FSAg_phiphi-self.FSAu2B2-self.G**2/FSAB2)/iota**2)
      else:
          self.r_eff_appr2=0.0
      self.Ntor=Nzeta
      self.Npol=Ntheta
      self.Dtor=Dzeta
      self.Dpol=Dtheta
      self.StelSym=Geom.StelSym
      
      if name=='Boozer':
          self.Dzetaphi  =Booz_Dzetaphi
          self.Dzetapzeta=Booz_Dzetapzeta
          self.zeta  =zeta
          self.theta =theta
          self.phi   =Booz_phi
          self.vthet =Booz_vthet
          self.pzeta =Booz_pzeta       #pzeta=cylphi
          self.ptheta=Booz_ptheta
          
          self.R=R
          self.Z=Z
          self.cylphi=cylphi
          self.R00=np.mean(R)
          if not(Geom.StelSym):
            self.Z00=np.mean(Z)
          self.X=X
          self.Y=Y
          self.B=B
          self.u_chi=u
          self.u_psi=u/iota
          self.h=h
          self.Bpsitilde   = Booz_B_psi_tilde
          if extended:
              self.Bpsi00 = B_psi00
          self.dBpsidtheta = Booz_dBpsidtheta
          self.dBpsidzeta  = Booz_dBpsidzeta
          self.Jacob_psi_theta_zeta = Booz_Jacob_psi_theta_zeta
          self.Jacob_psi_ptheta_pzeta=Booz_Jacob_psi_ptheta_pzeta
          if makeSmith:
              self.Jacob_psi_spol_stor=Booz_Jacob_psi_spol_stor
          self.g_thetatheta = g_thetatheta
          self.g_thetazeta  = g_thetazeta
          self.g_zetazeta   = g_zetazeta
          if extended:
              self.g_psiBpsiB  = g_psiBpsiB #Mark with Booz because it is important which coords are held constant
              self.g_psitheta  = g_psitheta
              self.g_zetapsi   = g_zetapsi
              self.gthetatheta = gthetatheta
              self.gthetazeta  = gthetazeta
              self.gzetazeta   = gzetazeta
              self.gzetapsi    = gzetapsi
              self.gpsitheta   = gpsitheta              
          
          self.gpsipsi      = gpsipsi
          self.g_phiphi     = g_phiphi
          self.g_vthetvthet = g_vthetvthet
          self.g_vthetphi   = g_vthetphi
          if makeSmith:
              self.g_spolspol = g_spolspol
              self.g_storstor = g_storstor
              self.g_spolstor = g_spolstor
          if extended:
              self.g_psiHpsiH  = g_psiHpsiH
              self.g_psivthet  = g_psivthet
              self.g_phipsi    = g_phipsi
              self.gvthetvthet = gvthetvthet
              self.gvthetphi   = gvthetphi
              self.gphiphi     = gphiphi
              self.gphipsi     = gphipsi
              self.gpsivthet   = gpsivthet              
          
          self.gradpsidotgradB = Booz_gradpsidotgradB
          self.curv_normal  = Booz_curv_normal
          
          self.B00=B00
          self.Bmn=Bmn
          
          #if Geom.StelSym %sine components exist
          #            self.parity=parity
          #else
          #            self.parity=Geom.parity{rind}
          #end

          #self.gradpsidotgradBmn=mnFourierlib.fft(Booz_gradpsidotgradB)
          #self.gradpsidotgradBmn=mnmat.fft(Booz_gradpsidotgradB)
          #self.dBpsidthetamn=mnFourierlib.fft(Booz_dBpsidtheta)
          #self.dBpsidzetamn=mnmat.fft(Booz_dBpsidzeta)
          #self.hmn=mnFourierlib.fft(h)

      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      #% Interpolate to a phi vthet or a pzeta ptheta grid 
      #% Here, phi is the Hamada tor. coord. and vthet (\vartheta) is the Hamada poloidal coordinate.
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if name=='Hamada' or name=='Pest' or name=='Smith':
        #Choose the same resolution
        newtor=zeta  #This is just to copy the uniform grid
        newpol=theta #This is just to copy the uniform grid
      
        if name=='Hamada':
          self.phi=newtor
          self.vthet=newpol
          Booz_newtor=Booz_phi
          Booz_newpol=Booz_vthet
          BoozDzetanewtor=Booz_Dzetaphi
          #Make the interpolation
          New_Dzetanewtor=griddatacyclic(Booz_newpol,Booz_newtor,BoozDzetanewtor,self.Nperiods)
          self.zeta=newtor+New_Dzetanewtor
          self.theta=newpol+iota*New_Dzetanewtor
          self.Dzetaphi=New_Dzetanewtor
          self.Dzetapzeta=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_Dzetapzeta,
                                                        self.theta,self.zeta,self.Nperiods)
          if makeSmith:
              self.Dzetastor=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_Dzetastor,
                                                        self.theta,self.zeta,self.Nperiods)
          self.pzeta =self.zeta-self.Dzetapzeta
          self.ptheta=self.theta-iota*self.Dzetapzeta
        elif name=='Pest':
          self.pzeta=newtor
          self.ptheta=newpol
          Booz_newtor=Booz_pzeta
          Booz_newpol=Booz_ptheta
          BoozDzetanewtor=Booz_Dzetapzeta         
          #Make the interpolation
          New_Dzetanewtor=griddatacyclic(Booz_newpol,Booz_newtor,BoozDzetanewtor,self.Nperiods)
          self.zeta=newtor+New_Dzetanewtor
          self.theta=newpol+iota*New_Dzetanewtor
          self.Dzetaphi=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_Dzetaphi,
                                                      self.theta,self.zeta,self.Nperiods)
          self.Dzetapzeta=New_Dzetanewtor
          if makeSmith:
              self.Dzetastor=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_Dzetastor,
                                                      self.theta,self.zeta,self.Nperiods)
          self.phi  =self.zeta-self.Dzetaphi
          self.vthet=self.theta-iota*self.Dzetaphi
          
        elif name=='Smith':
          self.stor=newtor
          self.spol=newpol
          Booz_newtor=Booz_stor
          Booz_newpol=Booz_spol
          BoozDzetanewtor=BoozDzetastor         
          self.zeta=Smith_zeta
          self.theta=Smith_theta
          self.Dzetaphi=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_Dzetaphi,
                                                      self.theta,self.zeta,self.Nperiods)
          self.Dzetapzeta=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_Dzetapzeta,
                                                        self.theta,self.zeta,self.Nperiods)
          self.Dzetastor=SmithDzetastor
          self.phi  =self.zeta-self.Dzetaphi
          self.vthet=self.theta-iota*self.Dzetaphi
          self.pzeta =self.zeta-self.Dzetapzeta
          self.ptheta=self.theta-iota*self.Dzetapzeta
          
        #Now we come to the interpolation of a lot of stuff...
      

        self.R=fluxcoorddiscr.interp2_cyclic(theta,zeta,R,self.theta,self.zeta,self.Nperiods)
        self.Z=fluxcoorddiscr.interp2_cyclic(theta,zeta,Z,self.theta,self.zeta,self.Nperiods)
        self.cylphi=self.pzeta
        self.R00=np.mean(self.R)
        if not(Geom.StelSym):
          self.Z00=np.mean(self.Z)
        self.X=fluxcoorddiscr.interp2_cyclic(theta,zeta,X,self.theta,self.zeta,self.Nperiods)
        self.Y=fluxcoorddiscr.interp2_cyclic(theta,zeta,Y,self.theta,self.zeta,self.Nperiods)
        self.B=fluxcoorddiscr.interp2_cyclic(theta,zeta,B,self.theta,self.zeta,self.Nperiods)
        self.u_chi=fluxcoorddiscr.interp2_cyclic(theta,zeta,u,self.theta,self.zeta,self.Nperiods)
        self.u_psi=fluxcoorddiscr.interp2_cyclic(theta,zeta,u/iota,self.theta,self.zeta,self.Nperiods)
        self.h=fluxcoorddiscr.interp2_cyclic(theta,zeta,h,self.theta,self.zeta,self.Nperiods)
        self.Bpsitilde=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_B_psi_tilde,
                                                     self.theta,self.zeta,self.Nperiods)# Boozer specific
        self.dBpsidtheta=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_dBpsidtheta,
                                                       self.theta,self.zeta,self.Nperiods) #Boozer specific
        self.dBpsidzeta =fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_dBpsidzeta,
                                                       self.theta,self.zeta,self.Nperiods) #Boozer specific
        self.Jacob_psi_theta_zeta=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_Jacob_psi_theta_zeta,
                                                                self.theta,self.zeta,self.Nperiods) #Boozer specific
        self.Jacob_psi_ptheta_pzeta=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_Jacob_psi_ptheta_pzeta,
                                                                self.theta,self.zeta,self.Nperiods) #Boozer specific
        self.Jacob_psi_spol_stor=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_Jacob_psi_spol_stor,
                                                                self.theta,self.zeta,self.Nperiods) #Boozer specific
        self.g_thetatheta   =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_thetatheta,self.theta,self.zeta,self.Nperiods)
        self.g_thetazeta    =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_thetazeta,self.theta,self.zeta,self.Nperiods)
        self.g_zetazeta     =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_zetazeta,self.theta,self.zeta,self.Nperiods)
        if extended:
            self.g_psiBpsiB =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_psiBpsiB,self.theta,self.zeta,self.Nperiods)
            self.g_psitheta =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_psitheta,self.theta,self.zeta,self.Nperiods)
            self.g_zetapsi  =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_zetapsi,self.theta,self.zeta,self.Nperiods)
            self.gthetatheta=fluxcoorddiscr.interp2_cyclic(theta,zeta,gthetatheta,self.theta,self.zeta,self.Nperiods)
            self.gthetazeta =fluxcoorddiscr.interp2_cyclic(theta,zeta,gthetazeta,self.theta,self.zeta,self.Nperiods)
            self.gzetazeta  =fluxcoorddiscr.interp2_cyclic(theta,zeta,gzetazeta,self.theta,self.zeta,self.Nperiods)
            self.gzetapsi   =fluxcoorddiscr.interp2_cyclic(theta,zeta,gzetapsi,self.theta,self.zeta,self.Nperiods)
            self.gpsitheta  =fluxcoorddiscr.interp2_cyclic(theta,zeta,gpsitheta ,self.theta,self.zeta,self.Nperiods)
            
        self.gpsipsi        =fluxcoorddiscr.interp2_cyclic(theta,zeta,gpsipsi,self.theta,self.zeta,self.Nperiods)
        self.g_phiphi       =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_phiphi,self.theta,self.zeta,self.Nperiods)
        self.g_vthetvthet   =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_vthetvthet,self.theta,self.zeta,self.Nperiods)
        self.g_vthetphi     =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_vthetphi,self.theta,self.zeta,self.Nperiods)
        if extended:
            self.g_psiHpsiH =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_psiHpsiH,self.theta,self.zeta,self.Nperiods)
            self.g_psivthet =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_psivthet,self.theta,self.zeta,self.Nperiods)
            self.g_phipsi   =fluxcoorddiscr.interp2_cyclic(theta,zeta,g_phipsi,self.theta,self.zeta,self.Nperiods)
            self.gvthetvthet=fluxcoorddiscr.interp2_cyclic(theta,zeta,gvthetvthet,self.theta,self.zeta,self.Nperiods)
            self.gvthetphi  =fluxcoorddiscr.interp2_cyclic(theta,zeta,gvthetphi,self.theta,self.zeta,self.Nperiods)
            self.gphiphi    =fluxcoorddiscr.interp2_cyclic(theta,zeta,gphiphi,self.theta,self.zeta,self.Nperiods)
            self.gphipsi    =fluxcoorddiscr.interp2_cyclic(theta,zeta,gphipsi,self.theta,self.zeta,self.Nperiods)
            self.gpsivthet  =fluxcoorddiscr.interp2_cyclic(theta,zeta,gpsivthet ,self.theta,self.zeta,self.Nperiods)
        if makeSmith:
            self.g_spolspol = fluxcoorddiscr.interp2_cyclic(theta,zeta,g_spolspol,self.theta,self.zeta,self.Nperiods)
            self.g_storstor = fluxcoorddiscr.interp2_cyclic(theta,zeta,g_storstor,self.theta,self.zeta,self.Nperiods)
            self.g_spolstor = fluxcoorddiscr.interp2_cyclic(theta,zeta,g_spolstor,self.theta,self.zeta,self.Nperiods)
        
        self.gradpsidotgradB=fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_gradpsidotgradB,self.theta,self.zeta,self.Nperiods)
        self.curv_normal    =fluxcoorddiscr.interp2_cyclic(theta,zeta,Booz_curv_normal,self.theta,self.zeta,self.Nperiods)

        self.B00=np.mean(self.B)

        self.Bmn=mnFourierlib.mnmat(self.B,Nperiods=self.Nperiods)

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  def interpFromOtherCoords(self,quant_in_oldcoord,oldcoordname):
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      #u,v are th old coordinates
      Nu=quant_in_oldcoord.shape[0]
      Nv=quant_in_oldcoord.shape[1]
      Dv=2.0*np.pi/Nv/self.Nperiods
      Du=2.0*np.pi/Nu
      #make a uniform grid in u,v
      u, v = np.mgrid[0.0:2.0*np.pi-Du:1j*Nu,0.0:2.0*np.pi/self.Nperiods-Dv:1j*Nv]

      if self.name==oldcoordname:
          return quant_in_oldcoord
      
      if oldcoordname=='Pest': #u,v=ptheta,pzeta
          return fluxcoorddiscr.interp2_cyclic(u,v,quant_in_oldcoord,self.ptheta,self.pzeta,self.Nperiods)
      elif oldcoordname=='Hamada':  #u,v=vthet,phi
          return fluxcoorddiscr.interp2_cyclic(u,v,quant_in_oldcoord,self.vthet,self.phi,self.Nperiods)
      elif oldcoordname=='Boozer':  #u,v=theta,zeta
          return fluxcoorddiscr.interp2_cyclic(u,v,quant_in_oldcoord,self.theta,self.zeta,self.Nperiods)
      else:
          sys.exit('Unknown coordinates '+oldcoordname+' !')


                
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  def disp(self):
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        frm='{:8.4f}'
        print('----------------------------------------------------------------------------------')
        print('Discretisation on a uniform '+self.name+' coordinate grid of the magnetic field')
        print('\\mathbf{B} = I\\nabla\\theta + G\\nabla\\zeta + B_\\psi(\\theta,\\zeta)\\nabla\\psi')
        print('\\mathbf{B} = I\\nabla\\vthet + G\\nabla\\phi + \\nabla H(\\psi,\\vthet,\\phi)')
        print('----------------------------------------------------------------------------------')
        print('name     = '+self.name+' : the uniform coordinate grid')
        print('StelSym  = '+str(self.StelSym))
        print('Nperiods = '+str(self.Nperiods))
        print('Npol     = '+str(self.Npol))
        print('Ntor     = '+str(self.Ntor))
        print('Dpol     = '+frm.format(self.Dpol))
        print('Dtor     = '+frm.format(self.Dtor))
        print('iota     = '+frm.format(self.iota))
        print('G        = '+frm.format(self.G))
        print('I        = '+frm.format(self.I))
        print('mu0dpdpsi= '+frm.format(self.mu0dpdpsi))
        print('FSAB2    = '+frm.format(self.FSAB2)+  '            : <B^2>')
        print('FSAu2B2  = '+frm.format(self.FSAu2B2)+'            : <u^2B^2>')
        print('FSAg_phiphi    = '+frm.format(self.FSAg_phiphi)+'      : <g_phiphi>')
        print('FSAgpsipsi     = '+frm.format(self.FSAgpsipsi)+ '      : <g^psipsi>')
        print('FSAsqrtgpsipsi = '+frm.format(self.FSAsqrtgpsipsi)+'      : <sqrt(g^psipsi)>')
        print('Jacob_psi_vthet_phi = '+frm.format(self.Jacob_psi_vthet_phi)+' : Hamada Jacobian (G+iota*I)/<B^2>')
        print('dpsidrGraz     = '+frm.format(self.dpsidrGraz)+'      : dpsi/dr_Graz')
        print('dsdrGraz       = '+frm.format(self.dsdrGraz)+'      : ds/dr_Graz')
        print('r_eff          = '+frm.format(self.r_eff)+'      : r_eff = sqrt(<g_phiphi><g^psipsi>)/|G|')
        print('R00 = '+frm.format(self.R00)+'                 : m=n=0 Fourier coeff. in '+self.name+' coordinates')
        if not(self.StelSym):
            print('Z00 = '+frm.format(self.Z00)+'                 : m=n=0 Fourier coeff. in '+self.name+' coordinates')
        print('B00 = '+frm.format(self.B00)+'                 : m=n=0 Fourier coeff. in '+self.name+' coordinates')
        print('------------------------------------------------------------------------------')
        print(' Npol x Ntor Scalar fields')
        print('-------------------------------------')
        print('zeta,theta     : Boozer coordinates')
        print('phi,vthet      : Hamada coordinates')
        print('pzeta,ptheta   : Pest coordinates (pzeta=cylphi)')
        print('Dzetaphi       : zeta-phi')
        print('Dzetapzeta     : zeta-pzeta')
        print('R,Z,cylphi     : Cylindrical coordinates (right handed in this order)')
        print('X,Y,Z          : Cartesian coordinates')
        print('B, h           : B and h=1/B^2')
        print('u_chi          : defined by <uB^2>=0 and (B dot grad) u = 2 iota B^-3 B x nabla psi dot nabla B')
        print('u_psi          : defined by <uB^2>=0 and (B dot grad) u = 2 B^-3 B x nabla psi dot nabla B')
        print('Bpsitilde      : B_psi - B_psi_00  (Boozer specific)')
        print('dBpsidtheta    : dB_psi/dtheta (Boozer specific)')
        print('dBpsidzeta     : dB_psi/dzeta (Boozer specific)')
        print('Jacob_psi_theta_zeta                  : Boozer Jacobian (G+iota*I)/B^2')
        print('g_thetatheta, g_thetazeta, g_zetazeta : Boozer metric')
        print('g_vthetvthet, g_vthetphi,  g_phiphi   : Hamada metric')
        print('gpsipsi         : |grad(psi)|^2')
        print('gradpsidotgradB : grad(psi) dot grad(B)')
        print('curv_normal     : normal curvature grad(psi)/|grad(psi)| dot (b dot grad)b')
        print('------------------------------------------------------------------------------')
        print(' mnmat objects (Fourier transformed Npol x Ntor Scalar fields)')
        print('---------------------------------------------------------------')
        print('Bmn')
        print('------------------------------------------------------------------------------')

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  def plot(self,toshow,title='',cmap=None):
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if cmap==None:
           #cmap='coolwarm'
           #cmap='RdBu'
           cmap='jet'

        if isinstance(toshow, str):
            fun=getattr(self,toshow)
        else:
            fun=toshow

        if self.name=='Boozer':
            x=self.zeta
            y=self.theta
        elif self.name=='Hamada':
            x=self.phi
            y=self.vthet
        elif self.name=='Pest':
            x=self.pzeta
            y=self.ptheta
        elif self.name=='Smith':
            x=self.stor
            y=self.spol
            
        z = fun[:-1, :-1]
        z_min, z_max = z.min(), z.max()

        fig = plt.figure()
        plt.subplot()
        plt.pcolor(x, y, z, cmap=cmap, vmin=z_min, vmax=z_max)

        plt.axis([x.min(), x.max(), y.min(), y.max()])
        ax = fig.gca()
        ax.set_xlabel(self.name+' toroidal coordinate')
        ax.set_ylabel(self.name+' poloidal coordinate')
        if title!='':
            ax.set_title(title)
        else:
            ax.set_title(toshow)
        plt.colorbar()
        #plt.show()
        return fig, ax
        
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  def plot3d(self,toshow,title='',torstride=1,polstride=1,cmap=None):
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        #help routine to plot routine
        def makefulltordata(input,N=0,type='periodic'):
          if isinstance(input,mnFourierlib.mnmat):
            if N!=0 and input.Nperiods!=N:
                sys.exit('Incompatible number of periods')
            inp=input.ifft()
            N=np.int(input.Nperiods)
          elif isinstance(input,np.ndarray):
            if N==0:
              sys.exit('Please provide number of field periods N!')
            inp=input
          else:
            sys.exit('Non-recognised input!')

          N=np.int(N)
            
          torPeriod=0
          polPeriod=0
          if type=='toroidal coordinate':
            torPeriod=2*np.pi/N
          elif type=='poloidal coordinate':
            polPeriod=2*np.pi

          out=np.zeros((inp.shape[0]+1,N*inp.shape[1]+1))
          for ind in range(N):
              out[0:inp.shape[0], ind*inp.shape[1]:(ind+1)*inp.shape[1]] = inp+ind*torPeriod

          out[0:inp.shape[0],N*inp.shape[1]]=inp[:,0]+N*torPeriod
          out[inp.shape[0],:]=out[0,:]+polPeriod

          return out
          
        #main plot routine
        if isinstance(toshow, str):
            fun=getattr(self,toshow)
        else:
            fun=toshow

        geomangf=-makefulltordata(self.cylphi,self.Nperiods,'toroidal coordinate')
        Rf=makefulltordata(self.R,self.Nperiods,'periodic')
        Zf=makefulltordata(self.Z,self.Nperiods,'periodic')
        Xf=Rf*np.cos(geomangf)
        Yf=Rf*np.sin(geomangf)
        if toshow=='zeta' or toshow=='phi'or toshow=='pzeta':
            funf=makefulltordata(fun,self.Nperiods,'toroidal coordinate')
        elif toshow=='theta' or toshow=='vthet'or toshow=='ptheta':
            funf=makefulltordata(fun,self.Nperiods,'poloidal coordinate')
        else:
            funf=makefulltordata(fun,self.Nperiods,'periodic') #assume periodic, change this

        #print Xf.shape
        #print Yf.shape
        #print Zf.shape
        #print funf.shape
            
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        #ax = fig.add_subplot(111, projection='3d')
        my_col=cm.jet(funf/np.amax(funf))

        color_dimension=funf
        minn, maxx = color_dimension.min(), color_dimension.max()
        norm = matplotlib.colors.Normalize(minn, maxx)
        if cmap==None:
           #m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
           m = plt.cm.ScalarMappable(norm=norm, cmap='coolwarm')
        else:
           m = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
           
        m.set_array([])
        fcolors = m.to_rgba(color_dimension)
        m.set_array(funf)

        funfnorm=(funf-minn)/(maxx-minn)
        
        #print fcolors.shape
        #surf=ax.plot_surface(Xf, Yf, Zf,linewidth=0,cmap=cm.coolwarm,rstride=1,cstride=1, color=funfnorm, vmin=minn, vmax=maxx, shade=False)
        surf=ax.plot_surface(Xf, Yf, Zf,rstride=polstride, cstride=torstride, facecolors=fcolors, vmin=minn, vmax=maxx, shade=False)#linewidth=0,rstride=1,cstride=1, facecolor=fcolors, vmin=minn, vmax=maxx, shade=False)
        plt.colorbar(m)#fig.colorbar(surf)#, shrink=0.5, aspect=5)
        #fig.colorbar(surf, shrink=0.5, aspect=5)
        #surf=ax.plot_wireframe(Xf, Yf, Zf,rstride=polstride,cstride=torstride)#,linewidth=0,cmap=cm.coolwarm,)#, color=funf)


        # Create cubic bounding box to simulate equal aspect ratio
        max_range = np.array([Xf.max()-Xf.min(), Yf.max()-Yf.min(), Zf.max()-Zf.min()]).max()
        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(Xf.max()+Xf.min())
        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Yf.max()+Yf.min())
        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Zf.max()+Zf.min())
        # Comment or uncomment following both lines to test the fake bounding box:
        for xb, yb, zb in zip(Xb, Yb, Zb):
            ax.plot([xb], [yb], [zb], 'w') 

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        if title!='':
            ax.set_title(title)
            
        #plt.show()
        return fig, ax



    
