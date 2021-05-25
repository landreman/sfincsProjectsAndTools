#!/usr/bin/env python
from __future__ import division
import numpy as np
import sys, time, multiprocessing
from scipy import interpolate
import matplotlib.pyplot as plt
import geomlib
import mnFourierlib
from fluxcoorddiscr import fluxcoorddiscr
##########################################################################
# Function crossectplot
# Plots an input quantity (toshow) on poloidal cross sections with
# selected toroidal Pest coordinates (Pest_zetas).
#
# Inputs
# ---------------
# toshow    : The quantity to plot. A list (of the same length as the input
#             rNs) where each list element is a numpy 2d array with values
#             on the respective flux surfaces on a uniform grid in the
#             poloidal and toroidal coordinates specified in the input
#             coordname.
# rNs       : Normalised radii of the flux surfaces in
# Geom      : A bcgeom object containing the geometry
# coordname : The types of poloidal and toroidal coordinates on which
#             toshow is given. Can be 'Boozer', 'Pest' or 'Hamada'.
# Pest_zetas: The toroidal coordinates specifying the cross sections to
#             plot. The number of plots produced is len(Pest_zetas).
# rNmin     : Nothing is shown inside this radius 
# rNmax     : Nothing is shown outside this radius
# rNcontours: A list of normalised radii where flux surface contours are
#             outlined
# savefilename : The base of the filename for the figure files. The
#                toroidal angle Pest_zeta (in degrees) will be added to
#                the filename for uniqueness.
#
# Output
# ---------------
# figures   : A list with the handles to the produced figures
#
##########################################################################
def oneradius(Geom,ridx,Npol,Ntor,myquant,coordname,pzweigh):
    #This is used for parallelizing the radius loop  
    Pest=fluxcoorddiscr(Geom,ridx,Npol,Ntor,name='Pest')
    QuantPest_fullsurf=Pest.interpFromOtherCoords(myquant,coordname)
    quantPest_ri=np.tensordot(pzweigh,QuantPest_fullsurf,axes=(1,1))
    R_ri=np.tensordot(pzweigh,Pest.R,axes=(1,1))
    Z_ri=np.tensordot(pzweigh,Pest.Z,axes=(1,1))
    return (quantPest_ri,R_ri,Z_ri)

def crossectplot(toshow,rNs,Geom,coordname='Boozer',Pest_zetas=0.0,
                 rNmin=0.0,rNmax=np.inf,rNcontours=None,savefilename=None):

  if not(isinstance(Geom,geomlib.bcgeom)):
    sys.exit('crossectplot has only been implemented for bcgeom input so far.')
  Nperiods=Geom.Nperiods
  if len(rNs)!=len(toshow):
      sys.exit('len(rNs) must be equal to len(toshow)')
  if np.any(np.diff(rNs)<=0):
      sys.exit('rNs must be an ascending array!')
  if rNmax<=rNmin:
      sys.exit('rNmax must be > rNmin!')
  if np.any(np.isnan(Pest_zetas)):
      sys.exit('Error: NaNs in the array Pest_zetas!')
      
  #Interpolation
  rNmin=max(rNmin,rNs[0])
  rNmax=min(rNmax,rNs[-1])
  rNmin=max(rNmin,Geom.rnorm[0])
  rNmax=min(rNmax,Geom.rnorm[-1])
  rnormminind=np.where(Geom.rnorm>=rNmin)[0][0]
  rnormmaxind=np.where(Geom.rnorm<=rNmax)[0][-1]
  if rnormminind.size==0:
      sys.exit('Decrease rNmin! No surfaces above it found in Geom.')
  if rnormmaxind.size==0:
      sys.exit('Increase rNmax! No surfaces below it found in Geom.')
  if np.isscalar(Pest_zetas):
      Pest_zetas=np.array([Pest_zetas])
      
  Npol=toshow[0].shape[0]
  Ntor=toshow[0].shape[1]
  for rind in range(1,len(toshow)):
      Npol=max(toshow[rind].shape[0],Npol)
      Ntor=max(toshow[rind].shape[1],Ntor)

  #Now expand all toshow[i]'s to this resolution
  tmpquant=np.zeros((len(toshow),Npol,Ntor)) #make an expanded copy
  for rind in range(len(toshow)):
      if toshow[rind].shape[0] == Npol and toshow[rind].shape[1] == Ntor:
          tmpquant[rind]=toshow[rind][:,:]
      else:
          tmpquant[rind]=mnFourierlib.mnmat(toshow[rind],Npol,Ntor,Nperiods).ifft()

  ridxs=np.zeros((len(toshow)),dtype=int)
  for rind in range(len(toshow)):
      ridxs[rind] = (np.abs(rNs[rind]-Geom.rnorm)).argmin()

  #check for myrNs outside rNmin,rNmax
  ridxs[np.where(ridxs<rnormminind)[0]]=rnormminind
  ridxs[np.where(ridxs>rnormmaxind)[0]]=rnormmaxind

  myrNs = Geom.rnorm[ridxs]

  #check for doubles
  mask = np.ones(len(myrNs), dtype=bool)
  mask[np.where(np.diff(myrNs)==0)[0]]=False
  myrNs=myrNs[mask]
  ridxs=ridxs[mask]

  interpfun=interpolate.interp1d(rNs,tmpquant,axis=0,kind='quadratic')

  myquant=interpfun(myrNs)

  uniformzcoord=np.linspace(0.0,2.0*np.pi/Nperiods,Ntor,endpoint=False)
  Dz=uniformzcoord[1]-uniformzcoord[0]

  pzweigh=np.zeros((len(Pest_zetas),Ntor))
  quantPest=np.zeros((len(Pest_zetas),len(myrNs),Npol))
  R=np.zeros((len(Pest_zetas),len(myrNs),Npol))
  Z=np.zeros((len(Pest_zetas),len(myrNs),Npol))

  for pzind in range(len(Pest_zetas)):
      Lind=np.where(uniformzcoord<=np.mod(Pest_zetas[pzind],2.0*np.pi/Nperiods))[0][-1]
      x=np.mod(Pest_zetas[pzind],2.0*np.pi/Nperiods)-uniformzcoord[Lind]
      pzweigh[pzind,Lind]=(Dz-x)/Dz
      pzweigh[pzind,np.mod(Lind+1,Ntor)]=x/Dz

  t0=time.time()
  parallelize=True
  if parallelize:
    Ncpu = multiprocessing.cpu_count()
    args=[(Geom,ridxs[ri],Npol,Ntor,
           myquant[ri],coordname,pzweigh) for ri in range(len(myrNs))]
    with multiprocessing.Pool(processes=Ncpu) as pool:
      out=pool.starmap(oneradius,args) #The distributed calculation
    for ri in range(len(myrNs)):
      quantPest[:,ri,:]=out[ri][0]  
      R[:,ri,:]=out[ri][1]  
      Z[:,ri,:]=out[ri][2]  
  else:       
    for ri in range(len(myrNs)):
      Pest=fluxcoorddiscr(Geom,ridxs[ri],Npol,Ntor,name='Pest')
      QuantPest_fullsurf=Pest.interpFromOtherCoords(myquant[ri],coordname)
      quantPest[:,ri,:]=np.tensordot(pzweigh,QuantPest_fullsurf,axes=(1,1))
      R[:,ri,:]=np.tensordot(pzweigh,Pest.R,axes=(1,1))
      Z[:,ri,:]=np.tensordot(pzweigh,Pest.Z,axes=(1,1))
      
  print('loop time in crossectplot: '+str(time.time()-t0)+' s')
              
  quantPest=np.append(quantPest,quantPest[:,:,[0]],axis=2)
  R=np.append(R,R[:,:,[0]],axis=2)
  Z=np.append(Z,Z[:,:,[0]],axis=2)

  #Now prevent skewness by using a half grid in theta
  Rext=(R[:,:,:-1]+R[:,:,1:])/2.0
  Zext=(Z[:,:,:-1]+Z[:,:,1:])/2.0
  Rext=np.append(Rext,R[:,:,[0]],axis=2)
  Rext=np.append(R[:,:,[0]],Rext,axis=2)
  Zext=np.append(Zext,Z[:,:,[0]],axis=2)
  Zext=np.append(Z[:,:,[0]],Zext,axis=2)
  #This last value is not shown in the plot:
  quantPestext=np.append(quantPest,quantPest[:,:,[0]],axis=2)

  rNM=np.outer(myrNs,np.ones((R.shape[2])))

  figs=[None]*len(Pest_zetas)
  axs=[None]*len(Pest_zetas)
  
  for pzind in range(len(Pest_zetas)):
    (figs[pzind],axs[pzind])=plt.subplots(1,1)
    #print('hhh') 
    #figures[pzind]=plt.figure()
    if not(rNcontours is None):
      if len(rNcontours)>0:
        CS=plt.contour(R[pzind],Z[pzind],rNM,colors='k',levels=rNcontours)
        plt.clabel(CS, CS.levels, inline=True, fmt='%.1f', fontsize=10)
    plt.pcolor(Rext[pzind],Zext[pzind],quantPestext[pzind])
    plt.axis('equal')
    plt.colorbar()

    if not(savefilename is None):
        degrees=str(int(round(Pest_zetas[pzind]*180/np.pi)))
        if len(savefilename)<5:
            filename=savefilename+degrees+'.eps'
        elif savefilename[-4:]=='.eps':
            filename=savefilename[:-4]+degrees+'.eps'
        elif savefilename[-4:]=='.pdfs':
            filename=savefilename[:-4]+degrees+'.pdf'
        else:
            filename=savefilename+degrees+'.eps'
        plt.savefig(filename)

  return figs, axs
    
    
