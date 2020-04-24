#!/usr/bin/env python
from __future__ import division
import sys, os, multiprocessing
import numpy as np

from mnFourierlib import mnmat
from mnFourierlib import mnlist
from geomlib import bcgeom
from geomlib import vmecgeom
from fluxcoorddiscr import fluxcoorddiscr
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from crossectplot import crossectplot

def wait_for_user(message=None):
    if message is None:
        message = 'Press the <ENTER> key to continue...'
    input(message)
    
def findpenetrated(Xf,Yf,Zf,P0,P1,fig=None,ax=None):
  # Finds surface elements penetrated by the line P0->P1.
  # Found elements are not allowed to be further away from P0 than P1,
  # allowing the user to exclude points on the other side of the torus,
  # when the line is approximately horisontal  
  p=np.array([P1[0]-P0[0],P1[1]-P0[1],P1[2]-P0[2]])
  absp=np.sqrt(np.sum(p**2))
  phat=p/absp
  q0f=np.zeros((3,Xf.shape[0],Xf.shape[1]))
  q0f[0]=Xf-P0[0]
  q0f[1]=Yf-P0[1]
  q0f[2]=Zf-P0[2]
  q0=q0f[:,:-1,:-1]
  absq0=np.sqrt(np.sum(q0**2,axis=0))
  
  pcrossq0f=np.zeros_like(q0f)
  pcrossq0f[0]=p[1]*q0f[2]-p[2]*q0f[1]
  pcrossq0f[1]=p[2]*q0f[0]-p[0]*q0f[2]
  pcrossq0f[2]=p[0]*q0f[1]-p[1]*q0f[0]
  pcrossq0=pcrossq0f[:,:-1,:-1]

  q1=np.zeros_like(q0)
  q2=np.zeros_like(q0)
  q3=np.zeros_like(q0)
  pcrossq1=np.zeros_like(q0)
  pcrossq2=np.zeros_like(q0)
  pcrossq3=np.zeros_like(q0)
  q1=q0f[:,:-1,1:]
  q3=q0f[:,1:,:-1]
  q2=q0f[:,1:,1:]
  pcrossq1=pcrossq0f[:,:-1,1:]
  pcrossq3=pcrossq0f[:,1:,:-1]
  pcrossq2=pcrossq0f[:,1:,1:]

  tripl01=np.sum(pcrossq0*q1,axis=0)
  tripl12=np.sum(pcrossq1*q2,axis=0)
  tripl23=np.sum(pcrossq2*q3,axis=0)
  tripl30=np.sum(pcrossq3*q0,axis=0)

  if False:
    thi=3
    zi=4
    if tripl01.shape[0]>3:
      print('q0='+str(q0[:,thi,zi]))
      print('q1='+str(q1[:,thi,zi]))
      print('q2='+str(q2[:,thi,zi]))
      print('q3='+str(q3[:,thi,zi]))
      #print('tripl01='+str(tripl01[thi,zi]))
      #print('tripl12='+str(tripl12[thi,zi]))
      #print('tripl23='+str(tripl23[thi,zi]))
      #print('tripl30='+str(tripl30[thi,zi]))
    else:
      if fig is not None:
        i=1
        j=1
        ax[0].plot(q0[0,i,j],q0[1,i,j],'d')
        fig.show()
        wait_for_user()
      print('q0='+str(q0))
      print('q1='+str(q1))
      print('q2='+str(q2))
      print('q3='+str(q3))
      #print('tripl01='+str(tripl01))
      #print('tripl12='+str(tripl12))
      #print('tripl23='+str(tripl23))
      #print('tripl30='+str(tripl30))

  goodID=np.where(np.logical_and(
                  np.logical_and(absq0<absp,tripl01*tripl12>0),
                  np.logical_and(tripl12*tripl23>0,tripl23*tripl30>0)))
  
  signpdotgradr=np.sign(tripl01[goodID])
  possign=np.where(tripl01[goodID]>=0)[0]
  negsign=np.where(tripl01[goodID]<0)[0]
  if np.any(absq0[goodID]>=absp):
    print('absq0[goodID]<absp='+str(absq0[goodID]<absp))
    print('absq0 rule used for exclusion!')

  signorderedID=[None]*2
  signorderedID[0]=(goodID[0][possign],goodID[1][possign])
  signorderedID[1]=(goodID[0][negsign],goodID[1][negsign])
  
  return goodID,signpdotgradr,signorderedID

def refine(P0,P1,corn,thetacorn,zetacorn,Rmn,Zmn,Dzetacylphimn,rind,ii,irefine,
           fig=None,ax=None):
  #Refine the position of the line crossing with the surface element with the corners corn   
  out_corn     =np.nan*np.zeros_like(corn)      #X,Y,Z values at the refined corners
  out_thetacorn=np.nan*np.zeros_like(thetacorn) #theta values at the refined corners
  out_zetacorn =np.nan*np.zeros_like(zetacorn)  #zeta values at the the refined corners
  if np.any(np.isnan(corn)):
    return out_corn,out_thetacorn,out_zetacorn
  #original corner numbering: [0 1]
  #                           [3 2]
  #numbering of double grid: [0   01   1 ]
  #                          [30  Mid  12]
  #                          [3   23   2 ] 
  theta01=(thetacorn[0]+thetacorn[1])/2.0
  theta12=(thetacorn[1]+thetacorn[2])/2.0
  theta23=(thetacorn[2]+thetacorn[3])/2.0
  theta30=(thetacorn[3]+thetacorn[0])/2.0
  thetaMid=(theta01+theta23)/2.0    
  zeta01=(zetacorn[0]+zetacorn[1])/2.0
  zeta12=(zetacorn[1]+zetacorn[2])/2.0
  zeta23=(zetacorn[2]+zetacorn[3])/2.0
  zeta30=(zetacorn[3]+zetacorn[0])/2.0
  zetaMid=(zeta01+zeta23)/2.0
  R01=Rmn.evalpoint(theta01,zeta01)    
  R12=Rmn.evalpoint(theta12,zeta12)    
  R23=Rmn.evalpoint(theta23,zeta23)    
  R30=Rmn.evalpoint(theta30,zeta30)    
  RMid=Rmn.evalpoint(thetaMid,zetaMid)    
  Z01=Zmn.evalpoint(theta01,zeta01)    
  Z12=Zmn.evalpoint(theta12,zeta12)    
  Z23=Zmn.evalpoint(theta23,zeta23)    
  Z30=Zmn.evalpoint(theta30,zeta30)    
  ZMid=Zmn.evalpoint(thetaMid,zetaMid)    
  cylphi01=zeta01-Dzetacylphimn.evalpoint(theta01,zeta01)    
  cylphi12=zeta12-Dzetacylphimn.evalpoint(theta12,zeta12)    
  cylphi23=zeta23-Dzetacylphimn.evalpoint(theta23,zeta23)    
  cylphi30=zeta30-Dzetacylphimn.evalpoint(theta30,zeta30)    
  cylphiMid=zetaMid-Dzetacylphimn.evalpoint(thetaMid,zetaMid)    
  X01=R01*np.cos(-cylphi01)
  X12=R12*np.cos(-cylphi12)
  X23=R23*np.cos(-cylphi23)
  X30=R30*np.cos(-cylphi30)
  XMid=RMid*np.cos(-cylphiMid)
  Y01=R01*np.sin(-cylphi01)
  Y12=R12*np.sin(-cylphi12)
  Y23=R23*np.sin(-cylphi23)
  Y30=R30*np.sin(-cylphi30)
  YMid=RMid*np.sin(-cylphiMid)
  
  Xf=np.array([[corn[0,0],X01,  corn[0,1]],
              [X30,    XMid, X12    ],
              [corn[0,3],X23,  corn[0,2]]])
  Yf=np.array([[corn[1,0],Y01,  corn[1,1]],
              [Y30,      YMid, Y12    ],
              [corn[1,3],Y23,  corn[1,2]]])
  Zf=np.array([[corn[2,0],Z01,  corn[2,1]],
              [Z30,      ZMid, Z12    ],
              [corn[2,3],Z23,  corn[2,2]]])
  XYZf=np.array([Xf,Yf,Zf])
  #print('XYZf=')
  #print(XYZf)
  #print(XYZf.shape)
  thetaf=np.array([[thetacorn[0],  theta01,  thetacorn[1]],
                   [theta30,       thetaMid, theta12    ],
                   [thetacorn[3],  theta23,  thetacorn[2]]])
  zetaf=np.array([[zetacorn[0],  zeta01,  zetacorn[1]],
                  [zeta30,       zetaMid, zeta12    ],
                  [zetacorn[3],  zeta23,  zetacorn[2]]])
  
  goodID,signpdotgradr,signorderedID=findpenetrated(Xf,Yf,Zf,P0,P1,fig=fig,ax=ax)
  #print(goodID)
  if goodID[0].size==0:
    print('found no crossing! rind='+str(rind)+', crossing='+str(ii)+'irefine='+str(irefine)+'Try increasing Ntheta and Nzeta!')
    return out_corn,out_thetacorn,out_zetacorn
  if len(goodID[0])>1:
    sys.exit('strange!')
    return out_corn,out_thetacorn,out_zetacorn
  thind=goodID[0][0]
  zind=goodID[1][0]
  #print(thind)
  #print(zind)
  out_corn[:,0]=XYZf[:,thind,zind]
  out_corn[:,1]=XYZf[:,thind+1,zind]
  out_corn[:,2]=XYZf[:,thind+1,zind+1]
  out_corn[:,3]=XYZf[:,thind,zind+1]
  out_thetacorn[0]=thetaf[thind,zind]
  out_thetacorn[1]=thetaf[thind+1,zind]
  out_thetacorn[2]=thetaf[thind+1,zind+1]
  out_thetacorn[3]=thetaf[thind,zind+1]
  out_zetacorn[0]=zetaf[thind,zind]
  out_zetacorn[1]=zetaf[thind+1,zind]
  out_zetacorn[2]=zetaf[thind+1,zind+1]
  out_zetacorn[3]=zetaf[thind,zind+1]
  return out_corn,out_thetacorn,out_zetacorn

def makefulltordata(input,N=0,type='periodic'):
  if isinstance(input,mnmat):
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

def onesurf_crossings(rind,Geom,P0,P1,Ntheta,Nzeta,Nrefinements=10,rind_plot=-1):
  p=np.array([P1[0]-P0[0],P1[1]-P0[1],P1[2]-P0[2]])
  #print('r/a='+str(Geom.rnorm[rind]))

  #make the grid
  Dtheta=2*np.pi/Ntheta
  Dzeta=2*np.pi/Nzeta/Geom.Nperiods
  theta, zeta = np.mgrid[0.0:2.0*np.pi-Dtheta:1j*Ntheta,0.0:2.0*np.pi/Geom.Nperiods-Dzeta:1j*Nzeta]

  Rmn=mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='R')
  Zmn=mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='Z')
  Dzetacylphimn=mnmat(Geom,Ntheta=Ntheta,Nzeta=Nzeta,rind=rind,quantity='Dphi')
  Dzetacylphi=Dzetacylphimn.ifft()
  cylphi=zeta-Dzetacylphi
  geomang=-cylphi #(R,Z,cylphi) and (R,geomang,Z) are right handed systems. 
  R=Rmn.ifft()
  Z=Zmn.ifft()
  X=R*np.cos(geomang)
  Y=R*np.sin(geomang)
  zetaf=makefulltordata(zeta,Geom.Nperiods,'toroidal coordinate')
  thetaf=makefulltordata(theta,Geom.Nperiods,'poloidal coordinate')
  geomangf=-makefulltordata(cylphi,Geom.Nperiods,'toroidal coordinate')
  Rf=makefulltordata(R,Geom.Nperiods,'periodic')
  Zf=makefulltordata(Z,Geom.Nperiods,'periodic')
  Xf=Rf*np.cos(geomangf)
  Yf=Rf*np.sin(geomangf)

  goodID,signpdotgradr,signorderedID=findpenetrated(Xf,Yf,Zf,P0,P1)

  #print('goodID='+str(goodID))
  #print('goodID[0]='+str(goodID[0]))
  #sys.exit()

  zetacorn =np.nan*np.zeros((len(goodID[0]),4))
  thetacorn=np.nan*np.zeros((len(goodID[0]),4))
  corn     =np.nan*np.zeros((len(goodID[0]),3,4))


  colors=['b','g','c','m','y','k','r']

  Booz = fluxcoorddiscr(Geom,rind,Ntheta,Nzeta,name='Boozer')

  if rind==rind_plot:
    fig1,ax1=Booz.plot3d('B',title='Magnetic field',torstride=1,polstride=1,cmap=None,wireframe=True)
    ax1.plot([P0[0],P1[0]],[P0[1],P1[1]],[P0[2],P1[2]],'r')

    xl=ax1.get_xlim()
    yl=ax1.get_ylim()
    zl=ax1.get_zlim()
    fig1.show()

    fig2,ax2=plt.subplots(1,1,sharex=True,figsize=(6,7),subplot_kw={'projection': '3d'})
    #fig3,ax3=plt.subplots(3,1,sharex=True,figsize=(6,7))

  for ii in range(len(goodID[0])):
    ipol=goodID[0][ii]
    itor=goodID[1][ii]
    zetacorn[ii,0]=zetaf[ipol,itor]
    thetacorn[ii,0]=thetaf[ipol,itor]
    zetacorn[ii,1]=zetaf[ipol+1,itor]
    thetacorn[ii,1]=thetaf[ipol+1,itor]
    zetacorn[ii,2]=zetaf[ipol+1,itor+1]
    thetacorn[ii,2]=thetaf[ipol+1,itor+1]
    zetacorn[ii,3]=zetaf[ipol,itor+1]
    thetacorn[ii,3]=thetaf[ipol,itor+1]
    corn[ii,0,0]=Xf[ipol,itor]
    corn[ii,1,0]=Yf[ipol,itor]
    corn[ii,2,0]=Zf[ipol,itor]
    corn[ii,0,1]=Xf[ipol+1,itor]
    corn[ii,1,1]=Yf[ipol+1,itor]
    corn[ii,2,1]=Zf[ipol+1,itor]
    corn[ii,0,2]=Xf[ipol+1,itor+1]
    corn[ii,1,2]=Yf[ipol+1,itor+1]
    corn[ii,2,2]=Zf[ipol+1,itor+1]
    corn[ii,0,3]=Xf[ipol,itor+1]
    corn[ii,1,3]=Yf[ipol,itor+1]
    corn[ii,2,3]=Zf[ipol,itor+1]
    Q0=corn[ii,:,0]
    Q1=corn[ii,:,1]
    Q2=corn[ii,:,2]
    Q3=corn[ii,:,3]

    if rind==rind_plot:
      ax2.plot([P0[0],P1[0]],[P0[1],P1[1]],[P0[2],P1[2]],'r')
      ax2.plot([Q0[0],Q1[0]],[Q0[1],Q1[1]],[Q0[2],Q1[2]],'r')
      ax2.plot([Q1[0],Q2[0]],[Q1[1],Q2[1]],[Q1[2],Q2[2]],'r')
      ax2.plot([Q2[0],Q3[0]],[Q2[1],Q3[1]],[Q2[2],Q3[2]],'r')
      ax2.plot([Q3[0],Q0[0]],[Q3[1],Q0[1]],[Q3[2],Q0[2]],'r')

    #ax3[0].plot(Q0[0],Q0[1],'+',Q1[0],Q1[1],'x',Q2[0],Q2[1],'o',Q3[0],Q3[1],'.')
    #fig3.show()

    #print('Qx='+str(np.array([[Q0[0],Q1[0]],[Q3[0],Q2[0]]])))
    #print('Qy='+str(np.array([[Q0[1],Q1[1]],[Q3[1],Q2[1]]])))
    #print('Qz='+str(np.array([[Q0[2],Q1[2]],[Q3[2],Q2[2]]])))
     
    for irefine in range(Nrefinements):
      #print('---CROSSING: '+str(ii)+'-------ITERATION: '+str(irefine+1)+'------------------')
      corn[ii],thetacorn[ii],zetacorn[ii]=refine(P0,P1,corn[ii],thetacorn[ii],zetacorn[ii],
                                                 Rmn,Zmn,Dzetacylphimn,rind,ii,irefine)#,fig=fig3,ax=ax3)
      #print('RIND'+str(rind)+' CROSS'+str(ii)+' ITER'+str(irefine+1)+
      #    ', th:'+str(thetacorn[:,0])+', z:'+str(zetacorn[:,0]))
      Q0=corn[ii,:,0]
      Q1=corn[ii,:,1]
      Q2=corn[ii,:,2]
      Q3=corn[ii,:,3]
      if rind==rind_plot:
        ax2.plot([Q0[0],Q1[0]],[Q0[1],Q1[1]],[Q0[2],Q1[2]],colors[irefine])
        ax2.plot([Q1[0],Q2[0]],[Q1[1],Q2[1]],[Q1[2],Q2[2]],colors[irefine])
        ax2.plot([Q2[0],Q3[0]],[Q2[1],Q3[1]],[Q2[2],Q3[2]],colors[irefine])
        ax2.plot([Q3[0],Q0[0]],[Q3[1],Q0[1]],[Q3[2],Q0[2]],colors[irefine])
      
    if rind==rind_plot:
      ax2.set_xlim(xl)
      ax2.set_ylim(yl)
      ax2.set_zlim(zl)

    #print('RIND'+str(rind)+' CROSS'+str(ii)+' ITER'+str(irefine+1)+
    #      ', th:'+str(thetacorn[:,0])+', z:'+str(zetacorn[:,0]))
  
  #fig2.show()
  #wait_for_user()
  #MidXYZ=np.sum(corn[-1],axis=1)/4
  
  return np.array([signpdotgradr,corn[:,0,0],corn[:,1,0],corn[:,2,0],thetacorn[:,0],zetacorn[:,0]])
  

def findSightLineCrossings(Geom,P0,P1):
  R0=np.sqrt(P0[0]**2+P0[1]**2)
  R1=np.sqrt(P1[0]**2+P1[1]**2)

  #find the toroidal angles
  asin0=np.arcsin(P0[1]/R0)
  asin1=np.arcsin(P1[1]/R1)
  if P0[0]<0:
    if asin0>0:
      asin0=np.pi-asin0
    else:
      asin0=-np.pi-asin0 
  if P1[0]<0:
    if asin1>0:
      asin1=np.pi-asin1
    else:
      asin1=-np.pi-asin1 
  P0pzeta=-asin0 #PEST zeta is =-geomang
  P1pzeta=-asin1
  #print('P0pzeta='+str(P0pzeta))
  #print('P1pzeta='+str(P1pzeta))
  Nr=Geom.nsurf
  Ntheta=115
  Nzeta=117

  args=[(rind,Geom,P0,P1,Ntheta,Nzeta) for rind in range(Nr)]#Geom.nsurf)]

  Ncpu = multiprocessing.cpu_count()
  with multiprocessing.Pool(processes=Ncpu) as pool:
    out=pool.starmap(onesurf_crossings,args) #The distributed calculation

  #there can be crossings on two sides of the cross section
  #positive (p) is where grad(r) dot (P1-P0) > 0, and negative (n) is where grad(r) dot (P1-P0) < 0
  nthet=np.nan*np.zeros((Nr,))
  pthet=np.nan*np.zeros((Nr,))
  nzeta=np.nan*np.zeros((Nr,))
  pzeta=np.nan*np.zeros((Nr,))
  nXYZ =np.nan*np.zeros((Nr,3))
  pXYZ =np.nan*np.zeros((Nr,3))
  nR   =np.nan*np.zeros((Nr,))
  pR   =np.nan*np.zeros((Nr,))
  pgradrXYZ=np.nan*np.zeros((Nr,3))
  ngradrXYZ=np.nan*np.zeros((Nr,3))
  nl   =np.nan*np.zeros((Nr,))
  pl   =np.nan*np.zeros((Nr,))
  ngradrdotdgradl=np.nan*np.zeros((Nr,))
  pgradrdotdgradl=np.nan*np.zeros((Nr,))
  Blist=[None]*Nr
  zetalist=[None]*Nr

  print('looping over radii')  
  for rind in range(Nr):
    drdpsi=2*Geom.rnorm[rind]/Geom.minorradiusW7AS/Geom.psi_a
    Booz = fluxcoorddiscr(Geom,rind,Ntheta,Nzeta,name='Boozer')
    Blist[rind]=Booz.B
    zetalist[rind]=Booz.zeta
    if out[rind].size>0:
      if len(out[rind][0])>2:
        sys.exit('More than two flux surface intersections. Please choose P0 and P1 differently!')
      posind=np.where(out[rind][0]==1)[0][0]
      negind=np.where(out[rind][0]==-1)[0][0]
      if posind.size>0:
        pXYZ[rind]=out[rind][1:4,posind]
        pR[rind]=np.sqrt(pXYZ[rind][0]**2+pXYZ[rind][1]**2)
        pthet[rind]=out[rind][4,posind]
        pzeta[rind]=out[rind][5,posind]
        pl[rind]=np.sqrt(np.sum((pXYZ[rind]-P0)**2))
        pgradrXYZ[rind,0]=fluxcoorddiscr.interp2_cyclic(Booz.theta,Booz.zeta,Booz.gradpsi_X,pthet[rind],pzeta[rind],Geom.Nperiods)*drdpsi
        pgradrXYZ[rind,1]=fluxcoorddiscr.interp2_cyclic(Booz.theta,Booz.zeta,Booz.gradpsi_Y,pthet[rind],pzeta[rind],Geom.Nperiods)*drdpsi
        pgradrXYZ[rind,2]=fluxcoorddiscr.interp2_cyclic(Booz.theta,Booz.zeta,Booz.gradpsi_Z,pthet[rind],pzeta[rind],Geom.Nperiods)*drdpsi
        pgradrdotdgradl[rind]=np.sum(pgradrXYZ[rind]*(pXYZ[rind]-P0))
      if negind.size>0:
        nXYZ[rind]=out[rind][1:4,negind]
        nR[rind]=np.sqrt(nXYZ[rind][0]**2+nXYZ[rind][1]**2)
        nthet[rind]=out[rind][4,negind]
        nzeta[rind]=out[rind][5,negind]
        nl[rind]=np.sqrt(np.sum((nXYZ[rind]-P0)**2))
        ngradrXYZ[rind,0]=fluxcoorddiscr.interp2_cyclic(Booz.theta,Booz.zeta,Booz.gradpsi_X,nthet[rind],nzeta[rind],Geom.Nperiods)*drdpsi
        ngradrXYZ[rind,1]=fluxcoorddiscr.interp2_cyclic(Booz.theta,Booz.zeta,Booz.gradpsi_Y,nthet[rind],nzeta[rind],Geom.Nperiods)*drdpsi
        ngradrXYZ[rind,2]=fluxcoorddiscr.interp2_cyclic(Booz.theta,Booz.zeta,Booz.gradpsi_Z,nthet[rind],nzeta[rind],Geom.Nperiods)*drdpsi
        ngradrdotdgradl[rind]=np.sum(ngradrXYZ[rind]*(nXYZ[rind]-P0))
  return (pXYZ,pR,pthet,pzeta,pl,pgradrdotdgradl,
          nXYZ,nR,nthet,nzeta,nl,ngradrdotdgradl,
          Blist,zetalist,P0pzeta,P1pzeta)

#####################################################################
# Here the code starts
######################################################################
        
Geom=bcgeom('/afs/ipp-garching.mpg.de/home/s/smithh/Forskning/Stellarator/sfincs/gitsfincs/equilibria/w7x-sc1+252.bc')
#Geom=bcgeom('/afs/ipp-garching.mpg.de/home/s/smithh/Forskning/Stellarator/sfincs/gitsfincs/equilibria/w7x-sc1.bc',max_m=4,maxabs_n=4)
P0=[3.0,0.3,0.18] #x,y,z of starting point
P1=[7.0,0.7,0.18] #x,y,z of ending point

(pXYZ,pR,pthet,pzeta,pl,pgradrdotdgradl,
 nXYZ,nR,nthet,nzeta,nl,ngradrdotdgradl,
 Blist,zetalist,P0pzeta,P1pzeta)=findSightLineCrossings(Geom,P0,P1)

figR,axR=plt.subplots(4,1,sharex=True,figsize=(6,7))#,subplot_kw={'projection': '3d'})

axR[0].plot(pl,pthet,'b',nl,nthet,'g')
axR[1].plot(pl,pzeta,'b',nl,nzeta,'g')
axR[2].plot(pl,pgradrdotdgradl,'b',nl,ngradrdotdgradl,'g')
axR[3].plot(pl,pR,'b',nl,nR,'g')

axR[0].set_ylabel(r'$\theta$')
axR[1].set_ylabel(r'$\zeta$')
axR[2].set_ylabel(r'$\nabla r \cdot \hat{l}$')
axR[3].set_ylabel(r'$R$')

axR[3].set_xlabel(r'$l$')

#figR.show()

#inputs to crossectplot:
rinds=range(Geom.nsurf)
rNs=Geom.rnorm[rinds]
rNmin=0.0
rNmax=np.inf
print('Pest zeta values for cross section plots: '+str(P0pzeta)+', '+str(P1pzeta)) 
Pest_zetas=[P0pzeta,P1pzeta]#[0.0,np.pi/Geom.Nperiods] #toroidal angles to make plots for
rNcontours=[0.1,0.2,0.3,0.5,0.7,0.9] #flux surface contours to draw

figs,axs=crossectplot(Blist,rNs,Geom,coordname='Boozer',
                        Pest_zetas=Pest_zetas,
                        rNmin=rNmin,rNmax=rNmax,
                        rNcontours=rNcontours,
                        savefilename='figuretest.eps')
axs[0].set_title(r'Cross section at $\zeta_{Pest}=$'+'{:.4f}'.format(P0pzeta))
axs[1].set_title(r'Cross section at $\zeta_{Pest}=$'+'{:.4f}'.format(P1pzeta))
axs[0].plot(pR,pXYZ[:,2],'b',nR,nXYZ[:,2],'g')
axs[1].plot(pR,pXYZ[:,2],'b',nR,nXYZ[:,2],'g')


figs[0].show()
figs[1].show()

figR.show()

wait_for_user()
