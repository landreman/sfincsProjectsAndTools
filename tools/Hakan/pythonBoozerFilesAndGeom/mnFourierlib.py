#!/usr/bin/env python
from __future__ import division
import numpy as np
import sys, copy

import geomlib

#####################################################################
# This library includes two classes, mnmat and mnlist,
# which are used to represent the m,n Fourier coefficients on a flux
# surface as a matrix or as a list.
#
#####################################################################

###########################################################################################
###########################################################################################
###########################################################################################
class mnlist(object):
###########################################################################################
###########################################################################################
###########################################################################################

    def __init__(self,input,n=None,data=None,cosparity=1,Nperiods=None,
                 rind=None,quantity='B',vmecgrid='half'):
        #if isinstance(input,mnmat):
        # instead of this, use e.g. Bmn.mnlist()
            
        if isinstance(input,geomlib.bcgeom):
            #quantity can be R,Z,B,Dphi
            geometry=input
            if rind is None:
                sys.exit('The radius index rind is required when taking input from a bcgeom!')

            if Nperiods != None:
                if geometry.Nperiods != Nperiods:
                    sys.exit('Input Nperiods does not match the value in the bcgeom!')
                    

            self.Nperiods=int(geometry.Nperiods)
            self.m=geometry.m[rind]
            self.n=geometry.n[rind]
            data=getattr(geometry,quantity)
            if quantity=='Dphi':
                self.data=data[rind]*2*np.pi/self.Nperiods
            else:
                self.data=data[rind]
            if geometry.StelSym:
                if quantity=='B' or quantity=='R':
                    self.cosparity=np.ones(len(self.m))
                elif quantity=='Z' or quantity=='Dphi':
                    self.cosparity=np.zeros(len(self.m))
                else:
                    sys.exit('Undefined quantity!')
            else:
                if quantity=='B' or quantity=='R':
                    self.cosparity=geometry.parity[rind]
                elif quantity=='Z' or quantity=='Dphi':
                    self.cosparity=1-geometry.parity[rind]
                else:
                    sys.exit('Undefined quantity!')
            
        elif isinstance(input,geomlib.vmecgeom):
            #quantity can be R,Z,B,lambda,B_u or B_w where (u,w,s) is RH and w=-v. (u,v,s) is the LH VMEC system.
            wout=input
            if rind is None:
                sys.exit('The radius index rind is required when taking input from a bcgeom!')

            if Nperiods != None:
                if wout.Nperiods != Nperiods:
                    sys.exit('Input Nperiods does not match the value in the vmecgeom!')
                    
            signchange=wout.signgs #is -1, because vmec is left handed
            if rind is None:
                sys.exit('The radius index rind is required when taking input from a vmecgeom!')
            if vmecgrid=='full':
                sys.exit('Extracting from full vmec grid is not implemented yet!')
            skip=wout.skip #=1, this is the skip in the beginning of the half grid vmec arrays
            skrind=skip+rind
            rindf_R = skrind
            rindf_L = skrind-1
            
            self.Nperiods=int(wout.nfp)

            if wout.StelSym:
                if quantity=='B':
                    self.m   =wout.xm_nyq
                    self.n   =wout.xn_nyq*signchange//self.Nperiods #signchange because of toroidal direction swap
                    self.data=wout.bmnc[skrind]
                    self.cosparity=np.ones(len(self.m))
                elif quantity=='R':
                    self.m   =wout.xm
                    self.n   =wout.xn*signchange//self.Nperiods
                    self.data=(wout.rmnc[rindf_L]+wout.rmnc[rindf_R])/2.0
                    self.cosparity=np.ones(len(self.m))
                elif quantity=='Z':
                    self.m   =wout.xm
                    self.n   =wout.xn*signchange//self.Nperiods
                    self.data=(wout.zmns[rindf_L]+wout.zmns[rindf_R])/2.0
                    self.cosparity=np.zeros(len(self.m))
                elif quantity=='lambda':
                    self.m   =wout.xm
                    self.n   =wout.xn*signchange//self.Nperiods
                    self.data=(wout.lmns[rindf_L]+wout.lmns[rindf_R])/2.0
                    self.cosparity=np.zeros(len(self.m))
                elif quantity=='B_u':
                    self.m   =wout.xm_nyq
                    self.n   =wout.xn_nyq*signchange//self.Nperiods
                    self.data=wout.bsubumnc[skrind]
                    self.cosparity=np.ones(len(self.m))                    
                elif quantity=='B_w':
                    self.m   =wout.xm_nyq
                    self.n   =wout.xn_nyq*signchange//self.Nperiods
                    self.data=wout.bsubvmnc[skrind]*signchange
                    self.cosparity=np.ones(len(self.m))                    
                else:
                    sys.exit('Undefined quantity!')
            else:
                if quantity=='B':
                    self.m   =np.concatenate((wout.xm_nyq,wout.xm_nyq))
                    self.n   =np.concatenate((wout.xn_nyq,wout.xn_nyq))*signchange//self.Nperiods
                    self.data=np.concatenate((wout.bmnc[skrind],wout.bmns[skrind]))
                    self.cosparity=np.concatenate((np.ones(wout.mnmax_nyq),np.zeros(wout.mnmax_nyq)))
                elif quantity=='R':
                    self.m   =np.concatenate((wout.xm,wout.xm))
                    self.n   =np.concatenate((wout.xn,wout.xn))*signchange//self.Nperiods
                    self.data=np.concatenate(((wout.rmnc[rindf_L]+wout.rmnc[rindf_R])/2.0,
                                              (wout.rmns[rindf_L]+wout.rmns[rindf_R])/2.0))
                    self.cosparity=np.concatenate((np.ones(wout.mnmax),np.zeros(wout.mnmax)))
                elif quantity=='Z':
                    self.m   =np.concatenate((wout.xm,wout.xm))
                    self.n   =np.concatenate((wout.xn,wout.xn))*signchange//self.Nperiods
                    self.data=np.concatenate(((wout.zmns[rindf_L]+wout.zmns[rindf_R])/2.0,
                                              (wout.zmnc[rindf_L]+wout.zmnc[rindf_R])/2.0))                    
                    self.cosparity=np.concatenate((np.zeros(wout.mnmax),np.ones(wout.mnmax)))
                elif quantity=='lambda':
                    self.m   =np.concatenate((wout.xm,wout.xm))
                    self.n   =np.concatenate((wout.xn,wout.xn))*signchange//self.Nperiods
                    self.data=np.concatenate(((wout.lmns[rindf_L]+wout.lmns[rindf_R])/2.0,
                                              (wout.lmnc[rindf_L]+wout.lmnc[rindf_R])/2.0))
                    self.cosparity=np.concatenate((np.zeros(wout.mnmax),np.ones(wout.mnmax)))
                elif quantity=='B_u':
                    self.m   =np.concatenate((wout.xm_nyq,wout.xm_nyq))
                    self.n   =np.concatenate((wout.xn_nyq,wout.xn_nyq))*signchange//self.Nperiods
                    self.data=np.concatenate((wout.bsubumnc[skrind],wout.bsubumns[skrind]))
                    self.cosparity=np.concatenate((np.ones(wout.mnmax_nyq),np.zeros(wout.mnmax_nyq)))
                elif quantity=='B_w':
                    self.m   =np.concatenate((wout.xm_nyq,wout.xm_nyq))
                    self.n   =np.concatenate((wout.xn_nyq,wout.xn_nyq))*signchange//self.Nperiods
                    self.data=np.concatenate((wout.bsubvmnc[skrind],wout.bsubvmns[skrind]))*signchange
                    self.cosparity=np.concatenate((np.ones(wout.mnmax_nyq),np.zeros(wout.mnmax_nyq)))
                else:
                    sys.exit('Undefined quantity!')
            self.m = self.m.astype(int)
            self.n = self.n.astype(int)

        #if isinstance(input,mnmat.mnmat):
        #    return input.mnlist()
        
        else: #input is supposed to be the m array
            m=input
            if Nperiods is None:
                self.Nperiods=None
            else: 
                self.Nperiods=int(Nperiods)
            if len(m)!=len(n):
                sys.exit("m and n have different lengths")
            if len(n)!=len(data):
                sys.exit("n and data have different lengths")
            if isinstance(cosparity,int):
                cosparity=float(cosparity)
            if isinstance(cosparity,float):
                self.cosparity=np.ones(len(m))*cosparity
            elif len(cosparity)!=len(m):
                sys.exit("m and cosparity have different lengths")
            else:
                self.cosparity=np.array(cosparity)
            self.m=np.array(m)
            self.n=np.array(n)
            self.data=np.array(data)

                
    def disp(self,threshold=-1.0):
        datafrm='{:8.4e}'
        modefrm='{:5d}'
        print('------------------------------------')
        if not(self.Nperiods is None):
            print('Nperiods='+str(self.Nperiods))
        print('parity  m    n    data')
        for ind in range(len(self.m)):
            if abs(self.data[ind])>threshold:
                if self.cosparity[ind]==1.0:
                    cos_sin_str='cos '
                else:
                    cos_sin_str='sin '
                print(cos_sin_str+modefrm.format(self.m[ind])+
                      modefrm.format(self.n[ind])+'    '+datafrm.format(self.data[ind]))
        print('------------------------------------')

###########################################################################################
###########################################################################################
###########################################################################################
class mnmat(object):
###########################################################################################
###########################################################################################
###########################################################################################
    
    def __init__(self,input,Ntheta=None,Nzeta=None,Nperiods=None,rind=None,
                 quantity='B',vmecgrid='half',warnings='off'):
        #Note that rind and quantity are only needed if input is a bcgeom

        if isinstance(input, float) or isinstance(input, int):
            if isinstance(input, int):
                input=float(input)
            if Ntheta is None:
                Ntheta=5 #default small number
            if Nzeta is None:
                Nzeta=5 #default small number
            if Ntheta%2 != 1 or Nzeta%2 != 1:
                sys.exit('sizes must be odd')
            Ntheta=int(Ntheta)
            Nzeta=int(Nzeta)
            self.Ntheta=Ntheta
            self.Nzeta=Nzeta

            maxm=(Ntheta-1)//2
            Nm=maxm+1
            maxabsn=(Nzeta-1)//2
            Nn=Nzeta
            [self.m,self.n]=np.mgrid[0:Nm,-maxabsn:maxabsn+1]
            #print self.m
            #print self.n

            
            m0ind=0
            n0ind=(Nzeta-1)//2
            self.m0ind=m0ind
            self.n0ind=n0ind
            Nperiods=int(Nperiods)
            self.Nperiods=Nperiods
            
            self.c=np.zeros((Nm,Nn))
            self.s=np.zeros((Nm,Nn))
            self.c[0,:n0ind]=np.nan
            self.s[0,:n0ind+1]=np.nan
            self.c[0,n0ind]=input
            #print('float input: '+str(self.c))

        elif isinstance(input, np.ndarray): #This is the fft routine, Nperiods is required input
            do_expand=False
            Ntheta_size=input.shape[0]
            newNtheta=Ntheta
            Ntheta=Ntheta_size
            if newNtheta is None:
                newNtheta=Ntheta_size
            elif newNtheta<Ntheta_size:
                print('Ntheta size of input: '+str(Ntheta_size))
                print('Ntheta parameter: '+str(newNtheta))
                sys.exit('Ntheta is smaller than input size in mnmat fft routine for ndarrays! Not implemented yet.')
            elif newNtheta>Ntheta_size:
                newNtheta=int(newNtheta)
                do_expand=True
            
            Nzeta_size=input.shape[1]
            newNzeta=Nzeta
            Nzeta=Nzeta_size
            if newNzeta is None:
                newNzeta=Nzeta_size
            elif newNzeta<Nzeta_size:
                print('Nzeta size of input: '+str(Nzeta_size))
                print('Nzeta parameter: '+str(newNzeta))
                sys.exit('Nzeta is smaller than input size in mnmat fft routine for ndarrays! Not implemented yet.')
            elif newNzeta>Nzeta_size:
                newNzeta =int(newNzeta)
                do_expand=True
            
            if Ntheta%2 != 1 or Nzeta%2 != 1:
                sys.exit('Sizes Ntheta and Nzeta of the input matrix must be odd!')
            self.Ntheta=Ntheta
            self.Nzeta=Nzeta
            if newNtheta%2 != 1 or newNzeta%2 != 1:
                sys.exit('New sizes Ntheta and Nzeta must be odd!')

            if Nperiods is None:
                sys.exit('Nperiods is needed, not for the fft but for later use!')
            Nperiods=int(Nperiods)
            
            maxm=(Ntheta-1)//2
            Nm=maxm+1
            maxabsn=(Nzeta-1)//2
            Nn=Nzeta
            self.m,self.n=np.mgrid[0:Nm,-maxabsn:maxabsn+1]
            
            self.m0ind=0
            self.n0ind=(Nzeta-1)//2
            self.Nperiods=Nperiods
            maxm=(self.Ntheta-1)//2
            n0ind=self.n0ind

            Fflip=np.fft.fftshift(np.fft.fft2(input),1)

            self.c =  2.0/(Ntheta*Nzeta) * np.real(Fflip[0:maxm+1,::-1]);
            self.s = -2.0/(Ntheta*Nzeta) * np.imag(Fflip[0:maxm+1,::-1]);

            self.c[0,:n0ind]=np.nan
            self.s[0,:n0ind+1]=np.nan
            self.c[0,n0ind]=self.c[0,n0ind]/2.0

            if do_expand:
                tmp = mnmat(self.mnlist(),newNtheta,newNzeta)
                self.Ntheta   = tmp.Ntheta
                self.Nzeta    = tmp.Nzeta
                self.m0ind    = tmp.m0ind
                self.n0ind    = tmp.n0ind
                self.Nperiods = tmp.Nperiods
                self.m        = tmp.m
                self.n        = tmp.n
                self.c        = tmp.c
                self.s        = tmp.s

        elif isinstance(input, mnlist):
            if input.Nperiods is None:
                if Nperiods is None:
                    sys.exit('Nperiods is missing!')
                else:    
                    self.Nperiods=int(Nperiods)
            elif not(Nperiods is None):
                if Nperiods==input.Nperiods:
                    self.Nperiods=int(Nperiods)
                else:
                    sys.exit('You are trying to change the number of field periods'+
                             'from '+str(input.Nperiods)+' to '+str(Nperiods)+'!')
            else:
                self.Nperiods=int(input.Nperiods)
            
            if Ntheta is None:
                Ntheta=int(max(input.m)*2+1)
            if Nzeta is None:
                Nzeta=int(max(abs(input.n))*2+1)
            if Ntheta%2 != 1 or Nzeta%2 != 1:
                sys.exit('sizes must be odd')
            self.Ntheta=Ntheta
            self.Nzeta=Nzeta
            maxm=int((Ntheta-1)//2)
            Nm=maxm+1
            maxabsn=int((Nzeta-1)//2)
            Nn=Nzeta
            self.m,self.n=np.mgrid[0:Nm,-maxabsn:maxabsn+1]
            
            m0ind=0
            n0ind=int((Nzeta-1)//2)
            self.m0ind=m0ind
            self.n0ind=n0ind
            self.c=np.zeros((Nm,Nn))
            self.s=np.zeros((Nm,Nn))

            if warnings=='on':
                if (input.m>maxm).any():
                    Ntheta_needed=int(max(input.m)*2+1)
                    sys.exit('Increase Ntheta to '+str(Ntheta_needed))
                if (abs(input.n)>maxabsn).any():
                    Nzeta_needed=int(max(abs(input.n))*2+1)
                    sys.exit('Increase Nzeta to '+str(Nzeta_needed))

            if (input.cosparity==1).all():
                #print(len(input.m))
                for ind in range(len(input.m)):
                    #print('ind='+str(ind)+', m='+str(input.m[ind])+', n='+str(input.n[ind]))
                    if input.m[ind]<=maxm and abs(input.n[ind])<=maxabsn:
                        if input.m[ind]==0 and input.n[ind]<=0:
                            self.c[0,n0ind-input.n[ind]]+=input.data[ind]            # (*)
                        else:
                            self.c[input.m[ind],n0ind+input.n[ind]]+=input.data[ind] #+= necessary because of (*)
            else:
                for ind in range(len(input.m)):
                    if input.m[ind]<=maxm and abs(input.n[ind])<=maxabsn:
                        if input.m[ind]==0 and input.n[ind]<0 and input.cosparity[ind]==1:
                            self.c[0,n0ind-input.n[ind]]+=input.data[ind]
                        elif input.m[ind]==0 and input.n[ind]<0 and input.cosparity[ind]==0:
                            self.s[0,n0ind-input.n[ind]]-=input.data[ind]
                        elif input.cosparity[ind]==1:
                            self.c[input.m[ind],n0ind+input.n[ind]]+=input.data[ind]
                        else: #sinus component
                            self.s[input.m[ind],n0ind+input.n[ind]]+=input.data[ind]
            self.c[0,:n0ind]=None
            self.s[0,:n0ind+1]=None
            #print('list input:')
            #print(self.c)
            #print(self.s)
        elif isinstance(input, geomlib.bcgeom):
            lista=mnlist(input,rind=rind,quantity=quantity)
            tmp=mnmat(lista,Ntheta,Nzeta)
            self.Nperiods=int(tmp.Nperiods)
            self.Ntheta  =int(tmp.Ntheta)
            self.Nzeta   =int(tmp.Nzeta)

            maxm=int((self.Ntheta-1)//2)
            Nm=maxm+1
            maxabsn=int((self.Nzeta-1)//2)
            Nn=self.Nzeta
            self.m,self.n=np.mgrid[0:Nm,-maxabsn:maxabsn+1]
             
            self.m0ind   =tmp.m0ind
            self.n0ind   =tmp.n0ind
            self.c=tmp.c
            self.s=tmp.s

        elif isinstance(input, geomlib.vmecgeom):
            #quantity can be R,Z,B,lambda,B_u or B_w where (u,w,s) is RH and w=-v. (u,v,s) is the LH VMEC system.
            lista=mnlist(input,rind=rind,quantity=quantity,vmecgrid=vmecgrid)
            #print('shape of lista: '+str(lista.m.shape))
            #print('Ntheta='+str(Ntheta))
            #print('Nzeta='+str(Nzeta))

            tmp=mnmat(lista,Ntheta=Ntheta,Nzeta=Nzeta)
            self.Nperiods=int(tmp.Nperiods)
            self.Ntheta  =int(tmp.Ntheta)
            self.Nzeta   =int(tmp.Nzeta)

            maxm=int((self.Ntheta-1)//2)
            Nm=maxm+1
            maxabsn=int((self.Nzeta-1)//2)
            Nn=self.Nzeta
            self.m,self.n=np.mgrid[0:Nm,-maxabsn:maxabsn+1]
            
            self.m0ind   =tmp.m0ind
            self.n0ind   =tmp.n0ind
            self.c=tmp.c
            self.s=tmp.s
            



                
    #-----------------------------------------------------------------------
    # Other functions
    #-----------------------------------------------------------------------
    def __add__(self,other):
        if isinstance(other, mnmat):
            summa=mnmat(0,Ntheta=self.Ntheta,Nzeta=self.Nzeta,Nperiods=self.Nperiods)
            summa.c=self.c+other.c            
            summa.s=self.s+other.s            
            return summa

    def __sub__(self,other):
        if isinstance(other, mnmat):
            difference=mnmat(0,Ntheta=self.Ntheta,Nzeta=self.Nzeta,Nperiods=self.Nperiods)
            difference.c=self.c-other.c            
            difference.s=self.s-other.s            
            return difference
        
    def __neg__(self):
        difference=mnmat(0,Ntheta=self.Ntheta,Nzeta=self.Nzeta,Nperiods=self.Nperiods)
        difference.c=-self.c
        difference.s=-self.s
        return difference

    def __mul__(self,other):
        if isinstance(other, mnmat):
            sys.exit('Cannot multiply two mnmats')
        out=copy.deepcopy(self)
        out.c=self.c*other
        out.s=self.s*other
        return out
        
    def __rmul__(self,other):
        if isinstance(other, mnmat):
            sys.exit('Cannot multiply two mnmats')
        out=copy.deepcopy(self)
        out.c=self.c*other
        out.s=self.s*other
        return out

    def expand(self,Ntheta,Nzeta,copy='True'): #Pad with zeros to a new larger size
        if Ntheta==self.Ntheta and Nzeta==self.Nzeta:
            if copy:
                return copy.deepcopy(self)
        elif Ntheta<self.Ntheta or Nzeta<self.Nzeta:
            sys.exit('In expand, the new Ntheta and Nzeta cannot be smaller than the originals!')
        else:
            if copy:
                return mnmat(self.mnlist(),Ntheta,Nzeta)
            else:
                tmp = mnmat(self.mnlist(),Ntheta,Nzeta)
                self.Ntheta   = tmp.Ntheta
                self.Nzeta    = tmp.Nzeta
                self.m0ind    = tmp.m0ind
                self.n0ind    = tmp.n0ind
                self.Nperiods = tmp.Nperiods
                self.m        = tmp.m
                self.n        = tmp.n
                self.c        = tmp.c
                self.s        = tmp.s


    def mnlist(self):
        maxm=int((self.Ntheta-1)//2)
        Nn=self.Nzeta
        Ns=int(Nn*maxm+(Nn-1)//2)
        Nc=Ns+1
        m_flat=self.m.flatten()
        n_flat=self.n.flatten()
        c_flat=self.c.flatten()
        s_flat=self.s.flatten()
        m=np.concatenate((m_flat[self.n0ind:],m_flat[self.n0ind+1:]))
        n=np.concatenate((n_flat[self.n0ind:],n_flat[self.n0ind+1:]))
        data=np.concatenate((c_flat[self.n0ind:],s_flat[self.n0ind+1:]))
        cosparity=np.concatenate((np.ones(Nc),np.zeros(Ns)))
        return mnlist(m,n,data,cosparity,self.Nperiods)
        
    def disp(self):
        if not(self.Nperiods is None):
            print('Nperiods='+str(self.Nperiods))
        print('------------ cosinus -------------')
        print(self.c)
        print('------------- sinus --------------')
        print(self.s)

    def evalpoint(self,u,v):
        if isinstance(u,int):
            u=float(u) 
        if isinstance(v,int):
            v=float(v) 
        if np.isscalar(u) and np.isscalar(v):
            c=self.c*np.cos(self.m * u - self.n * self.Nperiods * v)
            s=self.s*np.sin(self.m * u - self.n * self.Nperiods * v)
            c[0,:self.n0ind]=0
            s[0,:self.n0ind+1]=0
            return np.sum(c+s)
        if u.ndim==0:
            c=self.c*np.cos(self.m * u - self.n * self.Nperiods * v)
            s=self.s*np.sin(self.m * u - self.n * self.Nperiods * v)
            c[0,:self.n0ind]=0
            s[0,:self.n0ind+1]=0
            return np.sum(c+s)
        if u.ndim>1:
            ufl=u.flatten()
            vfl=v.flatten()
        elif u.ndim==1:
            ufl=u
            vfl=v
            
        c=self.c*np.cos(ufl[:,np.newaxis,np.newaxis]*self.m - vfl[:,np.newaxis,np.newaxis]*self.n*self.Nperiods)
        s=self.s*np.sin(ufl[:,np.newaxis,np.newaxis]*self.m - vfl[:,np.newaxis,np.newaxis]*self.n*self.Nperiods)

        c[:,0,:self.n0ind]=0
        s[:,0,:self.n0ind+1]=0
        
        return np.sum((c+s),axis=(1,2)).reshape(u.shape)
       
    def ifft(self):
        maxm=int((self.Ntheta-1)//2)
        n0ind=self.n0ind
        Nm=maxm+1
        Nn=self.Nzeta

        Fmncos=np.zeros((self.Ntheta,self.Nzeta))
        Fmnsin=np.zeros((self.Ntheta,self.Nzeta))

        Fmncos[0:maxm+1,:]=self.c
        Fmncos[0,n0ind] =Fmncos[0,n0ind]*2.0
        Fmncos[0,:n0ind]=Fmncos[0,:n0ind:-1]
        Fmnsin[0:maxm+1,:]=self.s
        Fmnsin[0,n0ind] =0.0
        Fmnsin[0,:n0ind]=-Fmnsin[0,:n0ind:-1]
        
        #Fmncos[maxm+1:,-1]=Fmncos[maxm:0:-1,-1]
        Fmncos[maxm+1:,n0ind]= Fmncos[maxm:0:-1,n0ind]
        Fmnsin[maxm+1:,n0ind]=-Fmnsin[maxm:0:-1,n0ind]

        Fmncos[maxm+1:,0:n0ind]  = Fmncos[maxm:0:-1,:n0ind:-1]
        Fmncos[maxm+1:,n0ind+1:] = Fmncos[maxm:0:-1,n0ind-1::-1]
        Fmnsin[maxm+1:,0:n0ind]  =-Fmnsin[maxm:0:-1,:n0ind:-1]
        Fmnsin[maxm+1:,n0ind+1:] =-Fmnsin[maxm:0:-1,n0ind-1::-1]

        F=self.Ntheta*self.Nzeta/2.0*(Fmncos-1j*Fmnsin);
        #print(np.real(F))

        f=np.real(np.fft.ifft2(np.fft.ifftshift(F[:,::-1],1)));

        return f

    def grad(self):
        #%Let us force the size of the data matrix to be the same
        #%In principle one should increase it because there is more sine 
        #%information than it can accommodate now, but that would be unpractical
        #%if we want to transform back to real space.

        dduFmn=copy.deepcopy(self)
        ddvFmn=copy.deepcopy(self)
        
        dduFmn.s = -self.c*self.m
        dduFmn.c =  self.s*self.m
        dduFmn.c[dduFmn.m0ind,dduFmn.n0ind]=0
        dduFmn.s[dduFmn.m0ind,dduFmn.n0ind]=np.nan
  
        ddvFmn.s =  self.c*self.n*self.Nperiods
        ddvFmn.c = -self.s*self.n*self.Nperiods
        ddvFmn.c[ddvFmn.m0ind,ddvFmn.n0ind]=0
        ddvFmn.s[ddvFmn.m0ind,ddvFmn.n0ind]=np.nan
 
        return dduFmn, ddvFmn

    def calcu(self,G,I,iota,radialvariable='poloidal flux',
              zeroout_Deltaiota=-1.0, reference_theta=None,reference_zeta=None):
        # This solves the equation
        #
        #   B\cdot\nabla u = - B\times\nabla\chi\cdot\nabla h
        #
        # given that <B^2 u> = 0
        #
        # Very often we are interested to solve the equation with 
        # h=1/B^2, in which case one puts fftmn(1./B.^2) as the 
        # first argument. For bootstrap coefficient calculations
        # one instead often solves the equation with h=1/sqrt(1-lambda B).
        #
        # The optional argument zeroout_Deltaiota can be used 
        # to zero out mn components of wmn
        # which have abs(n/m*NPeriods - iota) < zeroout_Deltaiota.
        #
        # The optional argument radialvariablecan be used to change from the
        # poloidal flux \chi as the flux surface label in the above
        # equation to the toroidal flux \psi by giving the strings
        # 'poloidal flux' (default) or 'toroidal flux'.
        #
        # The optional arguments reference_theta,reference_zeta can be used to set
        # u_00 so that u(reference_theta,reference_zeta)=0. If a reference point
        # is not given, then u_00=0.
        
        if radialvariable=='poloidal flux':
            iotaexp=1
        elif radialvariable=='toroidal flux':
            iotaexp=0
        else:
            sys.exit('Input radialvariable to calcu invalid!')

        umn=copy.deepcopy(self)
        N=self.Nperiods

        with np.errstate(invalid='ignore'): #don't warn about division by nan
            umn.c=iota**iotaexp*(G*self.m + I*self.n * N)/(self.n * N - iota*self.m)*self.c
            umn.s=iota**iotaexp*(G*self.m + I*self.n * N)/(self.n * N - iota*self.m)*self.s
        
        umn.c[umn.m0ind,umn.n0ind]=0
        umn.s[umn.m0ind,umn.n0ind]=np.nan

        m_tmp=np.ndarray.astype(umn.m,'float')
        m_tmp[umn.m0ind,:]=1e-10 #to avoid division by zero
        umn.c = np.where(np.logical_and(umn.m>0,abs(umn.n/m_tmp*N - iota)<zeroout_Deltaiota),0,umn.c)
        umn.s = np.where(np.logical_and(umn.m>0,abs(umn.n/m_tmp*N - iota)<zeroout_Deltaiota),0,umn.s)

        if not(reference_theta is None) and not(reference_zeta is None):
            umn.set00(-umn.evalpoint(reference_theta,reference_zeta))
                
        return umn
        
    def invJacBdotgrad(self,iota):
        #%Assumes Boozer coordinates!

        outmn=copy.deepcopy(self)
        with np.errstate(divide='ignore'): #don't warn about division by 0
           outmn.c = -self.s/(iota*self.m-self.Nperiods*self.n)
           outmn.s =  self.c/(iota*self.m-self.Nperiods*self.n)

        outmn.c[outmn.m0ind,outmn.n0ind]=0
        outmn.s[outmn.m0ind,outmn.n0ind]=np.nan

        return outmn

    def invgrad(self,ddvGmn,method=1):
        Gmn=copy.deepcopy(self)
        Gmn.c[:,:]=0
        Gmn.s[:,:]=0
        dduGmn=self
        if dduGmn.c[self.m0ind,self.n0ind]!=0:
            sys.exit('There is a 00 component in an mnmat to be integrated!')
        if ddvGmn.c[self.m0ind,self.n0ind]!=0:
            sys.exit('There is a 00 component in an mnmat to be integrated!')

        #Two options that should be equivalent
        if method==1:
            Gmn.s[1:,:] =  dduGmn.c[1:,:]/dduGmn.m[1:,:];
            Gmn.c[1:,:] = -dduGmn.s[1:,:]/dduGmn.m[1:,:];

            Gmn.s[0,Gmn.n0ind+1:] = -ddvGmn.c[0,Gmn.n0ind+1:]/[ddvGmn.n[0,Gmn.n0ind+1:]*Gmn.Nperiods];
            Gmn.c[0,Gmn.n0ind+1:] =  ddvGmn.s[0,Gmn.n0ind+1:]/[ddvGmn.n[0,Gmn.n0ind+1:]*Gmn.Nperiods];
        else:
            Gmn.s[:,Gmn.n0ind+1:] = -ddvGmn.c[:,Gmn.n0ind+1:]/[ddvGmn.n[:,Gmn.n0ind+1:]*Gmn.Nperiods];
            Gmn.c[:,Gmn.n0ind+1:] =  ddvGmn.s[:,Gmn.n0ind+1:]/[ddvGmn.n[:,Gmn.n0ind+1:]*Gmn.Nperiods];
            Gmn.s[1:,:Gmn.n0ind] = -ddvGmn.c[1:,:Gmn.n0ind]/[ddvGmn.n[1:,:Gmn.n0ind]*Gmn.Nperiods];
            Gmn.c[1:,:Gmn.n0ind] =  ddvGmn.s[1:,:Gmn.n0ind]/[ddvGmn.n[1:,:Gmn.n0ind]*Gmn.Nperiods];
    
            Gmn.s[1:,Gmn.n0ind] =  dduGmn.c[1:,Gmn.n0ind]/dduGmn.m[1:,Gmn.n0ind];
            Gmn.c[1:,Gmn.n0ind] = -dduGmn.s[1:,Gmn.n0ind]/dduGmn.m[1:,Gmn.n0ind];    
            
            
        return Gmn

    def get00(self):
        return self.c[self.m0ind,self.n0ind]
    
    def set00(self,value):
        self.c[self.m0ind,self.n0ind]=value
            
    def remove00(self,copy=False):
        if copy:
            outmn=copy.deepcopy(self)
            outmn.c[outmn.m0ind,outmn.n0ind]=0
            return outmn
        else:
            self.c[self.m0ind,self.n0ind]=0
            return self
