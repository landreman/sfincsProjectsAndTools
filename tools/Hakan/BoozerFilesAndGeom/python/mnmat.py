#!/usr/bin/env python
import numpy as np
import sys, copy
from mnlist import mnlist
from bcgeom import bcgeom

class mnmat:
    
    def __init__(self,input,Ntheta=None,Nzeta=None,Nperiods=None,rind=None,quantity='B'):
        #Note that rind and quantity are only needed if input is a bcgeom
        if isinstance(input, int):
            input=float(input)
            
        if isinstance(input, float):
            if Ntheta is None:
                Ntheta=5 #default small number
            if Nzeta is None:
                Nzeta=5 #default small number
            if Ntheta%2==0 or Nzeta%2==0:
                sys.exit('sizes must be odd')
            self.Ntheta=Ntheta
            self.Nzeta=Nzeta

            maxm=(Ntheta-1)/2
            Nm=maxm+1
            maxabsn=(Nzeta-1)/2
            Nn=Nzeta
            [self.m,self.n]=np.mgrid[0:Nm,-maxabsn:maxabsn+1]
            #print self.m
            #print self.n 
            
            m0ind=0
            n0ind=(Nzeta-1)/2
            self.m0ind=m0ind
            self.n0ind=n0ind
            self.Nperiods=Nperiods
            
            self.c=np.zeros((Nm,Nn))
            self.s=np.zeros((Nm,Nn))
            self.c[0,n0ind]=input
            print('float input: '+str(self.c))

        elif isinstance(input, np.ndarray): #This is the fft routine, Nperiods is required input
            Ntheta_size=input.shape[0]
            if Ntheta is None:
                Ntheta=Ntheta_size
            elif Ntheta_size!=Ntheta:
                sys.exit('Not implemented yet')
            else:
                Ntheta=Ntheta_size
            self.Ntheta=Ntheta
            
            Nzeta_size=input.shape[1]
            if Nzeta is None:
                Nzeta=Nzeta_size
            elif Nzeta_size!=Nzeta:
                sys.exit('Not implemented yet')
            else:
                Nzeta=Nzeta_size
            self.Nzeta=Nzeta

            if Nperiods is None:
                 sys.exit('Nperiods is needed, not for the fft but for later use!')
            
            maxm=(Ntheta-1)/2
            Nm=maxm+1
            maxabsn=(Nzeta-1)/2
            Nn=Nzeta
            self.m,self.n=np.mgrid[0:Nm,-maxabsn:maxabsn+1]
            
            self.m0ind=0
            self.n0ind=(Nzeta-1)/2
            self.Nperiods=Nperiods
            maxm=(self.Ntheta-1)/2
            n0ind=self.n0ind

            Fflip=np.fft.fftshift(np.fft.fft2(input),1)

            self.c =  2.0/(Ntheta*Nzeta) * np.real(Fflip[0:maxm+1,::-1]);
            self.s = -2.0/(Ntheta*Nzeta) * np.imag(Fflip[0:maxm+1,::-1]);

            self.c[0,0:n0ind]=None
            self.s[0,0:n0ind+1]=None
            self.c[0,n0ind]=self.c[0,n0ind]/2.0
            
            #print('ndarray input:')
            #print(self.c)
            #print(self.s)

        elif isinstance(input, mnlist):
            if input.Nperiods is None:
                self.Nperiods=Nperiods
            elif not(Nperiods is None):
                if Nperiods==input.Nperiods:
                    self.Nperiods=Nperiods
                else:
                    sys.exit('You are trying to change the number of field periods'+
                             'from '+str(input.Nperiods)+' to '+str(Nperiods)+'!')
            else:
                self.Nperiods=input.Nperiods
            
            if Ntheta is None:
                Ntheta=max(input.m)*2+1
            if Nzeta is None:
                Nzeta=max(abs(input.n))*2+1
            if Ntheta%2==0 or Nzeta%2==0:
                sys.exit('sizes must be odd')
            self.Ntheta=Ntheta
            self.Nzeta=Nzeta
            maxm=(Ntheta-1)/2
            Nm=maxm+1
            maxabsn=(Nzeta-1)/2
            Nn=Nzeta
            self.m,self.n=np.mgrid[0:Nm,-maxabsn:maxabsn+1]
            
            m0ind=0
            n0ind=(Nzeta-1)/2
            self.m0ind=m0ind
            self.n0ind=n0ind
            self.c=np.zeros((Nm,Nn))
            self.s=np.zeros((Nm,Nn))

            if (input.m>maxm).any():
                sys.exit('Increase Ntheta')
            if (abs(input.n)>maxabsn).any():
                sys.exit('Increase Nzeta')
            
            if input.cosparity.all():
                #print(len(input.m))
                for ind in range(len(input.m)):
                    #print('ind='+str(ind)+', m='+str(input.m[ind])+', n='+str(input.n[ind]))
                    if input.m[ind]==0 and input.n[ind]<=0:
                        self.c[0,n0ind-input.n[ind]]=input.data[ind]
                    else:
                        self.c[input.m[ind],n0ind+input.n[ind]]=input.data[ind]
            else:
                for ind in range(len(input.m)):
                    if input.m[ind]==0 and input.n[ind]<0 and input.cosparity[ind]==1:
                        self.c[0,n0ind-input.n[ind]]=input.data[ind]
                    elif input.m[ind]==0 and input.n[ind]<0 and input.cosparity[ind]==0:
                        self.s[0,n0ind-input.n[ind]]=-input.data[ind]
                    elif input.cosparity[ind]==1:
                        self.c[input.m[ind],n0ind+input.n[ind]]=input.data[ind]
                    else: #sinus component
                        self.s[input.m[ind],n0ind+input.n[ind]]=input.data[ind]
            #print('list input:')
            #print(self.c)
            #print(self.s)
        elif isinstance(input, bcgeom):
            lista=mnlist(input,rind=rind,quantity=quantity)
            tmp=mnmat(lista,Ntheta,Nzeta)
            self.Nperiods=tmp.Nperiods
            self.Ntheta  =tmp.Ntheta
            self.Nzeta   =tmp.Nzeta

            maxm=(self.Ntheta-1)/2
            Nm=maxm+1
            maxabsn=(self.Nzeta-1)/2
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
            summa=mnmat(0)
            summa.c=self.c+other.c            
            summa.s=self.c+other.s            
            return summa

    def disp(self):
        if not(self.Nperiods is None):
            print 'Nperiods='+str(self.Nperiods)
        print '------------ cosinus -------------'
        print self.c
        print '------------- sinus --------------'
        print self.s

    def evalpoint(self,u,v):
        if u.ndim==0:
            c=self.c*np.cos(self.m * u - self.n * self.Nperiods * v)
            s=self.s*np.sin(self.m * u - self.n * self.Nperiods * v)
            c[0,:self.n0ind]=0
            s[0,:self.n0ind+1]=0
            return (c+s).sum
        if u.ndim>1:
            ufl=u.flatten()
            vfl=v.flatten()
        elif u.ndim==1:
            ufl=u
            vfl=v
            
        c=self.c*np.cos(ufl[:,np.newaxis,np.newaxis]*self.m - vfl[:,np.newaxis,np.newaxis]*self.n*self.Nperiods)
        s=self.s*np.sin(ufl[:,np.newaxis,np.newaxis]*self.m - vfl[:,np.newaxis,np.newaxis]*self.n*self.Nperiods)
        
        return (c+s).sum((1,2)).reshape(u.shape)
       
    def ifft(self):
        maxm=(self.Ntheta-1)/2
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
        
        Fmncos[maxm+1:,-1]=Fmncos[maxm:0:-1,-1]
        Fmnsin[maxm+1:,n0ind]=-Fmnsin[maxm:0:-1,n0ind]

        Fmncos[maxm+1:,0:n0ind]  = Fmncos[maxm:0:-1,:n0ind:-1]
        Fmncos[maxm+1:,n0ind+1:] = Fmncos[maxm:0:-1,n0ind-1::-1]
        Fmnsin[maxm+1:,0:n0ind]  =-Fmnsin[maxm:0:-1,:n0ind:-1]
        Fmnsin[maxm+1:,n0ind+1:] =-Fmnsin[maxm:0:-1,n0ind-1::-1]

        F=self.Ntheta*self.Nzeta/2.0*(Fmncos-1j*Fmnsin);
        #print(F)

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

    def calcu(self,G,I,iota,radialvariable='poloidal flux',zeroout_Deltaiota=-1.0):
        #% This solves the equation
        #%
        #%   B\cdot\nabla u = - B\times\nabla\chi\cdot\nabla h
        #%
        #% given that <B^2 u> = 0
        #%
        #% Very often we are interested to solve the equation with 
        #% h=1/B^2, in which case one puts fftmn(1./B.^2) as the 
        #% first argument.
        #%
        #% The optional 6th argument zeroout_Deltaiota can be used 
        #% to zero out mn components of wmn
        #% which have abs(n/m*NPeriods - iota) < zeroout_Deltaiota.
        #%
        #% The optional 7th argument can be used to change from the
        #% poloidal flux \chi as the flux surface label in the above
        #% equation to the toroidal flux \psi by giving the strings
        #% 'poloidal flux' (default) or 'toroidal flux'.
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
