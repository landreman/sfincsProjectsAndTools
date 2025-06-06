#!/usr/bin/env python
from __future__ import division
import numpy as np
import sys, os, copy
import datetime

import mnFourierlib
import fluxcoorddiscr
import netCDF4
from netCDF4 import Dataset
#import timeit
import scipy.integrate as integrate

#####################################################################
# This library includes two main classes, bcgeom and vmecgeom,
# which are used to load and analyse Boozer .bc files and VMEC wout
# netcdf files.
#
#####################################################################

##########################################################################################
##########################################################################################
######################################## bcgeom ##########################################
##########################################################################################
##########################################################################################
#
# Simplest ways to call the initiation routine:
# Geom=bcgeom('examplefile.bc')
# Geom=bcgeom('examplefile.dat')
# Geom=bcgeom(wout), where wout is a vmecgeom object
#

class headertxt(object):
  def __init__(self):
      self.maincomment=''
      self.globalvars=''
      self.surfvars=''
      self.surfvarunits=''
      self.datavars=''

class Bfiltr(object):
  def __init__(self):
      self.min_Bmn=0
      self.max_m=np.inf
      self.maxabs_n=np.inf

class bcgeom(object):
    
  def __init__(self,input,min_Bmn=0,max_m=np.inf,maxabs_n=np.inf,
               symmetry='unknown',signcorr=2,verbose=1):
    if isinstance(input,netCDF4.Group) or isinstance(input,netCDF4.Dataset):
        readfromnc(input,self)
        if not(hasattr(self,'nsurf')):
            sys.exit('Could not load bcgeom')
    elif isinstance(input,str):
       filename=input
       # This reads the Boozer file of the .bc type
       # Note that the .bc file is left-handed and the output struct self
       # is right-handed.
       #
       # Only Fourier components with abs(Bmn)>min_Bmn, m<=max_m, |n|<=maxabs_n
       # are read.
       #
       # The input symmetry can be 'StelSym' to double-check that the 
       # input is stellarator symmetric
       #
       # If the .bc file was produced by JMC, one should supply an extra argument
       # choosing which type of correction of the sign of the total flux to make.
       # signcorr=1: Bphi and Btheta will get sign changes (as in SFINCS)
       # signcorr=2: Total toroidal flux will get a sign change (default, as Yuriy Turkin)

       def sscan(strng,dtype=float,count=-1):
         if isinstance(strng,list):
            out=[None]*len(strng)
            for i in range(len(strng)):
                out[i]=np.fromstring(strng[i],dtype,count,sep=' ')
            return tuple(out)  
         else:    
            return np.fromstring(strng,dtype,count,sep=' ')

       if signcorr==1:
           newsigncorrectionmethod=False
       elif signcorr==2:
           newsigncorrectionmethod=True
       else:
           sys.exit('Sign correction method for JG files not recognised!')


       if filename[-4:]=='.dat':
           #Henning Maassberg type file
           filetype='HM'
       elif filename[-3:]=='.bc':
           #Joachim Geiger (or Erika Strumberger) type file
           filetype='JG'
       else: #default to JG
           filetype='JG'

       if filetype=='JG': 
         with open(filename, 'r') as f:
           if sys.version_info[0] == 2:
               f.seek(-1,2)     # go to the file end.
           else:
               f.seek(0, os.SEEK_END)
               f.seek(f.tell() - 1, os.SEEK_SET)

           eof = f.tell()   # get the end of file location
           f.seek(0,0)      # go back to file beginning

           YTsign=1 #1 means no sign change, 
                   #-1 means sign change because YT has resaved JG's file
           tmp_str=f.readline()
           concat_str=tmp_str
           self.headertext=headertxt()
           if tmp_str[0:2]=='CC':
             if 'CStconfig' in tmp_str: 
                 YTsign=-1
             while tmp_str[0:2]=='CC':
                 tmp_str=f.readline()
                 if tmp_str[0:2]=='CC': #comment line
                     if concat_str[-1]=='\n':
                         concat_str = concat_str+tmp_str 
                     else:
                         concat_str = concat_str+'\n'+tmp_str 
                 if 'CStconfig' in tmp_str:
                     YTsign=-1
             if concat_str[-1]=='\n':
                 concat_str=concat_str[:-1]
             self.headertext.maincomment=concat_str
           else:
             self.headertext.maincomment=('CC ----------------------'+
                                        '------------------------')

           self.headertext.globalvars=tmp_str

     
           header_df=sscan(f.readline())
           self.headertext.surfvars=f.readline() #Variable name line  
           YTstyle=0
           if '[A]' in self.headertext.surfvars:
             #This is a file from Yuriy Turkin. Correct this line to the JG standard
             YTstyle=1 #indicates Yuriy Turkin style
             self.headertext.surfvars=(
               '       s         iota  curr_pol/nper    curr_tor    pprime   sqrt g(0,0)')
             self.headertext.surfvarunits=(
               '                            [A]            [A]   dp/ds,[Pa] (dV/ds)/nper')  
           self.m0b        = header_df[0]
           self.n0b        = header_df[1]
           self.nsurf      = int(header_df[2])
           self.Nperiods   = int(header_df[3])
           self.psi_a=np.nan #Insert psi_a at this place in the list, but set it later.
           self.torfluxtot = header_df[4]*YTsign #Note that this is not per pol. angle,
                                           #note: possible YTsign
           self.minorradiusW7AS      = header_df[5]
           self.majorradiusLastbcR00 = header_df[6]
           if len(header_df)>7:
             self.minorradiusVMEC=header_df[7]
           else:
             self.minorradiusVMEC=np.nan
           if len(header_df)>8:
             self.majorradiusVMEC=header_df[8]
           else:
             self.majorradiusVMEC=np.nan
           if len(header_df)>9:
             self.volumeVMEC=header_df[9]
           else:
             self.volumeVMEC=np.nan
         
           endoffile=False
           rind=-1
   
           torfluxnorm=np.array([])
           iota=np.array([])
           Bphi=np.array([])
           Btheta=np.array([])
           dpds=np.array([])
           dVdsoverNper=np.array([])
           B00=np.array([])
           R00=np.array([])
           no_of_modes=[]
           modesm=[]
           modesn=[]
           modesr=[]
           modesz=[]
           modesp=[]
           modesb=[]
           modesbnorm=[]
           modespar=[]

      
           while not endoffile:
             rind=rind+1
             if verbose>0:
               sys.stdout.write('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%5i/%5i'% (rind+1,self.nsurf))
             if not(YTstyle): #isempty(strfind(self.headertext.surfvars,'[A]'))
                 self.headertext.surfvarunits=f.readline() #unit line only in JG files
             #print(surfvarunits: '+self.headertext.surfvarunits)
             surfheader=sscan(f.readline()) #fscanf(fid,'%f',6)
             torfluxnorm=np.append(torfluxnorm, surfheader[0])
             iota =np.append(iota,surfheader[1])
             Bphi=np.append(Bphi,surfheader[2]*self.Nperiods/2.0/np.pi*(4*np.pi*1e-7)) #Tesla*meter
             Btheta=np.append(Btheta,surfheader[3]/2.0/np.pi*(4*np.pi*1e-7))              #Tesla*meter
             dpds=np.append(dpds,surfheader[4])
             dVdsoverNper=np.append(dVdsoverNper,surfheader[5])

             #f.readline() #just skip the return character
             tmpstrunits=f.readline() #units line
             if rind==0:
                 self.headertext.datavars=tmpstrunits
                 self.StelSym=True
                 if 'rmnc' in self.headertext.datavars:
                     self.StelSym=False #('rmnc' in self.headertext.datavars)
                     
                 if not(symmetry=='unknown'):
                     if self.StelSym and not(symmetry=='StelSym'):
                         sys.exit('Boozer file is stellarator symmetric, but input to'+
                                  'readBoozerfile said it should be non-stellarator symmetric!')
                     elif not(self.StelSym) and symmetry=='StelSym':
                         sys.exit('Boozer file is non-stellarator symmetric,'+
                                  'but input to readBoozerfile said it should be'+
                                  'stellarator symmetric!')
                    
             position=f.tell()
             if self.StelSym:
                 tmp=sscan(f.readline()) #,'%d %d %f %f %f %f',6)
                 while not(tmp[0]==0 and tmp[1]==0):
                     tmp=sscan(f.readline(),count=6) #sscanf(tmp_str1,'%d %d %f %f %f %f',6)    
                 B00=np.append(B00,tmp[5])
                 R00=np.append(R00,tmp[2])
             else:
                 tmp=sscan(f.readline())#sscanf(tmp_str1,'%d %d %f %f %f %f %f %f %f %f',10)
                 while not(tmp[0]==0 and tmp[1]==0):
                     tmp=sscan(f.readline()) #sscanf(tmp_str1,'%d %d %f %f %f %f %f %f %f %f',10)
         
                 B00=np.append(B00,tmp[8])
                 R00=np.append(R00,tmp[2])

             f.seek(position) #Rewind to beginning of surface data
       
             modesm.append(np.array([]))
             modesn.append(np.array([]))
             modesr.append(np.array([]))
             modesz.append(np.array([]))
             modesp.append(np.array([]))
             modesb.append(np.array([]))
             modesbnorm.append(np.array([]))
             modespar.append(np.array([]))

             proceed=True
             modeind=-1
             while proceed:
                 tmp_str=f.readline()
                 if f.tell() == eof+1:
                     #print('eof found '+tmp_str)
                     proceed=False
                     endoffile=True
                     #print('found end of file')
                 if ('s' in tmp_str): #Next flux surface has been reached
                     proceed=False
                     #print('found next surface')
                 else:
                     #print(tmp_str)
                     if self.StelSym:
                         tmp=sscan(tmp_str,count=6)
                         if ((abs(tmp[5])/B00[rind]>min_Bmn) and
                             (tmp[0]<=max_m) and (abs(tmp[1])<=maxabs_n)):
                             modeind=modeind+1
                             modesm[rind]=np.append(modesm[rind],tmp[0])
                             modesn[rind]=np.append(modesn[rind],tmp[1])
                             modesr[rind]=np.append(modesr[rind],tmp[2])
                             modesz[rind]=np.append(modesz[rind],tmp[3])
                             modesp[rind]=np.append(modesp[rind],tmp[4])
                             modesb[rind]=np.append(modesb[rind],tmp[5])
                             modesbnorm[rind]=np.append(modesbnorm[rind],tmp[5]/B00[rind])

                     else:
                         tmp=sscan(tmp_str,count=10)
                         if (tmp[0]<=max_m) and (abs(tmp[1])<=maxabs_n):
                             if (abs(tmp[8])/B00[rind]>min_Bmn):
                                 #Cosinus component
                                 modeind=modeind+1
                                 modesm[rind]=np.append(modesm[rind],tmp[0])
                                 modesn[rind]=np.append(modesn[rind],tmp[1])
                                 modesr[rind]=np.append(modesr[rind],tmp[2])
                                 modesz[rind]=np.append(modesz[rind],tmp[5])
                                 modesp[rind]=np.append(modesp[rind],tmp[7])
                                 modesb[rind]=np.append(modesb[rind],tmp[8])
                                 modesbnorm[rind]=np.append(modesbnorm[rind],tmp[8]/B00[rind])
                                 modespar[rind]=np.append(modespar[rind],1) #parity 1 <=> Cosinus component

                             if  (abs(tmp[9])/B00[rind]>min_Bmn):
                                 #Sinus component
                                 modeind=modeind+1
                                 modesm[rind]=np.append(modesm[rind],tmp[0])
                                 modesn[rind]=np.append(modesn[rind],tmp[1])
                                 modesr[rind]=np.append(modesr[rind],tmp[3])
                                 modesz[rind]=np.append(modesz[rind],tmp[4])
                                 modesp[rind]=np.append(modesp[rind],tmp[6])
                                 modesb[rind]=np.append(modesb[rind],tmp[9])
                                 modesbnorm[rind]=np.append(modesbnorm[rind],tmp[9]/B00[rind])
                                 modespar[rind]=np.append(modespar[rind],0) #parity 0 <=> Sinus component
                             #end
                         #end
                     #end
                 #end
             #end while proceed
             if modeind==-1:
                sys.exit('no modes found for rind='+str(rind))

             no_of_modes.append(int(modeind+1))

           #end while loop over radii
         #end with open(filename, 'r') as f
         if verbose>0:
           print('') #go to new line
         
         if any([a>0 for a in dVdsoverNper]):
           print('The coordinate system in the Boozer '+
                    'file should be left handed, but it has '+
                    'a positive Jacobian. Something is wrong! '+
                    'I have had this for LHD bc files before. '+
                    'The problem was perhaps introduced by Nikolai. '+
                    'The sign will be changed!')
           dVdsoverNper=-np.abs(dVdsoverNper)

         if self.torfluxtot*Bphi[0]>0:
           if not(self.StelSym):
               sys.exit('Unknown sign convention in non-stellarator-symmetric file!')
           if not(newsigncorrectionmethod):
               self.newsigncorr=0
               if not(self.StelSym):
                   sys.exit('Rotating 180 degrees only allowed for '+
                            'stellarator symmetric cases!')

               Bphi=-Bphi
               Btheta=-Btheta
           else: #Use newsigncorrectionmethod
               self.newsigncorr=1
               self.torfluxtot=-self.torfluxtot
         elif self.StelSym: #(and self.torfluxtot*Bphi[0]<0)
           sys.exit('Unknown sign convention in non-stellarator-symmetric file! '+
                    'Neither YT nor JG standard.')

         rthetazeta_righthanded=np.sign(self.torfluxtot*Bphi[0])
         #This is -1, because the Boozer file is supposed to be left handed
         if rthetazeta_righthanded==1:
           sys.exit('The coordinate system in the Boozer file was right handed')

         self.torfluxtot=self.torfluxtot*rthetazeta_righthanded  
         self.psi_a=self.torfluxtot/2.0/np.pi
         self.rnorm=np.sqrt(torfluxnorm)
         self.s=torfluxnorm
         self.Bphi=Bphi*rthetazeta_righthanded**2 #amperes law minus sign and direction switch sign
         self.Btheta=Btheta*rthetazeta_righthanded #amperes law minus sign
         self.iota=iota*rthetazeta_righthanded
         self.dpds=dpds
         if self.dpds[-1]==0: #Joachim sets the last point to zero for some reason
             self.dpds[-1]=2*self.dpds[-2]-self.dpds[-3] #Extrapolate!

         self.dVdsoverNper=dVdsoverNper*rthetazeta_righthanded
         self.FSAB2=abs(4*np.pi**2*self.psi_a/self.Nperiods*
                        (self.Bphi+self.iota*self.Btheta)/self.dVdsoverNper)
         self.nmodes=no_of_modes
         self.m=[e.astype(int) for e in modesm]
         self.n=[e.astype(int) for e in modesn] #sign is switched below
         self.B=modesb
         self.Bnorm=modesbnorm
         self.B00=B00
         self.Bfilter=Bfiltr()
         self.Bfilter.min_Bmn  = min_Bmn
         self.Bfilter.max_m    = max_m
         self.Bfilter.maxabs_n = maxabs_n
         self.R00=R00
         #%self.Z00=zeros(size(R00))
         self.R=modesr
         self.Z=modesz
         self.Dphi=modesp;
         #%for rind=1:self.nsurf
         #%  m0n0ind=find(self.m{rind}==0 & self.n{rind}==0);
         #%  self.Z00(rind)=self.Z{rind}(m0n0ind);
         #%This is always zero (also in Erika's files)

         if rthetazeta_righthanded==-1:
           for tmpind in range(len(self.n)):
               self.n[tmpind]=-self.n[tmpind]  #This assumes the argument  
                                               #(m theta - n N phi) in (r,theta,phi) left-handed 
                                               #Correct the Boozerfile first if another
                                               #convention was used there
               self.Dphi[tmpind]=-self.Dphi[tmpind] 
               #Dphi =N/(2*pi)*(zeta-(-geomang))

         if not(self.StelSym):
           self.parity=modespar

           #If non-stellarator symmetric, then assume that it is from Erika.
           #Then the minor and major radii are VMEC and not Joachims definitions
           self.minorradiusVMEC=self.minorradiusW7AS
           self.majorradiusVMEC=self.majorradiusLastbcR00
           self.volumeVMEC=np.pi*self.minorradiusVMEC**2*2*np.pi*self.majorradiusVMEC
           self.minorradiusW7AS=np.nan
           self.majorradiusLastbcR00=np.nan


         #%%%% Some security checks 
         #Ntheta=301;Nzeta=303;rind=1;
         #Bmn=mnmat(self,rind,'B',Ntheta,Nzeta,'forceSize');
         #B=ifftmn(Bmn);
         #FSAB2_test=Ntheta*Nzeta/sum(sum(B.^(-2)));
         #if abs(FSAB2_test./self.FSAB2(rind)-1)>0.05
         #  warning('dVdsoverNper not correct in bc file')
         #  self.dVdsoverNper=NaN*self.dVdsoverNper;
         #  self.FSAB2=NaN*self.FSAB2;
         #  abs(4*pi^2*self.psi_a/self.Nperiods*(self.Bphi+self.iota.*self.Btheta)./self.dVdsoverNper);
         #  %Calculate the correct values
         #   for rind=1:self.nsurf
         #    Bmn=mnmat(self,rind,'B',Ntheta,Nzeta,'forceSize');
         #    B=ifftmn(Bmn);
         #    self.FSAB2(rind)=Ntheta*Nzeta/sum(sum(B.^(-2)));
         #  end
         #  self.dVdsoverNper=abs(4*pi^2*self.psi_a/self.Nperiods*(self.Bphi+self.iota.*self.Btheta)./self.FSAB2);
         #end

       #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       #% Henning Maassberg type Boozer file  
       #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       if filetype=='HM': 
           self.StelSym=1
           self.newsigncorr=newsigncorrectionmethod

           with open(filename) as f:
               flines = f.readlines()
           
           #f = open(filename, 'r')
           #f.seek(-1,2)     # go to the file end.
           #eof = f.tell()   # get the end of file location
           #f.seek(0,0)      # go back to file beginning
           ind=0
           tmp_str=flines[ind]
           ind+=1
           while tmp_str[1]=='c' or tmp_str[1]=='C':  #Skip comment line
               tmp_str=flines[ind]#tmp_str=f.readline() 
               ind+=1
           torfluxnorm=np.array([])
           radii=np.array([])
           iota=np.array([])
           Bphi=np.array([])
           Btheta=np.array([])
           dpds=np.array([])
           dVdsoverNper=np.array([])
           B00=np.array([])
           R00=np.array([])
           no_of_modes=[]
           modesm=[]
           modesn=[]
           modesr=[]
           modesz=[]
           modesp=[]
           modesb=[]
           modesbnorm=[]
           modesdbdrnorm=[]
           modespar=[]
             
           endoffile=False
           rind=-1

           while not endoffile:
               rind=rind+1
               surfheader1     = sscan(tmp_str,count=8)
               radii           = np.append(radii,surfheader1[0]/100) #(cm->m)
               iota            = np.append(iota,surfheader1[1])
               Nperiods        = int(surfheader1[2])
               minorradiusW7AS = surfheader1[3]/100 #(cm)
               majorradiusLastbcR00 = surfheader1[4]/100 #(cm)
               if len(surfheader1)>5:
                   torfluxnorm = np.append(torfluxnorm,surfheader1[5])
                   R00         = np.append(R00,surfheader1[7]/100)
               else:
                   torfluxnorm = np.append(torfluxnorm,np.nan)
                   R00         = np.append(R00,np.nan)

                   
               tmp_str=flines[ind]
               ind+=1
               if tmp_str[1]=='>':
                   surfheader2=sscan(tmp_str[2:],count=3)
           
                   Bphi      = np.append(Bphi,surfheader2[0]*1.0e6*Nperiods/2.0/np.pi*(4*np.pi*1e-7)) #Tesla*meter
                   Btheta    = np.append(Btheta,surfheader2[1]*1.0e6/2.0/np.pi*(4*np.pi*1e-7))        #Tesla*meter
                   B00       = np.append(B00,surfheader2[2]) #Tesla
                   tmp_str=flines[ind]
                   ind+=1
               else:
                   Bphi   = np.append(Bphi,np.nan)
                   Btheta = np.append(Btheta,np.nan)
                   B00    = np.append(B00,np.nan)


               modesm.append(np.array([]))
               modesn.append(np.array([]))
               modesr.append(np.array([]))
               modesz.append(np.array([]))
               modesp.append(np.array([]))
               modesb.append(np.array([]))
               modesbnorm.append(np.array([]))
               modesdbdrnorm.append(np.array([]))
               modespar.append(np.array([]))
    
               proceed=True
               modeind=-1
           
               while proceed:
                   if sscan(tmp_str,count=1)==-1: #Next flux surface reached
                       proceed=False
                   else:
                       tmp=sscan(tmp_str,count=7)
                       if ((abs(tmp[2])>min_Bmn)and(tmp[0]<=max_m)and
                           (abs(tmp[1])<=maxabs_n)):
                           modeind=modeind+1
                           modesm[rind]       = np.append(modesm[rind],tmp[0])
                           modesn[rind]       = np.append(modesn[rind],tmp[1])
                           modesbnorm[rind]   = np.append(modesbnorm[rind],tmp[2])
                           modesb[rind]       = np.append(modesb[rind],tmp[3])
                           modesdbdrnorm[rind]= np.append(modesdbdrnorm[rind],tmp[4]*100) #conv. cm^-1 to m^-1
                           modesr[rind]= np.append(modesr[rind],tmp[5])
                           modesz[rind]= np.append(modesz[rind],tmp[6])

                   if ind<len(flines):
                       tmp_str=flines[ind]
                       ind+=1 #get the next line or surface header line
               #end while proceed
               if modeind==-1:
                   sys.exit('no modes found for rind='+str(rind))

               no_of_modes.append(int(modeind+1))

               if ind>=len(flines): #f.tell() == eof:
                   endoffile=True

           #end while not endoffile
           #f.close()
           
           #self.torfluxtot is not stored in Henning Maassberg's files. Because they are 
           #based on Joachim Geiger's .bc files, however, they are left-handed (r,pol,tor) and
           #inconsistent with the sign of torfluxtot in the .bc file. The same correction is
           #therefore made here as for the .bc files above.
           rthetazeta_righthanded=-1
     
           self.nsurf=len(radii)
           self.Nperiods=Nperiods
           self.minorradiusW7AS=minorradiusW7AS
           self.majorradiusLastbcR00=majorradiusLastbcR00
           self.minorradiusVMEC = np.nan
           self.majorradiusVMEC = np.nan
           self.volumeVMEC      = np.nan
           self.m0b             = np.nan
           self.n0b             = np.nan
           self.torfluxtot      = np.nan
           self.rnorm=radii/self.minorradiusW7AS
           self.s=torfluxnorm
           self.iota=iota*rthetazeta_righthanded
           self.Bphi=Bphi
           self.Btheta=Btheta*rthetazeta_righthanded
           self.nmodes=no_of_modes
           self.m=[e.astype(int) for e in modesm]
           self.n=[e.astype(int) for e in modesn] #sign is switched below
           self.B=modesb
           self.Bnorm=modesbnorm
           self.B00=B00
           self.Bfilter=Bfiltr()
           self.Bfilter.min_Bmn  = min_Bmn
           self.Bfilter.max_m    = max_m
           self.Bfilter.maxabs_n = maxabs_n
           self.R00=R00
           self.R=modesr
           self.Z=modesz
           if rthetazeta_righthanded==-1:
               for tmpind in range(len(self.n)):
                   self.n[tmpind]=-self.n[tmpind]
              
           if not(self.newsigncorr):
               self.Bphi=-self.Bphi
               self.Btheta=-self.Btheta
               #NB: toroidal flux is not stored in HM's .dat files.
               
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    #% Turn a loaded vmec dataset vmecgeom into a bcgeom
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    elif isinstance(input,vmecgeom):
        wout=input
        signchange=float(wout.signgs) #is -1, because vmec is left handed
        self.StelSym=input.StelSym
        self.headertext=headertxt()
        # Note that the following are only comments to what would appear in
        # a file stored with .write()
        # The internal representation in python is right handed and psi_a>0
        # in the direction of phi!
        if self.StelSym:
            # Use Joachim Geiger's JMC convention as default for stellarator symmetric files:
            # Toroidal flux is counted positive in the direction opposite to the toroidal coordinate.
            self.headertext.maincomment=('CC Converted from VMEC using HS''s python bcgeom.py routines\n'+
               'CC version_       = '+str(wout.version_)+'\n'+
               'CC input_extension= '+wout.input_extension+'\n'+
               'CC The phase convention (m theta - n phi) is used in this file.\n'+
               'CC The coordinate system (r,theta,phi) is left-handed.\n'+
               'CC Toroidal flux is counted positive in the direction opposite to '+
               'the toroidal coordinate (JG style).')
        else:
            # Use Erika Strumberger's convention as default for non-stellarator symmetric files:
            # Toroidal flux is counted positive in the direction of the toroidal coordinate.
            self.headertext.maincomment=('CC Converted from VMEC using HS''s python routines\n'+
            'CC VMEC version         = '+str(wout.version_)+'\n'+
            'CC VMEC input extension = '+wout.input_extension+'\n'+
            'CC The phase convention (m theta - n phi) is used in this file.\n'+
            'CC The coordinate system (r,theta,phi) is left-handed.\n'+
            'CC Toroidal flux is counted positive in the direction of the toroidal coordinate (ES,YT style).')
          
        self.m0b=wout.mpol
        self.n0b=wout.ntor
        self.nsurf= np.nan      #number of radial surfaces, Set this below
        self.Nperiods=int(wout.nfp)   #%!< number of field periods
        self.torfluxtot = wout.phi[wout.ns-1]*signchange
        self.psi_a=self.torfluxtot/2.0/np.pi
        self.minorradiusVMEC      = wout.Aminor_p  #minor plasma radius
        self.majorradiusLastbcR00 = np.nan         #Calculate this below (not necessary)
        self.minorradiusW7AS      = np.nan         #Calculate this below (not necessary)
        self.majorradiusVMEC      = wout.Rmajor_p  #major plasma radius
        self.volumeVMEC           = wout.volume_p  #plasma volume
        
        fullgrid_s           = np.array(wout.phi/wout.phi[wout.ns-1]) #full grid
        fullgrid_rnorm       = np.sqrt(fullgrid_s)     #full grid

        skip=wout.skip #=1,this is how many elements are skipped at low radii when going to half grid

        self.s      = (fullgrid_s[skip-1:-1]+fullgrid_s[skip:])/2.0  #half grid
        self.rnorm  = np.sqrt(self.s) #half grid
        self.nsurf  = len(self.s)
        self.dpds   = np.array(np.diff(wout.presf[skip-1:])/np.diff(fullgrid_s[skip-1:]))
        self.Bphi        = wout.bvco[skip:]*signchange #direction switch sign
        self.Btheta      = wout.buco[skip:] #sign change

        if self.StelSym:
            self.newsigncorr=True #default signcorr=2
            if signcorr==1:
                self.newsigncorr=False
            elif signcorr==2:
                self.newsigncorr=True #default signcorr=2
            if not(self.newsigncorr):
                self.Bphi      = -self.Bphi
                self.Btheta    = -self.Btheta
                self.torfluxtot= -self.torfluxtot
                self.psi_a     = -self.psi_a
                
        fullgrid_iota=wout.iotaf*signchange
        iota   = wout.iotas*signchange #half mesh
        self.iota=iota[skip:]
        
        self.Bfilter=Bfiltr()
        self.Bfilter.min_Bmn  = min_Bmn

        #Use a RH coordinate system (u,w,s) with w=-v. Here, (u,v) are the VMEC coordinates

        if max_m==np.inf:
          Ntheta = int(1+2*max(abs(wout.xm)))
          self.Bfilter.max_m    = int((Ntheta-1)//2)
        else:
          Ntheta = int(max_m*2+1)
          self.Bfilter.max_m    = int(max_m)
          
        if maxabs_n==np.inf:
          Nzeta=1+2*max(abs(wout.xn))//self.Nperiods
          self.Bfilter.maxabs_n = (Nzeta-1)//2
        else:
          Nzeta = int(maxabs_n*2+1)
          self.Bfilter.maxabs_n = int(maxabs_n)

        self.dVdsoverNper = np.zeros(len(self.s)) 
        self.B00          = np.zeros(len(self.s))
        self.R00          = np.zeros(len(self.s))
        self.nmodes       = [0] * len(self.s)  #np.zeros(len(self.s))
        self.FSAB2        = np.zeros(len(self.s))
        self.m=[]
        self.n=[]
        if not(self.StelSym):
            self.parity=[]
        self.B=[]
        self.R=[]
        self.Z=[]
        self.Bnorm=[]
        self.Dphi=[]

        #print('Nzeta,Ntheta = '+str(Nzeta)+', '+str(Ntheta)
        if verbose>0:
          print('Converting VMEC to Boozer coordinates...')
        for rind in range(len(self.s)):
            if verbose>0:
              sys.stdout.write('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bRadius%5i/%5i'% (rind+1,len(self.s)))
            Booz=fluxcoorddiscr.fluxcoorddiscr(wout,rind=rind,Npol=Ntheta,Ntor=Nzeta,name='Boozer')
            #(fig, ax)=Booz.plot('B')
            #fig.show()
            #print(Booz.B.shape)
            self.B00[rind]=Booz.B00
            self.R00[rind]=Booz.R00
            if self.StelSym:
                self.nmodes[rind]=int((Ntheta*Nzeta+1)//2)
            else:
                self.nmodes[rind]=int(Ntheta*Nzeta)
            self.dVdsoverNper[rind]=4*np.pi**2/self.Nperiods*abs(
                self.psi_a*(self.Bphi[rind]+self.iota[rind]*self.Btheta[rind]))/Booz.FSAB2

            self.FSAB2[rind]=Booz.FSAB2
            self.m.append(np.array([]))
            self.n.append(np.array([]))
            self.B.append(np.array([]))
            self.Bnorm.append(np.array([]))
            self.R.append(np.array([]))
            self.Z.append(np.array([0]))    #m=n=0 cos element is 0
            self.Dphi.append(np.array([0])) #m=n=0 cos element is 0
            if not(self.StelSym):
                self.parity.append(np.array([]))

            Blist=(mnFourierlib.mnmat(Booz.B,Nperiods=self.Nperiods)).mnlist()
            Rlist=(mnFourierlib.mnmat(Booz.R,Nperiods=self.Nperiods)).mnlist()
            Zlist=(mnFourierlib.mnmat(Booz.Z,Nperiods=self.Nperiods)).mnlist()
            Dphilist=(mnFourierlib.mnmat(Booz.Dzetapzeta*self.Nperiods/(2*np.pi),Nperiods=self.Nperiods)).mnlist()

            #NOTA BENE: The following works because Bmn.mnlist() gives first cos, then sin components of Bmn
            if self.StelSym:
                self.m[rind]   =np.append(self.m[rind],Blist.m[:(Ntheta*Nzeta+1)//2])
                self.n[rind]   =np.append(self.n[rind],Blist.n[:(Ntheta*Nzeta+1)//2])
                self.B[rind]   =np.append(self.B[rind],Blist.data[:(Ntheta*Nzeta+1)//2])
                self.R[rind]   =np.append(self.R[rind],Rlist.data[:(Ntheta*Nzeta+1)//2])
                self.Z[rind]   =np.append(self.Z[rind],Zlist.data[(Ntheta*Nzeta+1)//2:])
                self.Dphi[rind]=np.append(self.Dphi[rind],Dphilist.data[(Ntheta*Nzeta+1)//2:])
            else:
                self.m[rind]=np.append(self.m[rind],Blist.m)
                self.n[rind]=np.append(self.n[rind],Blist.n)
                self.parity[rind]=np.append(self.parity[rind],Blist.cosparity)
                self.B[rind]=np.append(self.B[rind],Blist.data)
                self.R[rind]=np.append(self.R[rind],Rlist.data)
                self.Z[rind]=np.append(self.Z[rind],np.concatenate((Zlist.data[(Ntheta*Nzeta+1)//2:],Zlist.data[1:(Ntheta*Nzeta+1)//2])))
                self.Dphi[rind]=np.append(self.Dphi[rind],np.concatenate((Dphilist.data[(Ntheta*Nzeta+1)//2:],Dphilist.data[1:(Ntheta*Nzeta+1)//2])))
                   
            self.Bnorm[rind]=np.append(self.Bnorm[rind],self.B[rind]/Booz.B00)
            #self.m=[e.astype(int) for e in modesm]
            #self.n=[e.astype(int) for e in modesn]

        #end radius loop
        if verbose>0:
          print('') #carriage return

        self.dVdsoverNper=np.abs(self.psi_a*4*np.pi**2.0/self.Nperiods*
                                 (self.Bphi+self.iota*self.Btheta)/self.FSAB2)

        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #% Calculate minorradiusW7AS
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        accum=0;
        for m in range(1,wout.xm[-1]+1):
            ii=[i for i,x in enumerate(wout.xm) if x==m]
            first_mode=ii[0]
            last_mode =ii[-1]
            ns=np.expand_dims(wout.xn[ii[0]:ii[-1]+1]//wout.nfp, axis=1) #row vector
            nnmat=(1.0+(-1.0)**(ns-ns.T))/2.0
            accum=accum+m*np.sum((np.expand_dims(wout.rmnc[-1][ii[0]:ii[-1]+1], axis=1)*
                                  np.expand_dims(wout.zmns[-1][ii[0]:ii[-1]+1], axis=0))*nnmat)

        self.minorradiusW7AS=np.sqrt(np.abs(accum))
        
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Calculate majorradiusLastbcR00
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Joachim defined his major radius to be
        # R00 in Boozer coordinates at s=0.995
        # which in his case is at the last half mesh radius
        # I choose to always take it at the last half mesh radius
        self.majorradiusLastbcR00=self.R00[-1]

  ##############################################################################
  # Calculate Volume-averaged beta
  ##############################################################################
  def averagebeta(self):
      p=integrate.cumtrapz(self.dpds,self.s,initial=0.0)
      p=p-p[-1]
      Ntheta=141 #Default values hardcoded here
      Nzeta=131  #Default values hardcoded here
      FSAinvB2=np.zeros(self.nsurf)
      
      for rind in range(self.nsurf):
        B=mnFourierlib.mnmat(self,rind=rind,Ntheta=Ntheta,Nzeta=Nzeta,quantity='B').ifft()
        FSAinvB2[rind]=np.sum(B**(-4))/(Ntheta*Nzeta)*self.FSAB2[rind]

      mu0=4*np.pi*1e-7
      FSAbeta=p*FSAinvB2*2.0*mu0
      beta=np.trapz(FSAbeta*self.dVdsoverNper,self.s)/np.trapz(self.dVdsoverNper,self.s)
      if beta<5.0e-4: #The calculation is rather imprecise
          beta=0.0
      return beta

  ##############################################################################
  # filter a bcgeom
  ##############################################################################
  def filter(self,varinput,max_m=np.inf,maxabs_n=np.inf,makecopy=True):
      #This can be called in the following ways:
      #out=self.filter(min_Bmn)
      #out=self.filter(min_Bmn,max_m,maxabs_n)
      #out=self.filter(Bfilter)
      # if makecopy=True the output is a new instance, and the input is untouched

      if makecopy:
         out=copy.deepcopy(self)
      
      if isinstance(varinput,Bfiltr):
          min_Bmn  = max(varinput.min_Bmn,self.Bfilter.min_Bmn)
          max_m    = min(varinput.max_m,self.Bfilter.max_m)
          maxabs_n = min(varinput.maxabs_n,self.Bfilter.maxabs_n)
      else:
          min_Bmn  = max(varinput,self.Bfilter.min_Bmn)
          max_m    = min(max_m,self.Bfilter.max_m)
          maxabs_n = min(maxabs_n,self.Bfilter.maxabs_n)

      if max_m < np.inf:
          max_m=int(max_m)
      if maxabs_n < np.inf:
          maxabs_n=int(maxabs_n)
          
      if makecopy:
          out.Bfilter.min_Bmn  = min_Bmn
          out.Bfilter.max_m    = max_m
          out.Bfilter.maxabs_n = maxabs_n
      else:
          self.Bfilter.min_Bmn  = min_Bmn
          self.Bfilter.max_m    = max_m
          self.Bfilter.maxabs_n = maxabs_n

          
      for rind in range(self.nsurf):
          inds=np.where(np.logical_and(np.abs(self.Bnorm[rind])>=min_Bmn,
                                       np.logical_and(self.m[rind]<=max_m,
                                                      np.abs(self.n[rind])<=maxabs_n)))[0]
          if makecopy:
              out.nmodes[rind]=len(inds)
              out.m[rind]=self.m[rind][inds]
              out.n[rind]=self.n[rind][inds]
              out.B[rind]=self.B[rind][inds]
              out.R[rind]=self.R[rind][inds]
              out.Z[rind]=self.Z[rind][inds]
              out.Bnorm[rind]=self.Bnorm[rind][inds]
              out.Dphi[rind]=self.Dphi[rind][inds]
              if not(self.StelSym):
                  out.parity[rind]=self.parity[rind][inds]

              #print(self.Bnorm[rind])
              #print(out.Bnorm[rind])
              #print(out.m[rind])
              #print(out.n[rind])
              #sys.exit('stop')
          else:
              self.nmodes[rind]=len(inds)
              self.m[rind]=self.m[rind][inds]
              self.n[rind]=self.n[rind][inds]
              self.B[rind]=self.B[rind][inds]
              self.R[rind]=self.R[rind][inds]
              self.Z[rind]=self.Z[rind][inds]
              self.Bnorm[rind]=self.Bnorm[rind][inds]
              self.Dphi[rind]=self.Dphi[rind][inds]
              if not(self.StelSym):
                  self.parity[rind]=self.parity[rind][inds]
      if makecopy:
          return out
      else:
          return self
      
  ##############################################################################
  # Write a bcgeom to a file
  ##############################################################################
  # One can choose between
  # minorradiusconvention='W7AS' (default) Joachim Geiger's convention
  # minorradiusconvention='VMEC'
  # and between
  # majorradiusconvention='lastR00' (default) Joachim Geiger's convention
  # majorradiusconvention='VMEC'
  def write(self,filename,Nradii=None,min_Bmn=0.0,nsortstyle=None,printheadercomment=True,
            minorradiusconvention='W7AS',majorradiusconvention='lastR00',HMstandardradii=None):
      #argument Nradii is only used for .dat files.
      #First check which type of file the user wants.
      if filename[-3:]=='.bc':
          filetype='bc'
      elif filename[-4:]=='.dat':
          filetype='dat'
          #190503: Changed so that Nradii=None now writes all available radii instead of the HM 7 choice.
          ##if Nradii is None:
          ##    Nradii=7 #7 is the default for HM's DKES runs
      else: #if no file extension is given, assume .bc
          filetype='bc' 
          filename=filename+'.bc'

      if filetype=='bc':
          if min_Bmn==0:
              selfie=copy.copy(self)
          else:
              selfie=self.filter(min_Bmn)
      else: #.dat
          #190503: Changed so that Nradii=None now writes all available radii instead of the HM 7 choice.
          if HMstandardradii is None:
            if Nradii is None or Nradii==self.nsurf:
              HMstandardradii=False
            else:       
              HMstandardradii=True
          if not(HMstandardradii):
              selfie=self.interp(self.s,'s','linear').filter(min_Bmn) #This is just to retrieve dBdrnorm
          else: # Nradii != self.nsurf:
              print("Nradii != nsurf. Interpolating to Nradii='+str(Nradii)+' radii according to HM's standard.")
              #HS's choice:
              #rnorms=np.linspace(self.rnorm[0],self.rnorm[-1],num=Nradii+1,endpoint=False)
              #HM's choice:
              rnorms=np.linspace(0.0,1.0,Nradii,endpoint=False)+1.0/Nradii/2.0
              ss=rnorms**2.0
              ss[ 0] = max(ss[ 0],self.s[ 0])
              ss[-1] = min(ss[-1],self.s[-1])
              if np.any(np.diff(ss)<0.0):
                  print('s='+str(ss))
                  print('s is non-monotonic!')
                  sys.exit('HM style choice of '+str(Nradii)+' radii failed. Please reduce Nradii!')
              selfie=self.interp(ss,'s','linear').filter(min_Bmn)
              
      f = open(filename, 'w')
      now = datetime.datetime.now()
      
      signchange=-1 #sign changer
      psi_a      = selfie.psi_a*signchange
      torfluxtot = selfie.torfluxtot*signchange
      Bphi=selfie.Bphi
      Btheta=selfie.Btheta*signchange
      if selfie.StelSym:
         if not(selfie.newsigncorr):
            Bphi   = -Bphi
            Btheta = -Btheta
         else:
            torfluxtot=-torfluxtot

      #print 'In write: torfluxtot='+str(torfluxtot)
      iota=selfie.iota*signchange
      dVdsoverNper=selfie.dVdsoverNper*signchange
      Dphi=[]
      Z=[]
      m=selfie.m
      n=[]
      for tmprind in range(len(selfie.n)):
          Dphi.append(np.array([]))
          Z.append(np.array([]))
          n.append(np.array([]))
          Z[tmprind]    = selfie.Z[tmprind] 
          Dphi[tmprind] =selfie.Dphi[tmprind]*signchange #swap coordinate direction
          n[tmprind]    =selfie.n[tmprind]*signchange #swap coordinate direction
              
          if selfie.StelSym: #for m=0, I choose to output only positive n in the bc file
              indsm0=np.where(np.logical_and(selfie.m[tmprind]==0,selfie.n[tmprind]!=0))[0]
              signn=np.sign(n[tmprind][indsm0])
              n[tmprind][indsm0]=signn*n[tmprind][indsm0]
              Dphi[tmprind][indsm0]=signn*Dphi[tmprind][indsm0] #Dphi and Z are sin-series so changing
              Z[tmprind][indsm0]=signn*Z[tmprind][indsm0]       #sign of n -> component sign change
          
          #n[tmprind]=np.where(selfie.m[tmprind]==0,np.abs(selfie.n[tmprind]),-selfie.n[tmprind])

      if minorradiusconvention=='VMEC':
          if np.isnan(selfie.minorradiusVMEC):
              sys.exit("Error: minorradiusconvention='VMEC' could not be used when saving"+
                       "the Boozer file, because the VMEC minor radius is unknown!")
          chosenminorradius=selfie.minorradiusVMEC
      else: #minorradiusconvention=='W7AS' (default)
          chosenminorradius=selfie.minorradiusW7AS
          
      if majorradiusconvention=='VMEC':
          if np.isnan(selfie.majorradiusVMEC):
              sys.exit("Error: majorradiusconvention='VMEC' could not be used when saving"+
                       "the Boozer file, because the VMEC major radius is unknown!")
          chosenmajorradius=selfie.majorradiusVMEC
      else: #majorradiusconvention=='lastR00' (default)
          chosenmajorradius=selfie.majorradiusLastbcR00
              
      #################################################################################
      if filetype=='bc':
      #################################################################################
         if printheadercomment:
            f.write('CC Saved %4d.%2d.%2d %2d:%2d by HS python routine bcgeom.write.\n' %
                    (now.year,now.month,now.day,now.hour,now.minute))
            f.write('CC ---------------------------------------------------------\n')
            f.write(selfie.headertext.maincomment+'\n')

         if selfie.StelSym: #JG includes a little more info in this case
            if (not(np.isnan(selfie.minorradiusVMEC)) and np.isnan(selfie.majorradiusVMEC)
                and np.isnan(selfie.volumeVMEC)):
                f.write(' m0b  n0b nsurf nper  flux/[Tm^2]     a/[m]     R/[m]   avol/[m]\n')
                f.write('%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f\n' %
                        (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                         chosenminorradius,chosenmajorradius,selfie.minorradiusVMEC))
            elif (not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC))
                  and np.isnan(selfie.volumeVMEC)):
                f.write(' m0b  n0b nsurf nper  flux/[Tm^2]     a/[m]     R/[m]   avol/[m]  '+
                        'Rvol/[m]\n')
                f.write('%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n' %
                        (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                         chosenminorradius,chosenmajorradius,
                         selfie.minorradiusVMEC,selfie.majorradiusVMEC,
                         np.pi*selfie.minorradiusVMEC**2.0*2.0*np.pi*selfie.majorradiusVMEC))
            elif (not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC))
                  and not(np.isnan(selfie.volumeVMEC))):
                f.write(' m0b  n0b nsurf nper  flux/[Tm^2]     a/[m]     R/[m]   avol/[m]  '+
                        'Rvol/[m] Vol/[m^3]\n')
                f.write('%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n' %
                        (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                         chosenminorradius,chosenmajorradius,
                         selfie.minorradiusVMEC,selfie.majorradiusVMEC,selfie.volumeVMEC))
            else: #all three are nan
                f.write(' m0b  n0b nsurf nper flux/[Tm^2]     a/[m]     R/[m]\n')
                f.write('%3d%5d%6d%4d%16.6E%10.5f%10.5f\n' %
                        (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                         chosenminorradius,chosenmajorradius))
         else: #non-StelSym, ES uses a little less info than JG
            f.write(' m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]\n')
            #Difficult! I Think Erika uses VMEC quantities here but I am not sure!
            if not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC)):
                f.write('%4d%6d%7d%6d%15.8f%15.8f%15.8f\n' %
                        (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                         selfie.minorradiusVMEC,selfie.majorradiusVMEC))
            else:
                sys.exit('Did not find the minorradiusVMEC and majorradiusVMEC!')

         for rind in range(0,selfie.nsurf):
           if selfie.StelSym:
               f.write('       s         iota  curr_pol/nper    curr_tor    pprime   sqrt '+
                        'g(0,0)\n')
               f.write('                            [A]            [A]   dp/ds,[Pa] (dV/ds)/nper\n')
               f.write('%12.4E%12.4E%12.4E%12.4E%12.4E%12.4E\n' %
                       (selfie.s[rind],iota[rind],
                        Bphi[rind]/(selfie.Nperiods/2.0/np.pi*(4.0*np.pi*1.0e-7)),
                        Btheta[rind]*1.0e7/(1.0/2.0/np.pi*(4.0*np.pi)),
                        selfie.dpds[rind],dVdsoverNper[rind]))
           else:
               f.write('        s               iota           Jpol/nper          Itor'+
                       '            pprime         sqrt g(0,0)\n')
               f.write('                                          [A]           [A] '+
                       '            [Pa]         (dV/ds)/nper\n')
               f.write('%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E\n' %
                       (selfie.s[rind],iota[rind],
                        Bphi[rind]/(selfie.Nperiods/2.0/np.pi*(4.0*np.pi*1.0e-7)),
                        Btheta[rind]*1.0e7/(1.0/2.0/np.pi*(4.0*np.pi)),
                        selfie.dpds[rind],dVdsoverNper[rind]))
         
           if selfie.StelSym:
               f.write('    m    n        r/[m]           z/[m] (phib-phi)*nper/twopi     '+
                       'bmn/[T]\n')
           else:
               f.write('    m    n      rmnc [m]         rmns [m]         zmnc [m]         '+
                        'zmns [m]         vmnc [ ]         vmns [ ]         '+
                        'bmnc [T]         bmns [T]\n')

           if selfie.StelSym:
              if not(nsortstyle is None): #Dont mn-sort the components because they usually come sorted
                  if nsortstyle=='ascend':
                      inds=np.lexsort((n[rind],m[rind]))
                  elif nsortstyle=='descend':
                      inds=np.lexsort((-n[rind],m[rind]))
                  for j in range(0,selfie.nmodes[rind]):
                      rmnc=selfie.R[rind][inds[j]]
                      zmns=Z[rind][inds[j]]
                      vmns=Dphi[rind][inds[j]]
                      bmnc=selfie.B[rind][inds[j]]
                      f.write('%5d%5d%16.8E%16.8E%16.8E%16.8E\n' %
                              (m[rind][inds[j]],n[rind][inds[j]],
                               selfie.R[rind][inds[j]],Z[rind][inds[j]],
                               Dphi[rind][inds[j]],selfie.B[rind][inds[j]]))
              else:
                   for ind in range(0,selfie.nmodes[rind]):
                       f.write('%5d%5d%16.8E%16.8E%16.8E%16.8E\n' %
                               (m[rind][ind],n[rind][ind],
                                selfie.R[rind][ind],Z[rind][ind],
                                Dphi[rind][ind],selfie.B[rind][ind]))

           else:
              #imn=np.concatenate((range(selfie.nmodes[rind]),selfie.m[rind],selfie.n[rind]),axis=1).T
              #imnsort=np.sort(imn,kind='mergesort',axis=[1,2])
              if nsortstyle=='ascend':
                  inds=np.lexsort((selfie.parity[rind],n[rind],m[rind]))
              elif nsortstyle=='descend':
                  inds=np.lexsort((selfie.parity[rind],-n[rind],m[rind]))
              else: #(nsortstyle=='HM')
                  inds_acsend=np.lexsort((selfie.parity[rind],n[rind],m[rind]))
                  indsdescend=np.lexsort((selfie.parity[rind],-n[rind],m[rind]))
                  inds=np.where(m[rind][inds_acsend]==0,inds_acsend,indsdescend)
              j=0
              while j<selfie.nmodes[rind]:
                  thism=m[rind][inds[j]]
                  thisn=n[rind][inds[j]]
                  rmnc=0.0
                  rmns=0.0
                  zmnc=0.0
                  zmns=0.0
                  vmnc=0.0
                  vmns=0.0
                  bmnc=0.0
                  bmns=0.0
                  #print(thism, thisn, j, selfie.nmodes[rind], selfie.parity[rind][inds[j]])
                  if selfie.parity[rind][inds[j]]==1: #j is cosinus components
                       rmnc=selfie.R[rind][inds[j]]
                       zmns=Z[rind][inds[j]]
                       vmns=Dphi[rind][inds[j]]
                       bmnc=selfie.B[rind][inds[j]]
                       j+=1              
                  else: #j is sinus components
                      rmns=selfie.R[rind][inds[j]]
                      zmnc=Z[rind][inds[j]]
                      vmnc=Dphi[rind][inds[j]]
                      bmns=selfie.B[rind][inds[j]]
                      if j+1==selfie.nmodes[rind]:
                          #print('A '+str(j))
                          j+=1
                      else: 
                          if thism==m[rind][inds[j+1]] and thisn==n[rind][inds[j+1]] and selfie.parity[rind][inds[j+1]]==1:
                              #print('B1 '+str(j))
                              rmnc=selfie.R[rind][inds[j+1]]
                              zmns=Z[rind][inds[j+1]]
                              vmns=Dphi[rind][inds[j+1]]
                              bmnc=selfie.B[rind][inds[j+1]]
                              j+=2
                          else:
                              #print('B2 '+str(j))
                              j+=1 
                  f.write('%5d%5d%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E\n' %
                          (thism,thisn,rmnc,rmns,zmnc,zmns,vmnc,vmns,bmnc,bmns))

      #################################################################################
      else: #save a .dat file instead
      #################################################################################
         if majorradiusconvention=='VMEC':
             print('Warning: You have chosen to make a Henning Maassberg style .dat file '+
                   'with the major radius taken from VMEC instead of using the standard "lastR00" definition.')
         else:
             if np.isnan(selfie.majorradiusLastbcR00):
                 chosenmajorradius=0.0


         if minorradiusconvention=='VMEC':
             print('Warning: You have chosen to make a Henning Maassberg style .dat file '+
                   'with the minor radius taken from VMEC instead of using the W7AS standard definition.')
             
         if printheadercomment:
             f.write('cc Saved %4d.%2d.%2d %2d:%2d by HS python routine bcgeom.write.\n' %
                     (now.year,now.month,now.day,now.hour,now.minute))
             f.write('cc ---------------------------------------------------------\n')
             f.write(selfie.headertext.maincomment+'\n')
             f.write('cc  m0b  n0b nsurf nper  flux/[Tm^2]     a/[m]     '+
                     'R/[m]   avol/[m]  Rvol/[m] Vol/[m^3]\n')

             if not(np.isnan(selfie.m0b)) and not(np.isnan(selfie.n0b)) and not(np.isnan(selfie.torfluxtot)):
                if (not(np.isnan(selfie.minorradiusVMEC)) and np.isnan(selfie.majorradiusVMEC)
                    and np.isnan(selfie.volumeVMEC)):
                   f.write('cc %3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f\n'%
                           (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                            chosenminorradius,chosenmajorradius,selfie.minorradiusVMEC))
                elif (not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC))
                      and np.isnan(selfie.volumeVMEC)):
                   f.write('cc %3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n'%
                           (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                            chosenminorradius,chosenmajorradius,
                            selfie.minorradiusVMEC,selfie.majorradiusVMEC,
                            np.pi*selfie.minorradiusVMEC**2.0*2.0*np.pi*selfie.majorradiusVMEC))
                elif (not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC))
                      and not(np.isnan(selfie.volumeVMEC))):
                   f.write('cc %3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n'%
                           (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                            chosenminorradius,chosenmajorradius,
                            selfie.minorradiusVMEC,selfie.majorradiusVMEC,selfie.volumeVMEC))
                else:
                   f.write('cc %3d%5d%6d%4d%16.6E%10.5f%10.5f\n'%
                           (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                            chosenminorradius,chosenmajorradius))
             else:
                if (not(np.isnan(selfie.minorradiusVMEC)) and np.isnan(selfie.majorradiusVMEC)
                    and np.isnan(selfie.volumeVMEC)):
                   f.write('cc %6d%4d%10.5f%10.5f%11.5f\n'%
                           (selfie.nsurf,selfie.Nperiods,
                            chosenminorradius,chosenmajorradius,selfie.minorradiusVMEC))
                elif (not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC))
                      and np.isnan(selfie.volumeVMEC)):
                   f.write('cc %6d%4d%10.5f%10.5f%11.5f%11.5f%11.5f\n'%
                           (selfie.nsurf,selfie.Nperiods,
                            chosenminorradius,chosenmajorradius,
                            selfie.minorradiusVMEC,selfie.majorradiusVMEC,
                            np.pi*selfie.minorradiusVMEC**2.0*2.0*np.pi*selfie.majorradiusVMEC))
                elif (not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC))
                      and not(np.isnan(selfie.volumeVMEC))):
                   f.write(' cc %6d%4d%10.5f%10.5f%11.5f%11.5f%11.5f\n'%
                           (selfie.nsurf,selfie.Nperiods,
                            chosenminorradius,chosenmajorradius,
                            selfie.minorradiusVMEC,selfie.majorradiusVMEC,selfie.volumeVMEC))
                else:
                   f.write('cc %6d%4d%10.5f%10.5f\n'%
                           (selfie.nsurf,selfie.Nperiods,
                            chosenminorradius,chosenmajorradius))
         #end if headercomment
         
         f.write(' c     m, n, B_mn/B_00, B_mn, (dB_mn/dr)/B_00 [1/cm], R_mn, z_mn\n')

         for rind in range(selfie.nsurf):
             f.write('%6.2f%8.4f%3i%8.2f%8.2f%7.4f%8.2f%8.2f  r,io,N,ra,Rm,s,r_o,R00\n'%
                     (chosenminorradius*selfie.rnorm[rind]*100.0,
                      iota[rind],selfie.Nperiods,chosenminorradius*100.0,
                      chosenmajorradius*100.0,selfie.s[rind],
                      chosenminorradius*selfie.rnorm[rind]*100.0,selfie.R00[rind]*100.0))
             f.write(' >  %13.4E%13.4E%13.4E     cur_pol=g, cur_tor=I, B_00\n'%
                     (Bphi[rind]/(1.0e6*selfie.Nperiods/2.0/np.pi*(4.0*np.pi*1.0e-7)),
                      Btheta[rind]/(1.0e6/2.0/np.pi*(4.0*np.pi*1.0e-7)),
                      selfie.B00[rind]))
             m0inds=np.where(m[rind]==0)[0]
             if nsortstyle=='ascend' or nsortstyle=='HM':
                 m0inds=m0inds[np.argsort(n[rind][m0inds])]
             else:
                 m0inds=m0inds[np.argsort(-n[rind][m0inds])]   
             for m0i in m0inds:
                 f.write('%4i%4i%11.6f%11.6f%12.3E%11.5f%11.5f\n'%
                         (m[rind][m0i],
                          n[rind][m0i],
                          selfie.Bnorm[rind][m0i],
                          selfie.B[rind][m0i],
                          selfie.dBdrnorm[rind][m0i]/selfie.B00[rind]/chosenminorradius/100.0,
                          selfie.R[rind][m0i]*100.0,
                          Z[rind][m0i]*100.0))
             for thism in range(1,int(np.max(m[rind]))+1):
                 thisminds=np.where(m[rind]==thism)[0]
                 if nsortstyle=='ascend':
                     thisminds=thisminds[np.argsort(n[rind][thisminds])]
                 else: #(including nsortstyle=='HM')
                     thisminds=thisminds[np.argsort(-n[rind][thisminds])]
                 for mi in thisminds:
                     f.write('%4i%4i%11.6f%11.6f%12.3E%11.5f%11.5f\n'%
                             (m[rind][mi],
                              n[rind][mi],
                              selfie.Bnorm[rind][mi],
                              selfie.B[rind][mi],
                              selfie.dBdrnorm[rind][mi]/selfie.B00[rind]/chosenminorradius/100.0,
                              selfie.R[rind][mi]*100.0,
                              Z[rind][mi]*100.0))
             #end for thism
             f.write('  -1   0 0.0 0.0 0.0 0.0 0.0     end of input\n')
         #end for rind
         
      #end if filetype=='bc' or 'dat'   
      f.close()

  ####################################################################
  # Interpolate a bcgeom to a set of new radii
  ####################################################################
  def interp(self,invar,invarname='rnorm',interptype='linear',padding=None):
      # This function interpolates the bcgeom to a new
      # set of radii. The chosen new values invar of the flux label 
      # given by the string invarname ('s', 'rnorm', 'rW7AS', 'rVMEC') are interpolated.
      # Default is 'rnorm'.
      # The input invarname can also be 'index', in which case no interpolation
      # takes place, the output is just self at the indices given in invar.
      #
      # The interptype can be 'linear' or 'nearest', default is 'linear'.


      if np.any(np.isnan(invar)):
          sys.exit('Input invar to bcgeom.interp contains nans!')
        
      smin=self.s[0]
      smax=self.s[-1]
      
      #Decide which variable to interpolate over
      if invarname=='s':
          s=invar
          if (padding is None) and (np.any(s<smin) or np.any(s>smax)):
            print('In bcgeom.interp: s='+str(s))
            print('smin='+str(smin)+', smax='+str(smax))
            sys.exit('Chosen s out of range!')
      elif invarname=='rnorm':
          s=invar**2
          if (padding is None) and (np.any(s<smin) or np.any(s>smax)):
            sys.exit('Chosen rnorm '+str(invar)+' out of range! rnorm_min='+str(np.sqrt(smin))+
                     ', rnorm_max='+str(np.sqrt(smax)))
      elif invarname=='rW7AS':
          s=(invar/self.minorradiusW7AS)**2
          if (padding is None) and (np.any(s<smin) or np.any(s>smax)):
            sys.exit('Chosen rW7AS out of range!')
      elif invarname=='rVMEC':
          s=(invar/self.minorradiusVMEC)**2
          if (padding is None) and (np.any(s<smin) or np.any(s>smax)):
            sys.exit('Chosen rVMEC out of range!')
      elif invarname=='index':
          if np.any(invar<0) or np.any(invar>(len(self.s)-1)):
            sys.exit('Chosen index out of range!')
          invar=np.sort(invar)
          s=self.s(invar)
          interptype='nearest'
      else:
          sys.exit('invarname not recognised!')
      if padding=='end values' and not(invarname=='index'):
          s=np.where(s<smin,smin,s)
          s=np.where(s>smax,smax,s)
          
      out=copy.deepcopy(self)

      s=np.sort(s)
      if np.any(np.diff(s)<0):
          print('s='+str(s))
          sys.exit('In bcgeom.interp: '+invarname+' is non-monotonic')
      if np.any(np.diff(s)==0):
          print('s='+str(s))
          sys.exit('In bcgeom.interp: '+invarname+' contains doublets')
      out.nsurf=len(s)

      if interptype=='linear': #This is the default case
          actualinterpolation=True
          if len(s)==len(self.s):
              actualinterpolation=np.any(s!=self.s)
          if actualinterpolation:
             out.headertext.maincomment=(out.headertext.maincomment+'\n'+
                    'CC The file has been linearly interpolated to a new set of radii.')
          
          out.s=s
          out.rnorm        = np.sqrt(s)
          out.iota         = np.interp(s,self.s,self.iota)
          out.Bphi         = np.interp(s,self.s,self.Bphi)
          out.Btheta       = np.interp(s,self.s,self.Btheta)
          out.dpds         = np.interp(s,self.s,self.dpds)
          out.dVdsoverNper = np.interp(s,self.s,self.dVdsoverNper)
          out.B00          = np.interp(s,self.s,self.B00)
          out.R00          = np.interp(s,self.s,self.R00)

          out.m=[]
          out.n=[]
          if not(self.StelSym):
             out.parity=[]
          out.B=[]
          out.dBdrnorm=[]
          out.R=[]
          out.Z=[]
          out.Bnorm=[]
          out.Dphi=[]

          for surfind in range(out.nsurf):
              #print('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%4i/%4i'% (surfind+1,out.nsurf))
              if not(s[surfind] in self.s):
                  Rind=np.where(np.sign(self.s-s[surfind])==1)[0][0]
              else:
                  exactind=np.where(self.s==s[surfind])[0][0]
                  Rind=exactind+1
                  if exactind==len(s)-1:
                      Rind=len(s)-1

              #print('self.s='+str(self.s))
              #print('s='+str(s[surfind]))
              #print('Rind='+str(Rind))
              Lind=Rind-1
              
              Ds=self.s[Rind]-self.s[Lind]
              Drnorm=self.rnorm[Rind]-self.rnorm[Lind]
              wR=(s[surfind]-self.s[Lind])/Ds
              wL=(self.s[Rind]-s[surfind])/Ds
              
              out.m.append(np.array([]))
              out.n.append(np.array([]))
              if not(self.StelSym):
                  out.parity.append(np.array([]))
              out.B.append(np.array([]))
              out.dBdrnorm.append(np.array([]))
              out.Bnorm.append(np.array([]))
              out.R.append(np.array([]))
              out.Z.append(np.array([]))    
              out.Dphi.append(np.array([])) 


              if not(self.StelSym):
                 mode=-1
                 for par in range(2):
                    pindsL=np.where(self.parity[Lind]==par)[0]
                    pindsR=np.where(self.parity[Rind]==par)[0]

                    mL=self.m[Lind][pindsL]
                    nL=self.n[Lind][pindsL]
                    mR=self.m[Rind][pindsR]
                    nR=self.n[Rind][pindsR]

                    allns=np.sort(np.concatenate((nL,nR)))
                    if len(allns)>1:    
                        ns=allns[np.append(np.nonzero(np.diff(allns))[0],len(allns)-1)]
                    else:
                        ns=allns
                        
                    Nn=len(ns)

                    for nind in range(Nn):
                       n=ns[ninds]
                       allms=np.sort(np.concatenate((mL[np.where(nL==n)[0]],
                                                     mR[np.where(nR==n)[0]])))
                       if len(allms)>1:
                           ms=allms[np.append(np.nonzero(np.diff(allms))[0],len(allms)-1)]
                       else:
                           ms=allms
                       Nm=len(ms)

                       for mind in range(Nm):
                          mode+=1
                          m=ms[mind]
                          Lmnind=np.where(np.logical_and(self.m[Lind]==m,
                                                         self.n[Lind]==n,
                                                         self.parity[Lind]==par))[0]
                          Lexists=(Lmnind.size!=0)
                          Rmnind=np.where(np.logical_and(self.m[Rind]==m,
                                                         self.n[Rind]==n,
                                                         self.parity[Rind]==par))[0]
                          Rexists=(Rmnind.size!=0)
                          out.m[surfind]=np.append(out.m[surfind],m)
                          out.n[surfind]=np.append(out.n[surfind],n)
                          out.parity[surfind]=np.append(out.parity[surfind],par)
                          if Lexists and Rexists:
                             out.B[surfind]=np.append(out.B[surfind],
                                                      wL*self.B[Lind][Lmnind]+
                                                      wR*self.B[Rind][Rmnind])
                             out.dBdrnorm[surfind]=np.append(out.dBdrnorm[surfind],
                                        (self.B[Rind][Rmnind]-self.B[Lind][Lmnind])/Drnorm)
                             out.R[surfind]=np.append(out.R[surfind],
                                                      wL*self.R[Lind][Lmnind]+
                                                      wR*self.R[Rind][Rmnind])
                             out.Z[surfind]=np.append(out.Z[surfind],
                                                      wL*self.Z[Lind][Lmnind]+
                                                      wR*self.Z[Rind][Rmnind])
                             out.Dphi[surfind]=np.append(out.Dphi[surfind],
                                                      wL*self.Dphi[Lind][Lmnind]+
                                                      wR*self.Dphi[Rind][Rmnind])
                          elif Lexists and not(Rexists):
                             out.B[surfind]=np.append(out.B[surfind],
                                                      wL*self.B[Lind][Lmnind])
                             out.dBdrnorm[surfind]=np.append(out.dBdrnorm[surfind],
                                                         -self.B[Lind][Lmnind]/Drnorm)
                             out.R[surfind]=np.append(out.R[surfind],
                                                      wL*self.R[Lind][Lmnind])
                             out.Z[surfind]=np.append(out.Z[surfind],
                                                      wL*self.Z[Lind][Lmnind])
                             out.Dphi[surfind]=np.append(out.Dphi[surfind],
                                                      wL*self.Dphi[Lind][Lmnind])
                          elif not(Lexists) and Rexists:
                             out.B[surfind]=np.append(out.B[surfind],
                                                      wR*self.B[Rind][Rmnind])
                             out.dBdrnorm[surfind]=np.append(out.dBdrnorm[surfind],
                                                         self.B[Rind][Rmnind]/Drnorm)
                             out.R[surfind]=np.append(out.R[surfind],
                                                      wR*self.R[Rind][Rmnind])
                             out.Z[surfind]=np.append(out.Z[surfind],
                                                      wR*self.Z[Rind][Rmnind])
                             out.Dphi[surfind]=np.append(out.Dphi[surfind],
                                                      wR*self.Dphi[Rind][Rmnind])
                          else:
                             sys.exit('This is impossible')

              else: #Stellarator symmetric case
                 mode=-1

                 mL=self.m[Lind]
                 nL=self.n[Lind]
                 mR=self.m[Rind]
                 nR=self.n[Rind]

                 allns=np.sort(np.concatenate((nL,nR)))
                 if len(allns)>1:
                     ns=allns[np.append(np.nonzero(np.diff(allns))[0],len(allns)-1)]
                 else:
                     ns=allns
                     
                 Nn=len(ns)

                 for nind in range(Nn):
                    n=ns[nind]
                    allms=np.sort(np.concatenate((mL[np.where(nL==n)[0]],
                                                  mR[np.where(nR==n)[0]])))
                    if len(allms)>1:
                        ms=allms[np.append(np.nonzero(np.diff(allms))[0],len(allms)-1)]
                    else:
                        ms=allms
                        
                    Nm=len(ms)

                    for mind in range(Nm):
                       mode+=1
                       m=ms[mind]
                       Lmnind=np.where(np.logical_and(self.m[Lind]==m,
                                                      self.n[Lind]==n))[0]
                       Lexists=(Lmnind.size!=0)
                       Rmnind=np.where(np.logical_and(self.m[Rind]==m,
                                                      self.n[Rind]==n))[0]
                       Rexists=(Rmnind.size!=0)
                       out.m[surfind]=np.append(out.m[surfind],m)
                       out.n[surfind]=np.append(out.n[surfind],n)
                       if Lexists and Rexists:
                          out.B[surfind]=np.append(out.B[surfind],
                                                   wL*self.B[Lind][Lmnind]+
                                                   wR*self.B[Rind][Rmnind])
                          out.dBdrnorm[surfind]=np.append(out.dBdrnorm[surfind],
                                        (self.B[Rind][Rmnind]-self.B[Lind][Lmnind])/Drnorm)
                          out.R[surfind]=np.append(out.R[surfind],
                                                   wL*self.R[Lind][Lmnind]+
                                                   wR*self.R[Rind][Rmnind])
                          out.Z[surfind]=np.append(out.Z[surfind],
                                                   wL*self.Z[Lind][Lmnind]+
                                                   wR*self.Z[Rind][Rmnind])
                          out.Dphi[surfind]=np.append(out.Dphi[surfind],
                                                   wL*self.Dphi[Lind][Lmnind]+
                                                   wR*self.Dphi[Rind][Rmnind])
                       elif Lexists and not(Rexists):
                          out.B[surfind]=np.append(out.B[surfind],
                                                   wL*self.B[Lind][Lmnind])
                          out.dBdrnorm[surfind]=np.append(out.dBdrnorm[surfind],
                                                      -self.B[Lind][Lmnind]/Drnorm)
                          out.R[surfind]=np.append(out.R[surfind],
                                                   wL*self.R[Lind][Lmnind])
                          out.Z[surfind]=np.append(out.Z[surfind],
                                                   wL*self.Z[Lind][Lmnind])
                          out.Dphi[surfind]=np.append(out.Dphi[surfind],
                                                   wL*self.Dphi[Lind][Lmnind])
                       elif not(Lexists) and Rexists:
                          out.B[surfind]=np.append(out.B[surfind],
                                                   wR*self.B[Rind][Rmnind])
                          out.dBdrnorm[surfind]=np.append(out.dBdrnorm[surfind],
                                                      self.B[Rind][Rmnind]/Drnorm)
                          out.R[surfind]=np.append(out.R[surfind],
                                                   wR*self.R[Rind][Rmnind])
                          out.Z[surfind]=np.append(out.Z[surfind],
                                                   wR*self.Z[Rind][Rmnind])
                          out.Dphi[surfind]=np.append(out.Dphi[surfind],
                                                   wR*self.Dphi[Rind][Rmnind])
                       else:
                          sys.exit('This is impossible')
                    #end for mind
                 #end for nind
              #end if not(StelSym)
              
              out.Bnorm[surfind] = out.B[surfind]/out.B00[surfind]
              out.nmodes[surfind]= mode+1
          #end for surfind
          #print('\n')
      elif interptype=='nearest': #Pick the nearest existing surface
        actualinterpolation=True
        if len(s)==len(self.s):
            actualinterpolation=np.any(s!=self.s)
        if actualinterpolation:
            out.headertext.maincomment=(out.headertext.maincomment+'\n'+
                    'CC The file has only a subset of the flux surfaces of the original file.')

        if invarname=='index':
           inds=invar
        else:
           inds=np.zeros((len(s),),dtype=np.int)
           if invarname=='s':
              for iind in range(len(s)):
                 inds[iind]=(np.abs(self.s-s[iind])).argmin()
           else:
              for iind in range(len(s)):
                 inds[iind]=(np.abs(np.sqrt(self.s)-np.sqrt(s[iind]))).argmin()

        out.s      = self.s[inds]
        out.rnorm  = np.sqrt(out.s)
        out.iota   = self.iota[inds]
        out.Bphi   = self.Bphi[inds]
        out.Btheta = self.Btheta[inds]
        out.dpds   = self.dpds[inds]
        out.dVdsoverNper = self.dVdsoverNper[inds]
        out.B00    = self.B00[inds]
        out.R00    = self.R00[inds]
        out.nmodes = np.array(self.nmodes)[inds]
        out.B      = np.array(self.B)[inds]
        out.R      = np.array(self.R)[inds]
        out.Z      = np.array(self.Z)[inds]
        out.Dphi   = np.array(self.Dphi)[inds]
        out.Bnorm  = np.array(self.Bnorm)[inds]
        out.m      = np.array(self.m)[inds]
        out.n      = np.array(self.n)[inds]
        if not(self.StelSym):
           out.parity = self.parity[inds]

      else:
        sys.exit('Unknown interpolation type!')

      return out

        
  ####################################################################
  # Display the content of a bcgeom
  ####################################################################
  
  def disp(self,verbose=0):

    print('-----------------------------------------------')
    print('Comment lines: headertext.maincomment=')
    #print('-----------------------------------------------')

    print(self.headertext.maincomment)
    print('-----------------------------------------------')
    print('Scalars:')
    #print('-----------------------------------------------')
    print('nsurf                = '+str(self.nsurf))
    print('Nperiods             = '+str(self.Nperiods))
    print('psi_a                = '+str(self.psi_a))
    print('torfluxtot           = '+str(self.torfluxtot))
    print('minorradiusW7AS      = '+str(self.minorradiusW7AS))
    print('majorradiusLastbcR00 = '+str(self.majorradiusLastbcR00))
    print('minorradiusVMEC      = '+str(self.minorradiusVMEC))
    print('majorradiusVMEC      = '+str(self.majorradiusVMEC))
    print('volumeVMEC           = '+str(self.volumeVMEC))
    print('Bfilter.min_Bmn      = '+str(self.Bfilter.min_Bmn))
    print('Bfilter.max_m        = '+str(self.Bfilter.max_m))
    print('Bfilter.maxabs_n     = '+str(self.Bfilter.maxabs_n))

    print('-----------------------------------------------')
    print('Radial arrays:')
    #print('-----------------------------------------------')
    if verbose==1:
      print('rnorm='+str(self.rnorm))
      print('s='+str(self.s))
      print('Bphi='+str(self.Bphi))
      print('Btheta='+str(self.Btheta))
      print('iota='+str(self.iota))
      print('dpds='+str(self.dpds))
  
      print('dVdsoverNper='+str(self.dVdsoverNper))
      print('FSAB2='+str(self.FSAB2))
      print('nmodes='+str(self.nmodes))
      print('B00='+str(self.B00))
      print('R00='+str(self.R00))
    elif verbose==0:
      print('rnorm, s, Bphi, Btheta, iota, dpds, dVdsoverNper, ')
      print('FSAB2, nmodes, B00, R00')

    print('-----------------------------------------------')
    if self.StelSym:
       print('Fourier indexes:    m, n')
    else:
       print('Fourier indexes:    m, n, parity')

    print('Fourier quantities: B, R, Z, Dphi')
    print('-----------------------------------------------')



##########################################################################################
##########################################################################################
######################################## vmecgeom ########################################
##########################################################################################
##########################################################################################
    

class vmecgeom(object):
    
  def __init__(self,filename,min_Bmn=0,max_m=np.inf,maxabs_n=np.inf,
               symmetry='unknown'):
    # This reads the VMEC file of the wout type
    # Note that the wout file is left-handed and the output struct self
    # is right-handed.
    #
    # Only Fourier components with abs(Bmn)>min_Bmn, m<=max_m, |n|<=maxabs_n
    # are read.
    #
    # The input symmetry can be 'StelSym' to double-check that the 
    # input is stellarator symmetric
    #


    def ncdump(nc_fid, verb=False):
        '''
        ncdump outputs dimensions, variables and their attribute information.
        The information is similar to that of NCAR's ncdump utility.
        ncdump requires a valid instance of Dataset.

        Parameters
        ----------
        nc_fid : netCDF4.Dataset
            A netCDF4 dateset object
        verb : Boolean
            whether or not nc_attrs, nc_dims, and nc_vars are printed

        Returns
        -------
        nc_attrs : list
            A Python list of the NetCDF file global attributes
        nc_dims : list
            A Python list of the NetCDF file dimensions
        nc_vars : list
            A Python list of the NetCDF file variables
        '''
        def print_ncattr(key):
            """
            Prints the NetCDF file attributes for a given key

            Parameters
            ----------
            key : unicode
                a valid netCDF4.Dataset.variables key
            """
            try:
                print("\t\ttype:", repr(nc_fid.variables[key].dtype))
                for ncattr in nc_fid.variables[key].ncattrs():
                    print('\t\t%s:' % ncattr,\
                          repr(nc_fid.variables[key].getncattr(ncattr)))
            except KeyError:
                print("\t\tWARNING: %s does not contain variable attributes" % key)

        # NetCDF global attributes
        nc_attrs = nc_fid.ncattrs()
        if verb:
            print("NetCDF Global Attributes:")
            for nc_attr in nc_attrs:
                print('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
        nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
        # Dimension shape information.
        if verb:
            print("NetCDF dimension information:")
            for dim in nc_dims:
                print("\tName:", dim )
                print("\t\tsize:", len(nc_fid.dimensions[dim]))
                print_ncattr(dim)
        # Variable information.
        nc_vars = [var for var in nc_fid.variables]  # list of nc variables
        if verb:
            print("NetCDF variable information:")
            for var in nc_vars:
                if var not in nc_dims:
                    print('\tName:', var)
                    print("\t\tdimensions:", nc_fid.variables[var].dimensions)
                    print("\t\tsize:", nc_fid.variables[var].size)
                    print_ncattr(var)
        return nc_attrs, nc_dims, nc_vars


    if filename.endswith('.nc'):
      dataset=Dataset(filename)
      nc_attrs, nc_dims, nc_vars = ncdump(dataset)
      #print(nc_vars
    
      self.StelSym=not('bmns' in dataset.variables)
      self.skip = 1 #=1,this is how many elements are skipped at low radii when going to half grid
      self.filename         =filename
      self.version_         =float(dataset.variables['version_'][:])
      self.input_extension  =''.join(dataset.variables['input_extension'][:])
      self.mgrid_file       =''.join(dataset.variables['mgrid_file'][:])
      if 'pcurr_type'.encode('utf-8') in nc_vars:
          self.pcurr_type       =''.join(dataset.variables['pcurr_type'][:])
      if 'pmass_type'.encode('utf-8') in nc_vars:
          self.pmass_type       =''.join(dataset.variables['pmass_type'][:])
      if 'piota_type'.encode('utf-8') in nc_vars:
          self.piota_type       =''.join(dataset.variables['piota_type'][:])
      self.wb               =float(dataset.variables['wb'][:])
      self.wp               =float(dataset.variables['wp'][:])
      self.gamma            =float(dataset.variables['gamma'][:])
      self.rmax_surf        =float(dataset.variables['rmax_surf'][:])
      self.rmin_surf        =float(dataset.variables['rmin_surf'][:])
      self.zmax_surf        =float(dataset.variables['zmax_surf'][:])
      self.nfp              =int(dataset.variables['nfp'][:])
      self.ns               =int(dataset.variables['ns'][:])
      self.mpol             =int(dataset.variables['mpol'][:])
      self.ntor             =int(dataset.variables['ntor'][:])
      self.mnmax            =int(dataset.variables['mnmax'][:])
      self.mnmax_nyq        =int(dataset.variables['mnmax_nyq'][:])
      self.niter            =int(dataset.variables['niter'][:])
      self.itfsq            =float(dataset.variables['itfsq'][:])
      self.lasym__logical__ =bool(dataset.variables['lasym__logical__'][:])
      self.lrecon__logical__=bool(dataset.variables['lrecon__logical__'][:])
      self.lfreeb__logical__=bool(dataset.variables['lfreeb__logical__'][:])
      if 'lrfp__logical__'.encode('utf-8') in nc_vars:
          self.lrfp__logical__  =bool(dataset.variables['lrfp__logical__'][:])
      self.ier_flag         =int(dataset.variables['ier_flag'][:])
      self.aspect           =float(dataset.variables['aspect'][:])
      self.betatotal        =float(dataset.variables['betatotal'][:])
      self.betapol          =float(dataset.variables['betapol'][:])
      self.betator          =float(dataset.variables['betator'][:])
      self.betaxis          =float(dataset.variables['betaxis'][:])
      self.b0               =float(dataset.variables['b0'][:])
      self.rbtor0           =float(dataset.variables['rbtor0'][:])
      self.rbtor            =float(dataset.variables['rbtor'][:])
      self.signgs           =int(dataset.variables['signgs'][:])
      self.IonLarmor        =float(dataset.variables['IonLarmor'][:])
      self.volavgB          =float(dataset.variables['volavgB'][:])
      self.ctor             =float(dataset.variables['ctor'][:])
      self.Aminor_p         =float(dataset.variables['Aminor_p'][:])
      self.Rmajor_p         =float(dataset.variables['Rmajor_p'][:])
      self.volume_p         =float(dataset.variables['volume_p'][:])
      if 'ftolv'.encode('utf-8') in nc_vars:
          self.ftolv            =float(dataset.variables['ftolv'][:])
      if 'fsql'.encode('utf-8') in nc_vars:
          self.fsql             =float(dataset.variables['fsql'][:])
      if 'fsqr'.encode('utf-8') in nc_vars:
          self.fsqr             =float(dataset.variables['fsqr'][:])
      if 'fsqz'.encode('utf-8') in nc_vars:
          self.fsqz             =float(dataset.variables['fsqz'][:])
      self.nextcur          =int(dataset.variables['nextcur'][:])
      self.extcur           =np.array(dataset.variables['extcur'][:]).astype(float)
      self.mgrid_mode       =''.join(dataset.variables['mgrid_mode'][:])
      self.xm               =np.array(dataset.variables['xm'][:]).astype(int)
      self.xn               =np.array(dataset.variables['xn'][:]).astype(int)
      self.xm_nyq           =np.array(dataset.variables['xm_nyq'][:]).astype(int)
      self.xn_nyq           =np.array(dataset.variables['xn_nyq'][:]).astype(int)
      self.raxis_cc         =np.array(dataset.variables['raxis_cc'][:]).astype(float)
      self.zaxis_cs         =np.array(dataset.variables['zaxis_cs'][:]).astype(float)
      if 'raxis_cs' in dataset.variables:
          self.raxis_cs         =np.array(dataset.variables['raxis_cs'][:]).astype(float)
          self.zaxis_cc         =np.array(dataset.variables['zaxis_cc'][:]).astype(float)
      else:
          self.raxis_cs=None
          self.zaxis_cc=None
      if 'am'.encode('utf-8') in nc_vars:
          self.am               =np.array(dataset.variables['am'][:]).astype(float)
      if 'ac'.encode('utf-8') in nc_vars:
          self.ac               =np.array(dataset.variables['ac'][:]).astype(float)
      if 'ai'.encode('utf-8') in nc_vars:
          self.ai               =np.array(dataset.variables['ai'][:]).astype(float)
      if 'am_aux_s'.encode('utf-8') in nc_vars:
          self.am_aux_s         =np.array(dataset.variables['am_aux_s'][:]).astype(float)
      if 'am_aux_f'.encode('utf-8') in nc_vars:
          self.am_aux_f         =np.array(dataset.variables['am_aux_f'][:]).astype(float)
      if 'am_aux_s'.encode('utf-8') in nc_vars:
          self.ai_aux_s         =np.array(dataset.variables['ai_aux_s'][:]).astype(float)
      if 'ai_aux_f'.encode('utf-8') in nc_vars:
          self.ai_aux_f         =np.array(dataset.variables['ai_aux_f'][:]).astype(float)
      if 'ac_aux_s'.encode('utf-8') in nc_vars:
          self.ac_aux_s         =np.array(dataset.variables['ac_aux_s'][:]).astype(float)
      if 'ac_aux_f'.encode('utf-8') in nc_vars:
          self.ac_aux_f         =np.array(dataset.variables['ac_aux_f'][:]).astype(float)
      self.iotaf            =np.array(dataset.variables['iotaf'][:]).astype(float)
      if 'q_factor' in dataset.variables:
        self.q_factor         =np.array(dataset.variables['q_factor'][:]).astype(float)
      else:
        self.q_factor         =np.nan*self.iotaf
      self.presf            =np.array(dataset.variables['presf'][:]).astype(float)
      self.phi              =np.array(dataset.variables['phi'][:]).astype(float)
      self.phipf            =np.array(dataset.variables['phipf'][:]).astype(float)
      if 'chi' in dataset.variables:
        self.chi              =np.array(dataset.variables['chi'][:]).astype(float)
        self.chipf            =np.array(dataset.variables['chipf'][:]).astype(float)
      else:
        self.chi              =np.nan*self.phi
        self.chipf            =np.nan*self.phi
      self.jcuru            =np.array(dataset.variables['jcuru'][:]).astype(float)
      self.jcurv            =np.array(dataset.variables['jcurv'][:]).astype(float)
      self.iotas            =np.array(dataset.variables['iotas'][:]).astype(float)
      self.mass             =np.array(dataset.variables['mass'][:]).astype(float)
      self.pres             =np.array(dataset.variables['pres'][:]).astype(float)
      self.beta_vol         =np.array(dataset.variables['beta_vol'][:]).astype(float)
      self.buco             =np.array(dataset.variables['buco'][:]).astype(float)
      self.bvco             =np.array(dataset.variables['bvco'][:]).astype(float)
      self.vp               =np.array(dataset.variables['vp'][:]).astype(float)
      self.specw            =np.array(dataset.variables['specw'][:]).astype(float)
      self.phips            =np.array(dataset.variables['phips'][:]).astype(float)
      self.over_r           =np.array(dataset.variables['over_r'][:]).astype(float)
      self.jdotb            =np.array(dataset.variables['jdotb'][:]).astype(float)
      self.bdotgradv        =np.array(dataset.variables['bdotgradv'][:]).astype(float)
      self.DMerc            =np.array(dataset.variables['DMerc'][:]).astype(float)
      self.DShear           =np.array(dataset.variables['DShear'][:]).astype(float)
      self.DWell            =np.array(dataset.variables['DWell'][:]).astype(float)
      self.DCurr            =np.array(dataset.variables['DCurr'][:]).astype(float)
      self.DGeod            =np.array(dataset.variables['DGeod'][:]).astype(float)
      self.equif            =np.array(dataset.variables['equif'][:]).astype(float)
      self.fsqt             =np.array(dataset.variables['fsqt'][:]).astype(float)
      self.wdot             =np.array(dataset.variables['wdot'][:]).astype(float)

      self.rmnc             =np.array(dataset.variables['rmnc'][:,:]).astype(float)
      self.zmns             =np.array(dataset.variables['zmns'][:,:]).astype(float)
      self.lmns             =np.array(dataset.variables['lmns'][:,:]).astype(float)
      self.gmnc             =np.array(dataset.variables['gmnc'][:,:]).astype(float)
      self.bmnc             =np.array(dataset.variables['bmnc'][:,:]).astype(float)
      self.bsubumnc         =np.array(dataset.variables['bsubumnc'][:,:]).astype(float)
      self.bsubvmnc         =np.array(dataset.variables['bsubvmnc'][:,:]).astype(float)
      self.bsubsmns         =np.array(dataset.variables['bsubsmns'][:,:]).astype(float)
      self.bsupumnc         =np.array(dataset.variables['bsupumnc'][:,:]).astype(float)
      self.bsupvmnc         =np.array(dataset.variables['bsupvmnc'][:,:]).astype(float)
      if self.StelSym:
        self.rmns             =np.nan
        self.zmnc             =np.nan
        self.lmnc             =np.nan
        self.gmns             =np.nan
        self.bmns             =np.nan
        self.bsubumns         =np.nan
        self.bsubvmns         =np.nan
        self.bsubsmnc         =np.nan
        self.bsupumns         =np.nan
        self.bsupvmns         =np.nan
      else:
        self.rmns             =np.array(dataset.variables['rmns'][:,:]).astype(float)
        self.zmnc             =np.array(dataset.variables['zmnc'][:,:]).astype(float)
        self.lmnc             =np.array(dataset.variables['lmnc'][:,:]).astype(float)
        self.gmns             =np.array(dataset.variables['gmns'][:,:]).astype(float)
        self.bmns             =np.array(dataset.variables['bmns'][:,:]).astype(float)
        self.bsubumns         =np.array(dataset.variables['bsubumns'][:,:]).astype(float)
        self.bsubvmns         =np.array(dataset.variables['bsubvmns'][:,:]).astype(float)
        self.bsubsmnc         =np.array(dataset.variables['bsubsmnc'][:,:]).astype(float)
        self.bsupumns         =np.array(dataset.variables['bsupumns'][:,:]).astype(float)
        self.bsupvmns         =np.array(dataset.variables['bsupvmns'][:,:]).astype(float)

    #######################################################################################    
    # Read VMEC .txt file
    #######################################################################################    
    elif filename.endswith('.txt'):
      self.filename=filename
      self.input_extension='txt'
      def sscan(strng,dtype=float,count=-1):
        if isinstance(strng,list):
            out=[None]*len(strng)
            for i in range(len(strng)):
                out[i]=np.fromstring(strng[i],dtype,count,sep=' ')
            return tuple(out)  
        else:    
            return np.fromstring(strng,dtype,count,sep=' ')
        
      eps_w = 1.0e-4 #a small number.
      with open(filename) as f:
          fl = f.readlines()
      c_vmec_version_string=fl[0]
      c_version_number=c_vmec_version_string[c_vmec_version_string.find('=')+1:]
      r_version_number=float(c_version_number)
      self.version_=r_version_number

      (self.wb, self.wp, self.gamma, self.pfac,
       self.rmax_surf, self.rmin_surf, self.zmax_surf)=sscan(fl[1])

      if r_version_number > 8.00:
          (self.nfp, self.ns, self.mpol, self.ntor, self.mnmax, self.mnmax_nyq,
           self.itfsq, self.niter, self.iasym, self.ireconstruct, self.ier_flag)=sscan(fl[2],dtype=int)
      else:    
          (self.nfp, self.ns, self.mpol, self.ntor, self.mnmax, self.itfsq, self.niter,
           self.iasym, self.ireconstruct, self.ier_flag)=sscan(fl[2],dtype=int)
          self.mnmax_nyq=self.mnmax #both are the same
      self.StelSym=(self.iasym!=1)
      if self.mnmax_nyq==self.mnmax:
        mpol_nyq=self.mpol
        ntor_nyq=self.ntor
      else:    
        mpol_nyq=self.mpol*2+1  
        ntor_nyq=self.ntor*2  
          
      self.lasym__logical__ =(self.iasym==1)
      self.lrecon__logical__=(self.ireconstruct==1)
      self.lfreeb__logical__=None
      self.lrfp__logical__  =None

      (self.imse2_minus_1,self.itse, self.nbsets, self.nobd,
       self.nextcur, self.nstore_seq)=sscan(fl[3],dtype=int)
      i=4 #line index
      if self.nbsets > 0:
        i+=1 #  these values are not stored in the data structure.
      self.mgrid_file=fl[i].rstrip()
      
      i+=1
      self.xm=np.zeros((self.mnmax),dtype=int)
      self.xn=np.zeros((self.mnmax),dtype=int)
      self.xm_nyq=np.zeros((self.mnmax_nyq),dtype=int)
      self.xn_nyq=np.zeros((self.mnmax_nyq),dtype=int)
      
      self.rmnc=np.zeros(((self.ntor*2+1)*self.mpol,self.ns))
      self.zmns=np.zeros(((self.ntor*2+1)*self.mpol,self.ns))
      self.lmns=np.zeros(((self.ntor*2+1)*self.mpol,self.ns))
      
      self.bmnc=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
      self.gmnc=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
      self.bsubumnc=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
      self.bsubvmnc=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
      self.bsubsmns=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
      self.bsupumnc=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
      self.bsupvmnc=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
      if r_version_number <= 8.00:
        self.currvmnc=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
        
      if self.iasym==1:
        self.rmns=np.zeros(((self.ntor*2+1)*self.mpol,self.ns))
        self.zmnc=np.zeros(((self.ntor*2+1)*self.mpol,self.ns))
        self.lmnc=np.zeros(((self.ntor*2+1)*self.mpol,self.ns))

        self.bmns=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
        self.gmns=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
        self.bsubumns=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
        self.bsubvmns=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
        self.bsubsmnc=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
        self.bsupumns=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))
        self.bsupvmns=np.zeros(((ntor_nyq*2+1)*mpol_nyq,self.ns))

      #Start with rmn,zmn,lmn
      if r_version_number > 8.00:
        for j in range(self.ns):
          lind=0
          for m in range(self.mpol):
            nmin=-self.ntor
            if m==0:
              nmin=0
            for n in range(nmin,self.ntor+1):
              if j==0:
                (md,nd)=sscan(fl[i])
                i+=1
              (self.rmnc[lind,j],self.zmns[lind,j],self.lmns[lind,j])=sscan(fl[i])
              i+=1
              self.xm[lind]=m
              self.xn[lind]=n*self.nfp
              if self.iasym == 1:
                (self.rmns[lind,j],self.zmnc[lind,j],self.lmnc[lind,j])=sscan(fl[i])
                i+=1  
              lind+=1
            #end for n
          #end for m

          #Now load the Nyqvist grid quantities
          lind=0
          for m in range(mpol_nyq):
            nmin=-ntor_nyq
            if m==0:
              nmin=0
            for n in range(nmin,ntor_nyq+1):
              if j==0:
                (md,nd)=sscan(fl[i])
                i+=1
              (self.bmnc[lind,j],self.gmnc[lind,j],self.bsubumnc[lind,j],self.bsubvmnc[lind,j],
               self.bsubsmns[lind,j],self.bsupumnc[lind,j],self.bsupvmnc[lind,j])=sscan(fl[i])
              i+=1
              self.xm_nyq[lind]=m
              self.xn_nyq[lind]=n*self.nfp
              if self.iasym == 1:
                (self.bmns[lind,j],self.gmns[lind,j],self.bsubumns[lind,j],self.bsubvmns[lind,j],
                 self.bsubsmnc[lind,j],self.bsupumns[lind,j],self.bsupvmns[lind,j])=sscan(fl[i])
                i+=1
              lind+=1
            #end for n
          #end for m
        #end for j  
      else: #r_version_number <= 8.00:
        for j in range(self.ns):
          lind=0
          for m in range(self.mpol):
            nmin=-self.ntor
            if m==0:
              nmin=0
            for n in range(nmin,self.ntor+1):
              if j==0:
                (md,nd)=sscan(fl[i])
                i+=1
              (self.rmnc[lind,j],self.zmns[lind,j],self.lmns[lind,j],
               self.bmnc[lind,j],self.gmnc[lind,j],self.bsubumnc[lind,j],self.bsubvmnc[lind,j],
               self.bsubsmns[lind,j],self.bsupumnc[lind,j],self.bsupvmnc[lind,j],self.currvmnc[lind,j])=sscan(fl[i])
              i+=1
              self.xm[lind]=m
              self.xn[lind]=n*self.nfp
              if self.iasym == 1:
                (self.rmns[lind,j],self.zmnc[lind,j],self.lmnc[lind,j],
                 self.bmns[lind,j],self.gmns[lind,j],self.bsubumns[lind,j],self.bsubvmns[lind,j],
                 self.bsubsmnc[lind,j],self.bsupumns[lind,j],self.bsupvmns[lind,j])=sscan(fl[i])
                i+=1  
              lind+=1
            #end for n
          #end for m
          self.xm_nyq=self.xm
          self.xn_nyq=self.xn
        #end for j  
      #end if version
                               
      self.iotaf   =np.zeros((self.ns))
      self.presf   =np.zeros((self.ns))
      self.phipf   =np.zeros((self.ns))
      self.phi     =np.zeros((self.ns))
      self.jcuru   =np.zeros((self.ns))
      self.jcurv   =np.zeros((self.ns))
      self.iotas   =np.zeros((self.ns))
      self.mass    =np.zeros((self.ns))
      self.pres    =np.zeros((self.ns))
      self.beta_vol=np.zeros((self.ns))
      self.phip    =np.zeros((self.ns))
      self.buco    =np.zeros((self.ns))
      self.bvco    =np.zeros((self.ns))
      self.vp      =np.zeros((self.ns))
      self.over_r  =np.zeros((self.ns))
      self.specw   =np.zeros((self.ns))

      # reading of profiles.
      # This part has been copied from the read_wout_mod.f-file of LIBSTELL.
      if r_version_number <= 6.05+eps_w:
        quants=['iotas','mass','pres','phip','buco','bvco','phi','vp','over_r','jcuru','jcurv','specw']
        tmp=sscan(fl[i])
        i+=1    
        for qi in range(len(quants)):
          setattr(self,quants[qi],np.insert(tmp[(self.ns-1)*qi:(self.ns-1)*(qi+1)],0,0.0))
        (self.aspect, self.betatot, self.betapol, self.betator, self.betaxis, self.b0)=sscan(fl[i])
        i+=1
      elif r_version_number <= 6.20+eps_w:
        quants=['iotas','mass','pres','beta_vol','phip','buco','bvco','phi','vp','over_r','jcuru','jcurv','specw']
        tmp=sscan(fl[i])
        i+=1    
        for qi in range(len(quants)):
          setattr(self,quants[qi],np.insert(tmp[(self.ns-1)*qi:(self.ns-1)*(qi+1)],0,0.0))
        (self.aspect, self.betatot, self.betapol, self.betator, self.betaxis, self.b0)=sscan(fl[i])
        i+=1
      elif r_version_number <= 6.95+eps_w:
        quants=['iotas','mass','pres','beta_vol','phip','buco','bvco','phi','vp','over_r','jcuru','jcurv','specw']
        tmp=sscan(fl[i])
        i+=1    
        for qi in range(len(quants)):
          setattr(self,quants[qi],np.insert(tmp[(self.ns-1)*qi:(self.ns-1)*(qi+1)],0,0.0))
        (self.aspect, self.betatot, self.betapol, self.betator, self.betaxis, self.b0)=sscan(fl[i])
        i+=1
      else:
        quants=['iotaf','presf','phipf','phi','jcuru','jcurv']
        tmp=sscan(fl[i])
        i+=1     
        for qi in range(len(quants)):
          setattr(self,quants[qi],tmp[self.ns*qi:self.ns*(qi+1)])
          
        quants=['iotas','mass','pres','beta_vol','phip','buco','bvco','vp','over_r','specw']
        tmp=sscan(fl[i])
        i+=1     
        for qi in range(len(quants)):
          setattr(self,quants[qi],np.insert(tmp[(self.ns-1)*qi:(self.ns-1)*(qi+1)],0,0.0))
          
        (self.aspect, self.betatot, self.betapol, self.betator, self.betaxis, self.b0)=sscan(fl[i])
        i+=1

      if r_version_number > 6.10+eps_w:
        self.signgs=sscan(fl[i])[0]
        i+=1
        #print('len='+str(len(fl[i])))
        self.input_extension=fl[i].rstrip()
        i+=1
        (self.IonLarmor, self.VolAvgB, self.rbtor0, self.rbtor, self.ctor,
         self.Aminor_p, self.Rmajor_p, self.volume_p)=sscan(fl[i])
        i+=1
        
      # Mercier criterion  
      self.DMerc =np.zeros((self.ns))
      self.DShear=np.zeros((self.ns))
      self.DWell =np.zeros((self.ns))
      self.DCurr =np.zeros((self.ns))
      self.DQeod =np.zeros((self.ns))
      self.equif =np.zeros((self.ns))
      quants=['DMerc','DShear','DWell','DCurr','DGeod','equif']
      tmp=sscan(fl[i])
      #print('lentot='+str(len(tmp)))
      #print('len/6='+str(len(tmp)/6.0))
      i+=1
      if (r_version_number > 5.10+eps_w) and (r_version_number< 6.20-eps_w):
        for qi in range(len(quants)):
          setattr(self,quants[qi],np.insert(tmp[(self.ns-1)*qi:(self.ns-1)*(qi+1)],0,0.0))
      elif r_version_number >= 6.20-eps_w:
        for qi in range(len(quants)):
          setattr(self,quants[qi],np.append(np.insert(tmp[(self.ns-2)*qi:(self.ns-2)*(qi+1)],0,0.0),[0.0]))

      if self.nextcur>0:
        self.extcur=sscan(fl[i])
        i+=1
        self.curlabel=fl[i]
        i+=1
      else:
        self.curlabel=None
        self.extcur=None
        
      self.unidentified=[]  
      while i<len(fl):
        self.unidentified.append(fl[i])
        i+=1
      #print('len(unidentified)='+str(len(self.unidentified)))
      #print('i='+str(i))
      #print('len(fl)='+str(len(fl)))
      #print('fl[-1]='+str(fl[-1]))
      #print('len(sscan(self.unidentified[0]))='+str(len(sscan(self.unidentified[0]))))      
      #print('len(sscan(self.unidentified[1]))='+str(len(sscan(self.unidentified[1]))))      

    else:
      sys.exit('Error: Unknown file extension in the file '+filename+' !') 
        
  ####################################################################

  def write(self,filename):
    if not(filename.endswith('.txt')):
      sys.exit('Error in writing vmecgeom: Only .txt files supported so far.')
    #if self.version_<=8.00:
    #  sys.exit('Error in writing vmecgeom: Only saving of .txt files with version > 8.00 implemented. '+
    #           'This file hase is version {:4.2f}'.format(self.version_))
    eps_w = 1.0e-4 #a small number.
    frme=' {:.15E}'
    frmf15=' {:.15f}'
    frmf16=' {:.16f}'
    frmd='{:d}'
    def frmfune(vals):
        lista=[None]*len(vals)
        for i in range(len(lista)):
          if vals[i]==0.0 or abs(vals[i])>=10.0**5 or abs(vals[i])<0.1:
            lista[i]=frme.format(vals[i])
          else:
            p=0
            while abs(vals[i])>=10.0**p:
              p+=1
            thisform=' {:.'+str(16-p)+'f}'  
            lista[i]=thisform.format(vals[i]) 
          #elif abs(vals[i])>=1.0 and abs(vals[i])<10.0:
          #  lista[i]=frmf15.format(vals[i])
          #elif abs(vals[i])<1.0:  
          #  lista[i]=frmf16.format(vals[i])            
          #elif abs(vals[i])<100.0: 
          #  lista[i]=' {:.14f}'.format(vals[i])
          #else:
          #  lista[i]=' {:.13f}'.format(vals[i])    
        return ''.join(lista)    
    def frmfunf(vals):
        lista=[None]*len(vals)
        for i in range(len(lista)):
            lista[i]=frmf15.format(vals[i])
        return ''.join(lista)    
    def frmfund(vals):
        lista=[None]*len(vals)
        for i in range(len(lista)):
            lista[i]=' '+frmd.format(vals[i])
        return ''.join(lista)    
    try:
        f=open(filename,'w')
    except:
        sys.error('Error in writing vmecgeom: Could not open the file '+filename+' for writing!')    
    
    fl=['VMEC VERSION = {:4.2f}'.format(self.version_)+'\n']
    fl.append(frmf15.format(self.wb)+frmfune([self.wp, self.gamma])+
                        frmfunf([self.pfac,self.rmax_surf, self.rmin_surf])+
                        frmf16.format(self.zmax_surf)+'\n')
    if self.version_ > 8.00:
      fl.append(frmfund([self.nfp, self.ns, self.mpol, self.ntor, self.mnmax, self.mnmax_nyq,
                         self.itfsq, self.niter, self.iasym, self.ireconstruct, self.ier_flag])+'\n')
    else:
      fl.append(frmfund([self.nfp, self.ns, self.mpol, self.ntor, self.mnmax, self.itfsq, self.niter,
                         self.iasym, self.ireconstruct, self.ier_flag])+'\n')    
    fl.append(frmfund([self.imse2_minus_1,self.itse, self.nbsets, self.nobd,
                       self.nextcur, self.nstore_seq])+'\n')
    if self.nbsets>0:
        sys.exit('Error in writing vmecgeom: The variable nbsets>0. Unknown case.')
    if self.mnmax_nyq==self.mnmax:
        mpol_nyq=self.mpol
        ntor_nyq=self.ntor
    else:    
        mpol_nyq=self.mpol*2+1  
        ntor_nyq=self.ntor*2  
    fl.append(self.mgrid_file+'\n')
    #print(self.xm)
    #print(type(self.xm))
    #fl.append(' '.join(frmfund([self.xm,self.xn])))
    f.writelines(fl)
    fl=[]
    #Start with rmn,zmn,lmn
    if self.version_ > 8.00:
      for j in range(self.ns):
        lind=0
        for m in range(self.mpol):
          nmin=-self.ntor
          if m==0:
            nmin=0
          for n in range(nmin,self.ntor+1):
            if j==0:
              fl.append(frmfund([self.xm[lind],self.xn[lind]])+'\n')
            fl.append(frmfune((self.rmnc[lind,j],self.zmns[lind,j],self.lmns[lind,j]))+'\n')
            #self.xm[lind]=m
            #self.xn[lind]=n*self.nfp
            if self.iasym == 1:
              fl.append(frmfune((self.rmns[lind,j],self.zmnc[lind,j],self.lmnc[lind,j]))+'\n')
            lind+=1
          #end for n
        #end for m
        f.writelines(fl)
        fl=[]
        
        #Now load the Nyqvist grid quantities
        lind=0
        for m in range(mpol_nyq):
          nmin=-ntor_nyq
          if m==0:
            nmin=0
          for n in range(nmin,ntor_nyq+1):
            if j==0:
              fl.append(frmfund([self.xm[lind],self.xn[lind]])+'\n')
            fl.append(frmfune((self.bmnc[lind,j],self.gmnc[lind,j],self.bsubumnc[lind,j],self.bsubvmnc[lind,j],
                               self.bsubsmns[lind,j],self.bsupumnc[lind,j],self.bsupvmnc[lind,j]))+'\n')
            #self.xm_nyq[lind]=m
            #self.xn_nyq[lind]=n*self.nfp
            if self.iasym == 1:
              fl.append(frmfune((self.bmns[lind,j],self.gmns[lind,j],self.bsubumns[lind,j],self.bsubvmns[lind,j],
                                 self.bsubsmnc[lind,j],self.bsupumns[lind,j],self.bsupvmns[lind,j]))+'\n')
            lind+=1
          #end for n
        #end for m
        f.writelines(fl)
        fl=[]        
      #end for j  

    else: #r_version_number <= 8.00:
      for j in range(self.ns):
        lind=0
        for m in range(self.mpol):
          nmin=-self.ntor
          if m==0:
            nmin=0
          for n in range(nmin,self.ntor+1):
            if j==0:
              fl.append(frmfund([self.xm[lind],self.xn[lind]])+'\n')
            fl.append(frmfune((self.rmnc[lind,j],self.zmns[lind,j],self.lmns[lind,j],
             self.bmnc[lind,j],self.gmnc[lind,j],self.bsubumnc[lind,j],self.bsubvmnc[lind,j],
             self.bsubsmns[lind,j],self.bsupumnc[lind,j],self.bsupvmnc[lind,j],self.currvmnc[lind,j]))+'\n')
            #self.xm[lind]=m
            #self.xn[lind]=n*self.nfp
            if self.iasym == 1:
              fl.append(frmfune((self.rmns[lind,j],self.zmnc[lind,j],self.lmnc[lind,j],
               self.bmns[lind,j],self.gmns[lind,j],self.bsubumns[lind,j],self.bsubvmns[lind,j],
               self.bsubsmnc[lind,j],self.bsupumns[lind,j],self.bsupvmns[lind,j]))+'\n')
            lind+=1
          #end for n
        #end for m
        f.writelines(fl)
        fl=[]        
      #end for j  
    #end if version
        
    if self.version_ <= 6.05+eps_w:
      quants=['iotas','mass','pres','phip','buco','bvco','phi','vp','over_r','jcuru','jcurv','specw']
      tmp=np.array([])
      for qi in range(len(quants)):
        tmp=np.append(tmp,getattr(self,quants[qi])[1:])#,np.insert(tmp[(self.ns-1)*qi:(self.ns-1)*(qi+1)],0,0.0))
      fl.append(frmfune(tmp)+'\n')  
      fl.append(frmfune((self.aspect, self.betatot, self.betapol, self.betator, self.betaxis, self.b0))+'\n')
    elif self.version_ <= 6.20+eps_w:
      quants=['iotas','mass','pres','beta_vol','phip','buco','bvco','phi','vp','over_r','jcuru','jcurv','specw']
      tmp=np.array([])
      for qi in range(len(quants)):
        tmp=np.append(tmp,getattr(self,quants[qi])[1:])#,np.insert(tmp[(self.ns-1)*qi:(self.ns-1)*(qi+1)],0,0.0))
      fl.append(frmfune(tmp)+'\n')  
      fl.append(frmfune((self.aspect, self.betatot, self.betapol, self.betator, self.betaxis, self.b0))+'\n')
    elif self.version_ <= 6.95+eps_w:
      quants=['iotas','mass','pres','beta_vol','phip','buco','bvco','phi','vp','over_r','jcuru','jcurv','specw']
      tmp=np.array([])
      for qi in range(len(quants)):
        tmp=np.append(tmp,getattr(self,quants[qi])[1:])#,np.insert(tmp[(self.ns-1)*qi:(self.ns-1)*(qi+1)],0,0.0))
      fl.append(frmfune(tmp)+'\n')
      fl.append(frmfune((self.aspect, self.betatot, self.betapol, self.betator, self.betaxis, self.b0))+'\n')
    else:
      quants=['iotaf','presf','phipf','phi','jcuru','jcurv']
      tmp=np.array([])
      for qi in range(len(quants)):
        tmp=np.append(tmp,getattr(self,quants[qi]))#,tmp[self.ns*qi:self.ns*(qi+1)])
      fl.append(frmfune(tmp)+'\n')
      quants=['iotas','mass','pres','beta_vol','phip','buco','bvco','vp','over_r','specw']
      tmp=np.array([])
      for qi in range(len(quants)):
        tmp=np.append(tmp,getattr(self,quants[qi])[1:])#,np.insert(tmp[(self.ns-1)*qi:(self.ns-1)*(qi+1)],0,0.0))
      fl.append(frmfune(tmp)+'\n')
        
      fl.append(frmfune((self.aspect, self.betatot, self.betapol, self.betator, self.betaxis, self.b0))+'\n')

    if self.version_ > 6.10+eps_w:
      fl.append(' {:d}'.format(int(self.signgs))+'\n')
      if len(self.input_extension)<=101:
        fl.append('{:100}'.format(self.input_extension)+'\n')
      else:
        fl.append(self.input_extension+'\n')    
      fl.append(frmfune((self.IonLarmor, self.VolAvgB, self.rbtor0, self.rbtor, self.ctor,
                         self.Aminor_p, self.Rmajor_p, self.volume_p))+'\n')
      
    f.writelines(fl)
    fl=[]
    
    # Mercier criterion  
    quants=['DMerc','DShear','DWell','DCurr','DGeod','equif']
    #tmp=sscan(fl[i])
    #print('len='+str(len(tmp)/6.0))
    if (self.version_ > 5.10+eps_w) and (self.version_< 6.20-eps_w):
      tmp=np.array([])  
      for qi in range(len(quants)):
        tmp=np.append(tmp,getattr(self,quants[qi])[1:])#,np.insert(tmp[(self.ns-1)*qi:(self.ns-1)*(qi+1)],0,0.0))
      fl.append(frmfune(tmp)+'\n')
    elif self.version_ >= 6.20-eps_w:
      tmp=np.array([])    
      for qi in range(len(quants)):
        tmp=np.append(tmp,getattr(self,quants[qi])[1:-1])
        #setattr(self,quants[qi],np.append(np.insert(tmp[(self.ns-2)*qi:(self.ns-2)*(qi+1)],0,0.0),[0.0]))
      fl.append(frmfune(tmp)+'\n')
    if self.nextcur>0:
      fl.append(frmfune((self.extcur))+'\n')  
      fl.append(self.curlabel)
    f.writelines(fl)
    if hasattr(self,'unidentified'):
      for i in range(len(self.unidentified)):
        f.write(self.unidentified[i])  
    f.close()      
    #sys.exit('Under construction!!')  

  #####################################################################################################    
  def disp(self,verbose=False):
    frm='{:8.4f}'
    frme='{:10.4e}'
    frmi='{:8d}'
    frs='{:11s}'
    if verbose:
        print('-----------------------------------------------')
        print('File info')
        print('-----------------------------------------------')
        print('filename       = '+self.filename)
        print('version_       = '+str(self.version_))
        print('input_extension= '+self.input_extension)
        print('mgrid_file     = '+self.mgrid_file)
        if hasattr(self,'pcurr_type'):
          print('pcurr_type     = '+self.pcurr_type)
        if hasattr(self,'pmass_type'):
          print('pmass_type     = '+self.pmass_type)
        if hasattr(self,'piota_type'):
          print('piota_type     = '+self.piota_type)
        if hasattr(self,'mgrid_mode'):
          print('mgrid_mode     = '+self.mgrid_mode)
        print(' ' )
    print('-----------------------------------------------')
    print('Scalars:')
    print('-----------------------------------------------')
    if verbose:
        print('wb                   = '+frm.format(self.wb)+' : Magnetic energy')
        print('wp                   = '+frm.format(self.wp)+' : Plasma energy')
        print('gamma                = '+frm.format(self.gamma)+' :')
        print('rmax_surf            = '+frm.format(self.rmax_surf)+' :')
        print('rmin_surf            = '+frm.format(self.rmin_surf)+' :')
        print('zmax_surf            = '+frm.format(self.zmax_surf)+' :')
    print('StelSym              =    '+str(self.StelSym))
    print('nfp                  = '+frm.format(self.nfp)+' : Number of field periods')
    print('ns                   = '+frmi.format(self.ns)+' : Number of flux surfaces')
    print('mpol                 = '+frmi.format(self.mpol)+' :')
    print('ntor                 = '+frmi.format(self.ntor)+' :')
    print('mnmax                = '+frmi.format(self.mnmax)+' :')
    print('mnmax_nyq            = '+frmi.format(self.mnmax_nyq)+' :')
    if verbose:
        print('niter                = '+frmi.format(self.niter)+' :')
        print('itfsq                = '+frm.format(self.itfsq)+' :')
        print('lasym__logical__     =     '+str(self.lasym__logical__)+' :')
        print('lrecon__logical__    =    '+str(self.lrecon__logical__)+' :')
        print('lfreeb__logical__    =    '+str(self.lfreeb__logical__)+' :')
        if hasattr(self,'lrfp__logical__'):
          print('lrfp__logical__      =    '+str(self.lrfp__logical__)+' :')
        print('ier_flag             = '+frmi.format(self.ier_flag)+' :')
        print('aspect               = '+frm.format(self.aspect)+' :')
    if hasattr(self,'betatotal'):
        print('betatotal            = '+frm.format(self.betatotal)+' :')
    print('betapol              = '+frm.format(self.betapol)+' :')
    print('betator              = '+frm.format(self.betator)+' :')
    print('betaxis              = '+frm.format(self.betaxis)+' :')
    print('b0                   = '+frm.format(self.b0)+' :')
    if verbose:
        if hasattr(self,'rbtor0'):
            print('rbtor0               = '+frm.format(self.rbtor0)+' :')
        if hasattr(self,'rbtor'):
            print('rbtor                = '+frm.format(self.rbtor)+' :')
        if hasattr(self,'signs'):
            print('signgs               = '+frmi.format(self.signgs)+' :')
        if hasattr(self,'IonLarmor'):
            print('IonLarmor            = '+frm.format(self.IonLarmor)+' :')
        if hasattr(self,'volavgB'):
            print('volavgB              = '+frm.format(self.volavgB)+' :')
        if hasattr(self,'ctor'):
            print('ctor                 = '+frm.format(self.ctor)+' :')
    if hasattr(self,'Aminor_p'):
        print('Aminor_p             = '+frm.format(self.Aminor_p)+' : Minor radius')
    if hasattr(self,'Rmajor_p'):
        print('Rmajor_p             = '+frm.format(self.Rmajor_p)+' : Major radius')
    if hasattr(self,'volume_p'):
        print('volume_p             = '+frm.format(self.volume_p)+' :')
    if verbose:
        if hasattr(self,'ftolv'):
          print('ftolv                = '+frm.format(self.ftolv)+' : Volumetric force residual')
        if hasattr(self,'fsql'):
          print('fsql                 = '+frm.format(self.fsql)+' : Lambda force residual')
        if hasattr(self,'fsqr'):
          print('fsqr                 = '+frm.format(self.fsqr)+' : R force residual')
        if hasattr(self,'fsqz'):
          print('fsqz                 = '+frm.format(self.fsqz)+' : Z force residual')
        if hasattr(self,'nextcur'):
          print('nextcur              = '+frmi.format(self.nextcur)+' :')

    print(' ')
    if verbose:
        print('-----------------------------------------------')
        print('Arrays:')
        print('-----------------------------------------------')
        print('raxis_cc        '+frs.format(self.xm_nyq.shape)+  ': raxis (cosnv)')
        if hasattr(self,'zaxis_cs'):
          print('zaxis_cs        '+frs.format(self.zaxis_cs.shape)+': zaxis (sinnv)')
        if not(self.StelSym):
          if hasattr(self,'raxis_cs'):
            print('raxis_cs        '+frs.format(self.raxis_cs.shape)+': raxis (sinnv)')
          if hasattr(self,'zaxis_cc'):
            print('zaxis_cc        '+frs.format(self.zaxis_cc.shape)+': zaxis (cosnv)')
        if hasattr(self,'am'):
          print('am              '+frs.format(self.am.shape)+': mass profile [Pa]')
        if hasattr(self,'ac'):
          print('ac              '+frs.format(self.ac.shape)+': normalised current density [dimensionless] (used if NCURR=1)')
        if hasattr(self,'ai'):
          print('ai              '+frs.format(self.ai.shape)+': iota profile (used if NCURR=0)')
        if hasattr(self,'am_aux_s'):
          print('am_aux_s        '+frs.format(self.am_aux_s.shape)+': Spline knots')
        if hasattr(self,'am_aux_f'):
          print('am_aux_f        '+frs.format(self.am_aux_f.shape)+': Spline values')
        if hasattr(self,'ai_aux_s'):
          print('ai_aux_s        '+frs.format(self.ai_aux_s.shape)+': Spline knots')
        if hasattr(self,'ai_aux_f'):
          print('ai_aux_f        '+frs.format(self.ai_aux_f.shape)+': Spline values')
        if hasattr(self,'ac_aux_s'):
          print('ac_aux_s        '+frs.format(self.ac_aux_s.shape)+': Spline knots')
        if hasattr(self,'ac_aux_f'):
          print('ac_aux_f        '+frs.format(self.ac_aux_f.shape)+': Spline values')
        if hasattr(self,'fsqt'):
          print('fsqt            '+frs.format(self.fsqt.shape)+':')
        if hasattr(self,'wdot'):
          print('wdot            '+frs.format(self.wdot.shape)+': rate of change of energy')
        if hasattr(self,'extcur'):
          if not(self.extcur is None):
            print('extcur          '+frs.format(self.extcur.shape)+':')
        if hasattr(self,'curlabel'):
          if not(self.curlabel is None):
            print('curlabel = '+self.curlabel)
        
    print('-----------------------------------------------')
    print('Radial arrays:')
    print('-----------------------------------------------')
    print('iotaf        '+frs.format(self.iotaf.shape)+   ': q-factor on full mesh')
    if hasattr(self,'q_factor'):
      if not((np.isnan(self.q_factor)).all()):
        print('q_factor     '+frs.format(self.q_factor.shape)+':')
    print('presf        '+frs.format(self.presf.shape)+': pressure on full mesh [Pa]')
    print('phi          '+frs.format(self.phi.shape)+  ': Toroidal flux on full mesh [Wb]')
    print('phipf        '+frs.format(self.phipf.shape)+': d(phi)/ds: Toroidal flux deriv on full mesh')
    if hasattr(self,'chi'):
      if not((np.isnan(self.chi)).all()):
        print('chi          '+frs.format(self.chi.shape)+  ': Poloidal flux on full mesh [Wb]')
        print('chipf        '+frs.format(self.chipf.shape)+': d(chi)/ds: Poroidal flux deriv on full mesh')
    print('jcuru        '+frs.format(self.jcuru.shape)+':')
    print('jcurv        '+frs.format(self.jcurv.shape)+':')
    print('iotas        '+frs.format(self.iotas.shape)+': iota half')
    print('mass         '+frs.format(self.mass.shape)+ ': mass half')
    print('pres         '+frs.format(self.pres.shape)+ ': pressure half [Pa]')
    print('beta_vol     '+frs.format(self.beta_vol.shape)+':')
    print('buco         '+frs.format(self.buco.shape)+':')
    print('bvco         '+frs.format(self.bvco.shape)+':')
    print('vp           '+frs.format(self.vp.shape)+':')
    print('specw        '+frs.format(self.specw.shape)+':')
    if hasattr(self,'phips'):
      print('phips        '+frs.format(self.phips.shape)+':')
    if hasattr(self,'over_r'):
      print('over_r       '+frs.format(self.over_r.shape)+':')
    if hasattr(self,'jdotb'):
      print('jdotb        '+frs.format(self.jdotb.shape)+':')
    if hasattr(self,'bdotgradv'):
      print('bdotgradv    '+frs.format(self.bdotgradv.shape)+':')
    if verbose:
        print('DMerc        '+frs.format(self.DMerc.shape)+':')
        print('DShear       '+frs.format(self.DShear.shape)+':')
        print('DWell        '+frs.format(self.DWell.shape)+':')
        print('DCurr        '+frs.format(self.DCurr.shape)+':')
        print('DGeod        '+frs.format(self.DGeod.shape)+':')
        print('equif        '+frs.format(self.equif.shape)+':')
    
    print(' ')
    print('-----------------------------------------------')
    print('Fourier quantities:')
    print('-----------------------------------------------')
    
    print('xm        '+frs.format(self.xm.shape)+      ': Poloidal mode numbers')
    print('xn        '+frs.format(self.xm.shape)+      ': Toroidal mode numbers')
    print('xm_nyq    '+frs.format(self.xm_nyq.shape)+  ': Poloidal mode numbers (Nyquist)')
    print('xn_nyq    '+frs.format(self.xm_nyq.shape)+    ': Toroidal mode numbers (Nyquist)')
    print('rmnc      '+frs.format(self.rmnc.shape)+': cosmn component of cylindrical R, full mesh [m]')
    print('zmns      '+frs.format(self.zmns.shape)+': sinmn component of cylindrical Z, full mesh [m]')
    print('lmns      '+frs.format(self.lmns.shape)+': sinmn component of lambda, half mesh')
    print('gmnc      '+frs.format(self.gmnc.shape)+': cosmn component of jacobian, half mesh')
    print('bmnc      '+frs.format(self.bmnc.shape)+': cosmn component of mod-B, half mesh')
    print('bsubumnc  '+frs.format(self.bsubumnc.shape)+': cosmn covariant u-component of B, half mesh')
    print('bsubvmnc  '+frs.format(self.bsubvmnc.shape)+': cosmn covariant v-component of B, half mesh')
    print('bsubsmns  '+frs.format(self.bsubsmns.shape)+': sinmn covariant s-component of B, full mesh')
    print('bsupumnc  '+frs.format(self.bsupumnc.shape)+':')
    print('bsupvmnc  '+frs.format(self.bsupvmnc.shape)+':')
    if not(self.StelSym):
      print('rmns      '+frs.format(self.rmns.shape)+': sinmn component of cylindrical R, full mesh [m]')
      print('zmnc      '+frs.format(self.zmnc.shape)+': cosmn component of cylindrical Z, full mesh [m]')
      print('lmnc      '+frs.format(self.lmnc.shape)+': cosmn component of lambda, half mesh')
      print('gmns      '+frs.format(self.gmns.shape)+': sinmn component of jacobian, half mesh')
      print('bmns      '+frs.format(self.bmns.shape)+': sinmn component of mod-B, half mesh')
      print('bsubumns  '+frs.format(self.bsubumns.shape)+': sinmn covariant u-component of B, half mesh')
      print('bsubvmns  '+frs.format(self.bsubvmns.shape)+': sinmn covariant v-component of B, half mesh')
      print('bsubsmnc  '+frs.format(self.bsubsmnc.shape)+': cosmn covariant s-component of B, full mesh')
      print('bsupumns  '+frs.format(self.bsupumns.shape)+':')
      print('bsupvmns  '+frs.format(self.bsupvmns.shape)+':')

    print('-----------------------------------------------')



    
 







      
   
