#!/usr/bin/env python
from __future__ import division
import numpy as np
import sys, copy
import datetime

import mnlist
import mnmat
from vmecgeom import vmecgeom
import fluxcoorddiscr


class headertxt:
  def __init__(self):
      self.maincomment=''
      self.globalvars=''
      self.surfvars=''
      self.surfvarunits=''
      self.datavars=''

class Bfiltr:
  def __init__(self):
      self.min_Bmn=0
      self.max_m=np.inf
      self.maxabs_n=np.inf

class bcgeom:
    
  def __init__(self,input,min_Bmn=0,max_m=np.inf,maxabs_n=np.inf,
               symmetry='unknown',signcorr=2):
    if isinstance(input,str):
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
           return np.fromstring(strng,dtype,count,sep=' ')

       if signcorr==1:
           newsigncorrectionmethod=False
       elif signcorr==2:
           newsigncorrectionmethod=True
       else:
           sys.exit('Sign correction method for JG files not recognised!')


       if filename[-4:-1]=='.dat':
           #Henning Maassberg type file
           filetype='HM'
       elif filename[-3:-1]=='.bc':
           #Joachim Geiger (or Erika Strumberger) type file
           filetype='JG'
       else: #default to JG
           filetype='JG'

       if filetype=='JG': 
         f = open(filename, 'r')
         f.seek(-1,2)     # go to the file end.
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
         self.Nperiods   = header_df[3]
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
           print '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%5i/%5i'% (rind+1,self.nsurf),
           if not(YTstyle): #%isempty(strfind(self.headertext.surfvars,'[A]'))
               self.headertext.surfvarunits=f.readline() #unit line only in JG files
           #print surfvarunits: '+self.headertext.surfvarunits
           surfheader=sscan(f.readline()) #fscanf(fid,'%f',6)
           torfluxnorm=np.append(torfluxnorm, surfheader[0])
           iota =np.append(iota,surfheader[1])
           Bphi=np.append(Bphi,surfheader[2]*self.Nperiods/2.0/np.pi*(4*np.pi*1e-7)) #Tesla*meter
           Btheta=np.append(Btheta,surfheader[3]/2.0/np.pi*(4*np.pi*1e-7))              #Tesla*meter
           dpds=np.append(dpds,surfheader[4])
           dVdsoverNper=np.append(dVdsoverNper,surfheader[5])

           #f.readline() #%just skip the return character
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
                   #print 'eof found '+tmp_str
                   proceed=False
                   endoffile=True
                   #print 'found end of file'
               if ('s' in tmp_str): #Next flux surface has been reached
                   proceed=False
                   #print 'found next surface'
               else:
                   #print tmp_str
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
         f.close()
         print '' #go to new line
         
         if any([a>0 for a in dVdsoverNper]):
           sys.exit('The coordinate system in the Boozer '+
                    'file should be left handed, but it has '+
                    'a positive Jacobian. Something is wrong!')

         if self.torfluxtot*Bphi[0]>0:
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
     
           f = open(filename, 'r')
           f.seek(-1,2)     # go to the file end.
           eof = f.tell()   # get the end of file location
           f.seek(0,0)      # go back to file beginning
           tmp_str=f.readline()
           while tmp_str[1]=='c' or tmp_str[1]=='C':  #Skip comment line
               tmp_str=f.readline() 

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
             
           endoffile=False
           rind=-1

           while not endoffile:
               rind=rind+1
               surfheader1     = sscan(tmp_str,count=8)
               radii           = np.append(radii,surfheader1[0]/100) #(cm)
               iota            = np.append(iota,surfheader1[1])
               Nperiods        = surfheader1[2]
               minorradiusW7AS = surfheader1[3]/100 #(cm)
               majorradiusLastbcR00 = surfheader1[4]/100 #(cm)
               if len(surfheader1)>5:
                   torfluxnorm = np.append(torfluxnorm,surfheader1[5])
                   R00         = np.append(R00,surfheader1[7]/100)
               else:
                   torfluxnorm = np.append(torfluxnorm,np.nan)
                   R00         = np.append(R00,np.nan)
       
               tmp_str=f.readline()
               if tmp_str[1]=='>':
                   surfheader2=sscan(tmp_str[2:],count=3)
           
                   Bphi      = np.append(Bphi,surfheader2[0]*1.0e6*Nperiods/2.0/pi*(4*pi*1e-7)) #Tesla*meter
                   Btheta    = np.append(Btheta,surfheader2[1]*1.0e6/2.0/pi*(4*pi*1e-7))        #Tesla*meter
                   B00       = np.append(B00,surfheader2[2]) #Tesla
                   tmp_str=f.readline()
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
         
                   tmp_str=f.readline() #get the next line or surface header line
               #end while proceed
               if modeind==-1:
                   sys.exit('no modes found for rind='+str(rind))

               no_of_modes.append(int(modeind+1))

               if f.tell() == eof:
                   endoffile=True

           #end while not endoffile
           f.close()
           
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
               
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    #% Turn a loaded vmec dataset vmecgeom into a bcgeom
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    elif isinstance(input,vmecgeom):
        wout=input
        signchange=float(wout.signgs) #is -1, because vmec is left handed
        self.StelSym=input.StelSym
        self.newsigncorr=True
        self.headertext=headertxt()
        # Note that the following are only comments to what would appear in a file stored with .write()
        # The internal representation in python is right handed and psi_a>0 in the direction of phi!
        if self.StelSym:
            # Use Joachim Geiger's JMC convention as default for stellarator symmetric files:
            # Toroidal flux is counted positive in the direction opposite to the toroidal coordinate.
            self.headertext.maincomment=('CC Converted from VMEC using HS''s python bcgeom.py routines\n'+
               'CC version_       = '+str(wout.version_)+'\n'+
               'CC input_extension= '+wout.input_extension+'\n'+
               'CC The phase convention (m theta - n phi) is used in this file.\n'+
               'CC The coordinate system (r,theta,phi) is left-handed.\n'+
               'CC Toroidal flux is counted positive in the direction opposite to the toroidal coordinate.')
        else:
            # Use Erika Strumberger's convention as default for non-stellarator symmetric files:
            # Toroidal flux is counted positive in the direction of the toroidal coordinate.
            self.headertext.maincomment=('CC Converted from VMEC using HSs routines\n'+
            'CC VMEC version         = '+str(wout.version_)+'\n'+
            'CC VMEC input extension = '+wout.input_extension+'\n'+
            'CC The phase convention (m theta - n phi) is used in this file.\n'+
            'CC The coordinate system (r,theta,phi) is left-handed.\n'+
            'CC Toroidal flux is counted positive in the direction of the toroidal coordinate.')
          
        self.m0b=wout.mpol
        self.n0b=wout.ntor
        self.nsurf= np.nan      #number of radial surfaces, Set this below
        self.Nperiods=wout.ns   #%!< number of field periods
        self.torfluxtot = wout.phi[wout.ns-1]*signchange
        self.psi_a=self.torfluxtot/2.0/np.pi
        self.minorradiusVMEC      = wout.Aminor_p  #minor plasma radius
        self.majorradiusLastbcR00 = np.nan         #Calculate this below (not necessary)
        self.minorradiusW7AS      = np.nan         #Calculate this below (not necessary)
        self.majorradiusVMEC      = wout.Rmajor_p  #major plasma radius
        self.volumeVMEC           = wout.volume_p  #plasma volume
        
        fullgrid_s           = np.array(wout.phi/wout.phi[wout.ns-1]) #full grid
        fullgrid_rnorm       = np.sqrt(fullgrid_s);     #full grid

        skip=wout.skip #=1,this is how many elements are skipped at low radii when going to half grid

        self.s      = (fullgrid_s[skip-1:-1]+fullgrid_s[skip:])/2.0  #half grid
        self.rnorm  = np.sqrt(self.s) #half grid
        self.nsurf  = len(self.s)
        self.dpds   = np.array(np.diff(wout.presf[skip-1:])/np.diff(fullgrid_s[skip-1:]))
        self.Bphi        = wout.bvco[skip:]*signchange #direction switch sign
        self.Btheta      = wout.buco[skip:] #sign change

        fullgrid_iota=wout.iotaf*signchange
        iota   = wout.iotas*signchange #half mesh
        self.iota=iota[skip:]
        
        self.Bfilter=Bfiltr()
        self.Bfilter.min_Bmn  = min_Bmn

        #Use a RH coordinate system (u,w,s) with w=-v. Here, (u,v) are the VMEC coordinates

        if max_m==np.inf:
          Ntheta = 1+2*max(abs(wout.xm))
          self.Bfilter.max_m    = (Ntheta-1)//2
        else:
          Ntheta = max_m*2+1
          self.Bfilter.max_m    = max_m
          
        if maxabs_n==np.inf:
          Nzeta=1+2*max(abs(wout.xn))
          self.Bfilter.maxabs_n = (Nzeta-1)//2
        else:
          Nzeta = maxabs_n*2+1
          self.Bfilter.maxabs_n = maxabs_n

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

        print 'Converting VMEC to Boozer coordinates...'
        for rind in range(len(self.s)):
            print '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bRadius%5i/%5i'% (rind+1,len(self.s)),
            Booz=fluxcoorddiscr.fluxcoorddiscr(wout,rind=rind,Ntheta=Ntheta,Nzeta=Nzeta,name='Boozer')

            self.B00[rind]=Booz.B00
            self.R00[rind]=Booz.R00
            if self.StelSym:
                self.nmodes[rind]=int((Ntheta*Nzeta+1)//2)
            else:
                self.nmodes[rind]=int(Ntheta*Nzeta)
            self.dVdsoverNper[rind]=4*np.pi**2/self.Nperiods*abs(
                self.psi_a*(self.Bphi[rind]+self.iota[rind]*self.Btheta[rind]))/Booz.FSAB2

            FSAB2[rind]=Booz.FSAB2
            self.m.append(np.array([]))
            self.n.append(np.array([]))
            self.B.append(np.array([]))
            self.Bnorm.append(np.array([]))
            self.R.append(np.array([]))
            self.Z.append(np.array([0]))    #m=n=0 cos element is 0
            self.Dphi.append(np.array([0])) #m=n=0 cos element is 0
            if not(self.StelSym):
                self.parity.append(np.array([]))
            
            Blist=(mnmat.mnmat(Booz.B,Nperiods=self.Nperiods)).mnlist()
            Rlist=(mnmat.mnmat(Booz.R,Nperiods=self.Nperiods)).mnlist()
            Zlist=(mnmat.mnmat(Booz.Z,Nperiods=self.Nperiods)).mnlist()
            Dphilist=(mnmat.mnmat(Booz.Dzetapzeta,Nperiods=self.Nperiods)).mnlist()

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

        #end radius loop
        print '' #carriage return

        self.dVdsoverNper=np.abs(self.psi_a*4*np.pi**2.0/self.Nperiods*
                                 (self.Bphi+self.iota*self.Btheta)/self.FSAB2)

        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #% Calculate minorradiusW7AS
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        accum=0;
        for m in range(1,wout.xm[-1]+1):
            ii=[i for i,x in numerate(wout.xm) if x==m]
            first_mode=ii[0]
            last_mode =ii[-1]
            ns=np.expand_dims(wout.xn[ii[0]:ii[-1]]/wout.nfp, axis=1) #row vector
            nnmat=(1.0+(-1.0)**ns-ns.T)/2.0
            accum=accum+m*np.sum((np.expand_dims(wout.rmnc[-1][ii[0]:ii[-1]], axis=1)*
                                  np.expand_dims(wout.zmns[-1][ii[0]:ii[-1]], axis=0))*nnmat)

        self.minorradiusW7AS=np.sqrt(np.abs(accum))
        
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Calculate majorradiusLastbcR00
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Joachim defined his major radius to be
        # R00 in Boozer coordinates at s=0.995
        # which in his case is at the last half mesh radius
        # I choose to always take it at the last half mesh radius
        self.majorradiusLastbcR00=self.R00[-1]

  ####################################################################
  # filter a bcgeom
  ####################################################################
  def filter(self,varinput,max_m=np.inf,maxabs_n=np.inf,makecopy=True):
      #This can be called in the following ways:
      #out=self.filter(min_Bmn)
      #out=self.filter(min_Bmn,max_m,maxabs_n)
      #out=self.filter(Bfilter)
      # if makecopy=True the output is a copy of the input

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
              out.m[rind]=self.m[rind][inds]
              out.n[rind]=self.n[rind][inds]
              out.B[rind]=self.B[rind][inds]
              out.R[rind]=self.R[rind][inds]
              out.Z[rind]=self.Z[rind][inds]
              out.Bnorm[rind]=self.Bnorm[rind][inds]
              out.Dphi[rind]=self.Dphi[rind][inds]
              if not(self.StelSym):
                  out.parity[rind]=self.parity[rind][inds]

              #print self.Bnorm[rind]
              #print out.Bnorm[rind]
              #print out.m[rind]
              #print out.n[rind]
              #sys.exit('stop')
          else:
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
      
  ####################################################################
  # Write a bcgeom to a file
  ####################################################################
  def write(self,filename,Nradii=None,min_Bmn=0,nsortstyle='ascend',printheadercomment=True):
      #argument Nradii is only used for .dat files.
      #First check which type of file the user wants.
      if filename[-3:]=='.bc':
          filetype='bc'
      elif filename[-4:]=='.dat':
          filetype='dat'
          if Nradii is None:
              Nradii=7 #7 is the default for HM's DKES runs
      else: #if no file extension is given, assume .bc
          filetype='bc' 
          filename=filename+'.bc'

      if filetype=='bc':
          if min_Bmn==0:
              selfie=copy.copy(self)
          else:
              selfie=self.filter(min_Bmn)
      else: #.dat
          if Nradii != self.nsurf:
              #HS's choice:
              #rnorms=np.linspace(self.rnorm[0],self.rnorm[-1],num=Nradii+1,endpoint=False)
              #HM's choice:
              rnorms=np.linspace(0.0,1.0,Nradii,endpoint=False)+1.0/Nradii/2.0
              rnorms[ 0] = max(rnorms[ 0],self.rnorm[ 0])
              rnorms[-1] = min(rnorms[-1],self.rnorm[-1])
              selfie=self.interp(rnorms**2.0,'s','linear').filter(min_Bmn)
          else:
              selfie=self.interp(self.s,'s','linear').filter(min_Bmn) #This is just to retrieve dBdrnorm
              
      f = open(filename, 'w')
      now = datetime.datetime.now()
      
      signchange=-1; #sign changer
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

      iota=selfie.iota*signchange
      dVdsoverNper=selfie.dVdsoverNper*signchange
      Dphi=[]
      m=selfie.m
      n=[]
      for tmpind in range(len(selfie.n)):
          Dphi.append(np.array([]))
          n.append(np.array([]))
          Dphi[tmpind]=-selfie.Dphi[tmpind] 
          n[tmpind]=np.where(selfie.m[tmpind]==0,np.abs(selfie.n[tmpind]),-selfie.n[tmpind])

      #################################################################################
      if filetype=='bc':
      #################################################################################
         if printheadercomment:
            f.write('CC Saved %4d.%2d.%2d %2d:%2d by HS python routine bcgeom.write.\n' %
                    (now.year,now.month,now.day,now.hour,now.minute))
            f.write('CC ---------------------------------------------------------\n')
            f.write(selfie.headertext.maincomment+'\n')

         if selfie.StelSym: #JG includes a little more info in this case
            if not(np.isnan(selfie.minorradiusVMEC)) and np.isnan(selfie.majorradiusVMEC) and np.isnan(selfie.volumeVMEC):
                f.write(' m0b  n0b nsurf nper  flux/[Tm^2]     a/[m]     R/[m]   avol/[m]\n')
                f.write('%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f\n' %
                        (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                         selfie.minorradiusW7AS,selfie.majorradiusLastbcR00,selfie.minorradiusVMEC))
            elif not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC)) and np.isnan(selfie.volumeVMEC):
                f.write(' m0b  n0b nsurf nper  flux/[Tm^2]     a/[m]     R/[m]   avol/[m]  '+
                        'Rvol/[m]\n')
                f.write('%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n' %
                        (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                         selfie.minorradiusW7AS,selfie.majorradiusLastbcR00,
                         selfie.minorradiusVMEC,selfie.majorradiusVMEC,
                         np.pi*selfie.minorradiusVMEC**2.0*2.0*np.pi*selfie.majorradiusVMEC))
            elif not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC)) and not(np.isnan(selfie.volumeVMEC)):
                f.write(' m0b  n0b nsurf nper  flux/[Tm^2]     a/[m]     R/[m]   avol/[m]  '+
                        'Rvol/[m] Vol/[m^3]\n')
                f.write('%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n' %
                        (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                         selfie.minorradiusW7AS,selfie.majorradiusLastbcR00,
                         selfie.minorradiusVMEC,selfie.majorradiusVMEC,selfie.volumeVMEC))
            else: #all three are nan
                f.write(' m0b  n0b nsurf nper flux/[Tm^2]     a/[m]     R/[m]\n')
                f.write('%3d%5d%6d%4d%16.6E%10.5f%10.5f\n' %
                        (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                         selfie.minorradiusW7AS,selfie.majorradiusLastbcR00))
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
              do_sort=False #If False, I choose not to mn-sort the components because they usually come sorted.
              if do_sort:
                  inds=np.lexsort((n[rind],m[rind]))
                  for j in range(0,selfie.nmodes[rind]):
                      rmnc=selfie.R[rind][inds[j]]
                      zmns=selfie.Z[rind][inds[j]]
                      vmns=Dphi[rind][inds[j]]
                      bmnc=selfie.B[rind][inds[j]]
                      f.write('%5d%5d%16.8E%16.8E%16.8E%16.8E\n' %
                              (m[rind][inds[j]],n[rind][inds[j]],
                               selfie.R[rind][inds[j]],selfie.Z[rind][inds[j]],
                               Dphi[rind][inds[j]],selfie.B[rind][inds[j]]))
              else:
                   for ind in range(0,selfie.nmodes[rind]):
                       f.write('%5d%5d%16.8E%16.8E%16.8E%16.8E\n' %
                               (m[rind][ind],n[rind][ind],
                                selfie.R[rind][ind],selfie.Z[rind][ind],
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
                  #print thism, thisn, j, selfie.nmodes[rind], selfie.parity[rind][inds[j]]
                  if selfie.parity[rind][inds[j]]==1: #j is cosinus components
                       rmnc=selfie.R[rind][inds[j]]
                       zmns=selfie.Z[rind][inds[j]]
                       vmns=Dphi[rind][inds[j]]
                       bmnc=selfie.B[rind][inds[j]]
                       j+=1              
                  else: #j is sinus components
                      rmns=selfie.R[rind][inds[j]]
                      zmnc=selfie.Z[rind][inds[j]]
                      vmnc=Dphi[rind][inds[j]]
                      bmns=selfie.B[rind][inds[j]]
                      if j+1==selfie.nmodes[rind]:
                          #print 'A '+str(j)
                          j+=1
                      else: 
                          if thism==m[rind][inds[j+1]] and thisn==n[rind][inds[j+1]] and selfie.parity[rind][inds[j+1]]==1:
                              #print 'B1 '+str(j)
                              rmnc=selfie.R[rind][inds[j+1]]
                              zmns=selfie.Z[rind][inds[j+1]]
                              vmns=Dphi[rind][inds[j+1]]
                              bmnc=selfie.B[rind][inds[j+1]]
                              j+=2
                          else:
                              #print 'B2 '+str(j)
                              j+=1 
                  f.write('%5d%5d%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E\n' %
                          (thism,thisn,rmnc,rmns,zmnc,zmns,vmnc,vmns,bmnc,bmns))

      #################################################################################
      else: #save a .dat file instead
      #################################################################################
         if np.isnan(selfie.majorradiusLastbcR00):
             majorradiusLastbcR00=0.0
         else:
             majorradiusLastbcR00=selfie.majorradiusLastbcR00

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
                            selfie.minorradiusW7AS,selfie.majorradiusLastbcR00,selfie.minorradiusVMEC))
                elif (not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC))
                      and np.isnan(selfie.volumeVMEC)):
                   f.write('cc %3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n'%
                           (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                            selfie.minorradiusW7AS,selfie.majorradiusLastbcR00,
                            selfie.minorradiusVMEC,selfie.majorradiusVMEC,
                            np.pi*selfie.minorradiusVMEC**2.0*2.0*np.pi*selfie.majorradiusVMEC))
                elif (not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC))
                      and not(np.isnan(selfie.volumeVMEC))):
                   f.write('cc %3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n'%
                           (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                            selfie.minorradiusW7AS,selfie.majorradiusLastbcR00,
                            selfie.minorradiusVMEC,selfie.majorradiusVMEC,selfie.VolumeVMEC))
                else:
                   f.write('cc %3d%5d%6d%4d%16.6E%10.5f%10.5f\n'%
                           (selfie.m0b,selfie.n0b,selfie.nsurf,selfie.Nperiods,torfluxtot,
                            selfie.minorradiusW7AS,selfie.majorradiusLastbcR00))
             else:
                if (not(np.isnan(selfie.minorradiusVMEC)) and np.isnan(selfie.majorradiusVMEC)
                    and np.isnan(selfie.volumeVMEC)):
                   f.write('cc %6d%4d%10.5f%10.5f%11.5f\n'%
                           (selfie.nsurf,selfie.Nperiods,
                            selfie.minorradiusW7AS,selfie.majorradiusLastbcR00,selfie.minorradiusVMEC))
                elif (not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC))
                      and np.isnan(selfie.volumeVMEC)):
                   f.write('cc %6d%4d%10.5f%10.5f%11.5f%11.5f%11.5f\n'%
                           (selfie.nsurf,selfie.Nperiods,
                            selfie.minorradiusW7AS,selfie.majorradiusLastbcR00,
                            selfie.minorradiusVMEC,selfie.majorradiusVMEC,
                            np.pi*selfie.minorradiusVMEC**2.0*2.0*np.pi*selfie.majorradiusVMEC))
                elif (not(np.isnan(selfie.minorradiusVMEC)) and not(np.isnan(selfie.majorradiusVMEC))
                      and not(np.isnan(selfie.volumeVMEC))):
                   f.write(' cc %6d%4d%10.5f%10.5f%11.5f%11.5f%11.5f\n'%
                           (selfie.nsurf,selfie.Nperiods,
                            selfie.minorradiusW7AS,selfie.majorradiusLastbcR00,
                            selfie.minorradiusVMEC,selfie.majorradiusVMEC,selfie.VolumeVMEC))
                else:
                   f.write('cc %6d%4d%10.5f%10.5f\n'%
                           (selfie.nsurf,selfie.Nperiods,
                            selfie.minorradiusW7AS,selfie.majorradiusLastbcR00))
         #end if headercomment
         
         f.write(' c     m, n, B_mn/B_00, B_mn, (dB_mn/dr)/B_00 [1/cm], R_mn, z_mn\n')

         for rind in range(selfie.nsurf):
             f.write('%6.2f%8.4f%3i%8.2f%8.2f%7.4f%8.2f%8.2f  r,io,N,ra,Rm,s,r_o,R00\n'%
                     (selfie.minorradiusW7AS*selfie.rnorm[rind]*100.0,
                      iota[rind],selfie.Nperiods,selfie.minorradiusW7AS*100.0,
                      selfie.majorradiusLastbcR00*100.0,selfie.s[rind],
                      selfie.minorradiusW7AS*selfie.rnorm[rind]*100.0,selfie.R00[rind]*100.0))
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
                          selfie.dBdrnorm[rind][m0i]/selfie.B00[rind]/selfie.minorradiusW7AS/100.0,
                          selfie.R[rind][m0i]*100.0,
                          selfie.Z[rind][m0i]*100.0))
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
                              selfie.dBdrnorm[rind][mi]/selfie.B00[rind]/selfie.minorradiusW7AS/100.0,
                              selfie.R[rind][mi]*100.0,
                              selfie.Z[rind][mi]*100.0))
             #end for thism
             f.write('  -1   0 0.0 0.0 0.0 0.0 0.0     end of input\n')
         #end for rind
         
      #end if filetype=='bc' or 'dat'   
      f.close()

  ####################################################################
  # Interpolate a bcgeom to a set of new radii
  ####################################################################
  def interp(self,invar,invarname='rnorm',interptype='linear'):
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
          if np.any(s<smin) or np.any(s>smax):
            sys.exit('Chosen s out of range!')
      elif invarname=='rnorm':
          s=invar**2
          if np.any(s<smin) or np.any(s>smax):
            sys.exit('Chosen rnorm out of range!')
      elif invarname=='rW7AS':
          s=(invar/self.minorradiusW7AS)**2
          if np.any(s<smin) or np.any(s>smax):
            sys.exit('Chosen rW7AS out of range!')
      elif invarname=='rVMEC':
          s=(invar/self.minorradiusVMEC)**2
          if np.any(s<smin) or np.any(s>smax):
            sys.exit('Chosen rVMEC out of range!')
      elif invarname=='index':
          if np.any(invar<0) or np.any(invar>(len(self.s)-1)):
            sys.exit('Chosen index out of range!')
          invar=np.sort(invar)
          s=self.s(invar)
          interptype='nearest'
      else:
          sys.exit('invarname not recognised!')
        
      out=copy.deepcopy(self)

      s=np.sort(s)
      if np.any(np.diff(s)<=0):
          sys.exit('invar contains doublets')
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
              #print '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%4i/%4i'% (surfind+1,out.nsurf)
              if not(s[surfind] in self.s):
                  Rind=np.where(np.sign(self.s-s[surfind])==1)[0][0]
              else:
                  exactind=np.where(self.s==s[surfind])[0][0]
                  Rind=exactind+1
                  if exactind==len(s):
                      Rind=len(s)-1
                      
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
          #print '\n'
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
        out.nmodes = self.nmodes[inds]
        out.B      = self.B[inds]
        out.R      = self.R[inds]
        out.Z      = self.Z[inds]
        out.Dphi   = self.Dphi[inds]
        out.Bnorm  = self.Bnorm[inds]
        out.m      = self.m[inds]
        out.n      = self.n[inds]
        if not(self.StelSym):
           out.parity = self.parity[inds]

      else:
        sys.exit('Unknown interpolation type!')

      return out

        
  ####################################################################
  # Display the content of a bcgeom
  ####################################################################
  
  def disp(self,verbose=0):

    print '-----------------------------------------------'
    print 'Comment lines: headertext.maincomment='
    #print '-----------------------------------------------'

    print self.headertext.maincomment
    print '-----------------------------------------------'
    print 'Scalars:'
    #print '-----------------------------------------------'
    print 'nsurf                = '+str(self.nsurf)
    print 'Nperiods             = '+str(self.Nperiods)
    print 'psi_a                = '+str(self.psi_a)
    print 'torfluxtot           = '+str(self.torfluxtot)
    print 'minorradiusW7AS      = '+str(self.minorradiusW7AS)
    print 'majorradiusLastbcR00 = '+str(self.majorradiusLastbcR00)
    print 'minorradiusVMEC      = '+str(self.minorradiusVMEC)
    print 'majorradiusVMEC      = '+str(self.majorradiusVMEC)
    print 'volumeVMEC           = '+str(self.volumeVMEC)
    print 'Bfilter.min_Bmn      = '+str(self.Bfilter.min_Bmn)
    print 'Bfilter.max_m        = '+str(self.Bfilter.max_m)
    print 'Bfilter.maxabs_n     = '+str(self.Bfilter.maxabs_n)

    print '-----------------------------------------------'
    print 'Radial arrays:'
    #print '-----------------------------------------------'
    if verbose==1:
      print 'rnorm='+str(self.rnorm)
      print 's='+str(self.s)
      print 'Bphi='+str(self.Bphi)
      print 'Btheta='+str(self.Btheta)
      print 'iota='+str(self.iota)
      print 'dpds='+str(self.dpds)
  
      print 'dVdsoverNper='+str(self.dVdsoverNper)
      print 'FSAB2='+str(self.FSAB2)
      print 'nmodes='+str(self.nmodes)
      print 'B00='+str(self.B00)
      print 'R00='+str(self.R00)
    elif verbose==0:
      print 'rnorm, s, Bphi, Btheta, iota, dpds, dVdsoverNper, '
      print 'FSAB2, nmodes, B00, R00'

    print '-----------------------------------------------'
    if self.StelSym:
       print 'Fourier indexes:    m, n'
    else:
       print 'Fourier indexes:    m, n, parity'

    print 'Fourier quantities: B, R, Z, Dphi'
    print '-----------------------------------------------'




    

    
 







      
   
