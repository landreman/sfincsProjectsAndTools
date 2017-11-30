#!/usr/bin/env python
import numpy as np
import sys
#from mnlist import mnlist
#from mnmat import mnmat

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
       # signcorr=2: Total toroidal flux will get a sign change (default, as YT)

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
               if tmp_str[0:2]=='CC':
                   concat_str = concat_str+'\n'+tmp_str #comment line
               if 'CStconfig' in tmp_str:
                   YTsign=-1

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
         self.nsurf      = header_df[2]
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
           self.VolumeVMEC=header_df[9]
         else:
           self.VolumeVMEC=np.nan
         
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
         no_of_modes=np.array([])
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
           Bphi=np.append(Bphi,surfheader[2]*self.Nperiods/2/np.pi*(4*np.pi*1e-7)) #Tesla*meter
           Btheta=np.append(Btheta,surfheader[3]/2/np.pi*(4*np.pi*1e-7))              #Tesla*meter
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
                   proceed=False
                   endoffile=True
                   #print 'found end of file'
               elif ('s' in tmp_str): #Next flux surface has been reached
                   proceed=False
                   #print 'found next surface'
               else:
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
           no_of_modes=np.append(no_of_modes,modeind)
           if no_of_modes[rind]==0:
               sys.exit('no modes found for rind='+str(rind))

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
         self.psi_a=self.torfluxtot/2/np.pi
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
         self.m=modesm
         self.n=modesn #sign is switched below
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
           self.VolumeVMEC=pi*self.minorradiusVMEC**2*2*pi*self.majorradiusVMEC
           self.minorradiusW7AS=NaN
           self.majorradiusLastbcR00=NaN


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

           endoffile=False
           rind=-1

           while not endoffile:
               rind=rind+1
               surfheader1=sscan(tmp_str,count=8)
               radii[rind] = surfheader1[0]/100 #(cm)
               iota[rind]  = surfheader1[1]
               Nperiods=surfheader1[2]
               minorradiusW7AS = surfheader1[3]/100 #(cm)
               majorradiusLastbcR00  = surfheader1[4]/100 #(cm)
               if len(surfheader1)>5:
                   torfluxnorm[rind] = surfheader1[5]
                   R00[rind]   = surfheader1[7]/100
               else:
                   torfluxnorm[rind] = np.nan
                   R00[rind] =np.nan
       
               tmp_str=f.readline()
               if tmp_str[1]=='>':
                   surfheader2=sscan(tmp_str[2:],count=3)
           
                   Bphi[rind]  =surfheader2[0]*1e6*Nperiods/2/pi*(4*pi*1e-7) #Tesla*meter
                   Btheta[rind]=surfheader2[1]*1e6/2/pi*(4*pi*1e-7)          #Tesla*meter
                   B00[rind]   =surfheader2[2] #Tesla
                   tmp_str=f.readline()
               else:
                   Bphi[rind]  =np.nan
                   Btheta[rind]=np.nan
                   B00[rind]   =np.nan

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
                           modesm[rind][modeind]=tmp[0]
                           modesn[rind][modeind]=tmp[1]
                           modesbnorm[rind][modeind]=tmp[2]
                           modesb[rind][modeind]=tmp[3]
                           modesdbdrnorm[rind][modeind]=tmp[4]*100 #convert cm^-1 -> m^-1
                           modesr[rind][modeind]=tmp[5]
                           modesz[rind][modeind]=tmp[6]
         
                   tmp_str=f.readline() #get the next line or surface header line
               #end while proceed
               no_of_modes[rind]=modeind  

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
           self.rnorm=radii/self.minorradiusW7AS
           self.s=torfluxnorm
           self.iota=iota*rthetazeta_righthanded
           self.Bphi=Bphi
           self.Btheta=Btheta*rthetazeta_righthanded
           self.nmodes=no_of_modes
           self.m=modesm
           self.n=modesn #sign is switched below
           self.B=modesb
           self.Bnorm=modesbnorm
           self.B00=B00
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
    elif isinstance(input,vmecgeom):
        wout=input
        signchange=float(wout.signgs) #is -1, because vmec is left handed
        self.headertext=headertxt()
        self.headertext.maincomment=('CC Converted from VMEC using HSs routines\n'+
        'CC version_       = '+str(self.version_)+'\n'+
        'CC input_extension= '+self.input_extension)
        #'CC The phase convention (m theta - n phi) is used in this file.\n'+
        #'CC The coordinate system (r,theta,phi) is left-handed.\n'+
        #'CC Toroidal flux is counted positive in the direction of the toroidal coordinate.'

        self.m0b=wout.mpol
        self.n0b=wout.ntor
        self.nsurf= np.nan      #number of radial surfaces, Set this below
        self.Nperiods=wout.ns   #%!< number of field periods
        self.torfluxtot = wout.phi[wout.ns-1]*signchange
        self.psi_a=Geom.torfluxtot/2/np.pi
        self.minorradiusVMEC      = wout.Aminor_p  #minor plasma radius
        self.majorradiusLastbcR00 = np.nan         #Calculate this below (not necessary)
        self.minorradiusW7AS      = np.nan         #Calculate this below (not necessary)
        self.majorradiusVMEC      = wout.Rmajor_p  #major plasma radius
        fullgrid_s           = wout.phi/wout.phi[wout.ns-1] #full grid
        fullgrid_rnorm       = sqrt(fullgrid_s);     #full grid

        skip=wout.skip #=1,this is how many elements are skipped at low radii when going to half grid

        self.s      = fullgrid_s(skip-1:-2)+fullgrid_s(skip:))/2  #half grid
        self.rnorm  = sqrt(self.s) #half grid
        self.nsurf  = len(self.s)
        self.dpds   = np.array(np.diff(wout.presf[skip-1:])./np.diff(fullgrid_s[skip-1:])
        self.Bphi        = wout.bvco[wout.skip:]*signchange #direction switch sign
        self.Btheta      = wout.buco[wout.skip:] #sign change

        fullgrid_iota=wout.iotaf*signchange
        iota   = wout.iotas*signchange #half mesh
        self.iota=iota[skip:]
        
        self.Bfilter.min_Bmn  = min_Bmn

        #Use a RH coordinate system (u,w,s) with w=-v. Here, (u,v) are the VMEC coordinates

        if max_m==np.inf:
          Ntheta = 1+2*max(abs(float(wout.xm)))
          self.Bfilter.max_m    = (Ntheta-1)/2
        else:
          Ntheta = max_m*2+1
          self.Bfilter.max_m    = max_m
          
        if maxabs_n==np.inf:
          Nzeta=1+2*max(abs(wout.xn))
          self.Bfilter.maxabs_n = (Nzeta-1)/2
        else:
          Nzeta = maxabs_n*2+1
          self.Bfilter.maxabs_n = maxabs_n

        dpds=np.array([])
        dVdsoverNper=np.array([])
        B00=np.array([])
        R00=np.array([])
        no_of_modes=np.array([])
        modesm=[]
        modesn=[]
        modesr=[]
        modesz=[]
        modesp=[]
        modesb=[]
        modesbnorm=[]
        modespar=[]
         
        for rind in range(len(self.s)):
            Booz=fluxcoorddiscr(wout,rind=rind,Ntheta=Ntheta,Nzeta=Nzeta,name='Boozer')
        
        self.nmodes=no_of_modes
        self.m=modesm
        self.n=modesn #sign is switched below
        self.B=modesb
        self.Bnorm=modesbnorm
        self.B00=B00
        self.R00=R00
        self.R=modesr
        self.Z=modesz
          
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
    print 'VolumeVMEC           = '+str(self.VolumeVMEC)
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




    

    
 







      
   
