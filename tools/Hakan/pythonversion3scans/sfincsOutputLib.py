#!/usr/bin/env python
import numpy as np 
import os, sys, inspect, math, h5py, copy
import subprocess
import matplotlib.pyplot as plt

def inp(promptstr):
        if sys.version_info[0] > 2:
            return input(promptstr)
        else:
            return raw_input(promptstr)

#################################################################################################################
class sfincsScan:
#################################################################################################################

  def initiate_variables(self,Nruns,Nspecies):
      #Quantities that should be the same for all runs
      self.RHSMode               = 1
      self.NPeriods              = None
      self.psiAHat               = None
      self.Nspecies              = Nspecies
      self.Zs                    = np.nan*np.zeros((Nspecies))
      self.mHats                 = np.nan*np.zeros((Nspecies))
      self.includePhi1           = None
      self.withAdiabatic         = None
      self.adiabaticZ            = None
      self.adiabaticMHat         = None
      self.withNBIspec           = None
      self.NBIspecZ              = None

      #initiate a lot
      self.finished                        = [None]*Nruns
      self.didNonlinearCalculationConverge = [None]*Nruns 
      self.Ntheta   = np.zeros((Nruns))
      self.Nzeta    = np.zeros((Nruns))
      self.Nxi      = np.zeros((Nruns))
      self.Nx       = np.zeros((Nruns))
      self.NL       = np.zeros((Nruns))
      self.solverTolerance = np.zeros((Nruns))
      self.theta    = [None]*Nruns
      self.zeta     = [None]*Nruns
      self.psiHat   = np.zeros((Nruns))
      self.rN       = np.zeros((Nruns))
      self.rHat     = np.zeros((Nruns))
      self.GHat     = np.zeros((Nruns))
      self.IHat     = np.zeros((Nruns))
      self.B0OverBBar  = np.zeros((Nruns))
      self.iota     = np.zeros((Nruns))
      self.VPrimeHat= np.zeros((Nruns))
      self.FSABHat2 = np.zeros((Nruns))
      self.alpha    = np.zeros((Nruns))
      self.Delta    = np.zeros((Nruns))
      self.nu_n     = np.zeros((Nruns))
      self.BHat                     = [None]*Nruns
      self.dBHatdtheta              = [None]*Nruns
      self.dBHatdzeta               = [None]*Nruns
      self.BHat_sub_psi             = [None]*Nruns
      self.BHat_sup_theta           = [None]*Nruns
      self.BHat_sup_zeta            = [None]*Nruns
      self.dBHat_sub_psi_dtheta     = [None]*Nruns
      self.dBHat_sub_psi_dzeta      = [None]*Nruns
      self.dIHat_dpsiHat            = [None]*Nruns
      self.dGHat_dpsiHat            = [None]*Nruns
      self.dBHat_sup_theta_dpsiHat  = [None]*Nruns
      self.dBHat_sup_theta_dzeta    = [None]*Nruns
      self.dBHat_sup_zeta_dpsiHat   = [None]*Nruns
      self.dBHat_sup_zeta_dtheta    = [None]*Nruns
      self.dBHatdpsiHat             = [None]*Nruns
      self.gpsiHatpsiHat            = [None]*Nruns

      self.nHats         = np.zeros((Nruns,Nspecies))
      self.THats         = np.zeros((Nruns,Nspecies))
      self.dnHatdpsiN    = np.zeros((Nruns,Nspecies))
      self.dnHatdrN      = np.zeros((Nruns,Nspecies))
      self.dnHatdrHat    = np.zeros((Nruns,Nspecies))
      self.dTHatdpsiN    = np.zeros((Nruns,Nspecies))
      self.dTHatdrN      = np.zeros((Nruns,Nspecies))
      self.dTHatdrHat    = np.zeros((Nruns,Nspecies))
      self.dPhiHatdpsiN  = np.zeros((Nruns))
      self.dPhiHatdrN    = np.zeros((Nruns))
      self.dPhiHatdrHat  = np.zeros((Nruns))
      self.EParallelHat  = np.zeros((Nruns))
      self.adiabaticNHat = np.nan*np.zeros((Nruns))
      self.adiabaticTHat = np.nan*np.zeros((Nruns))
      self.NBIspecNHat   = np.nan*np.zeros((Nruns))

      self.particleFlux_vm_rHat=np.zeros((Nruns,Nspecies))
      self.particleFlux_vm_rN  =np.zeros((Nruns,Nspecies))
      self.particleFlux_vm_psiN=np.zeros((Nruns,Nspecies))
      self.particleFlux_vd_rHat=np.zeros((Nruns,Nspecies))
      self.particleFlux_vd_rN  =np.zeros((Nruns,Nspecies))
      self.particleFlux_vd_psiN=np.zeros((Nruns,Nspecies))
      self.heatFlux_vm_rHat    =np.zeros((Nruns,Nspecies))
      self.heatFlux_vm_rN      =np.zeros((Nruns,Nspecies))
      self.heatFlux_vm_psiN    =np.zeros((Nruns,Nspecies))
      self.heatFlux_vd_rHat    =np.zeros((Nruns,Nspecies))
      self.heatFlux_vd_rN      =np.zeros((Nruns,Nspecies))
      self.heatFlux_vd_psiN    =np.zeros((Nruns,Nspecies))
      self.momentumFlux_vm_rHat=np.zeros((Nruns,Nspecies))
      self.momentumFlux_vm_rN  =np.zeros((Nruns,Nspecies))
      self.momentumFlux_vm_psiN=np.zeros((Nruns,Nspecies))
      self.momentumFlux_vd_rHat=np.zeros((Nruns,Nspecies))
      self.momentumFlux_vd_rN  =np.zeros((Nruns,Nspecies))
      self.momentumFlux_vd_psiN=np.zeros((Nruns,Nspecies))
      self.FSABFlow            =np.zeros((Nruns,Nspecies))
      self.FSABjHat            =np.zeros((Nruns))
      self.NTV                 =np.zeros((Nruns,Nspecies))

      self.classicalParticleFlux_rHat=np.zeros((Nruns,Nspecies))
      self.classicalParticleFlux_rN  =np.zeros((Nruns,Nspecies))
      self.classicalParticleFlux_psiN=np.zeros((Nruns,Nspecies))

      #These may not always need to be loaded:
      self.flow                     = [None]*Nruns
      self.densityPerturbation      = [None]*Nruns
      self.pressurePerturbation     = [None]*Nruns
      self.pressureAnisotropy       = [None]*Nruns
      self.Phi1Hat                  = [None]*Nruns
      self.NTVBeforeSurfaceIntegral = [None]*Nruns

  def __init__(self,mainDirectory,sortafter='rN',verbose=0,collapseErScans=False):
    ########################################################
    # Begin __init__()
    ########################################################
    if not(collapseErScans): #This is the normal case, loading one scan
      #print('normal case: collapseErScans=False')
      if mainDirectory is None:
        print('Missing input mainDirectory!')
        sys.exit('Missing input mainDirectory!')
      if mainDirectory[-1]=='/':
         mainDirectory=mainDirectory[:-1]
      self.mainDir=mainDirectory
      if verbose>0:
        print("Extracting data from radial scan in " + mainDirectory)
      # Get a list of the subdirectories:
      dirList=os.listdir(mainDirectory)
      CandidateDirs=[]
      for ind in range(len(dirList)):
        if os.path.isdir(mainDirectory+'/'+dirList[ind]):
          CandidateDirs.append(dirList[ind])

      CandidateDirs = sorted(CandidateDirs)
      if len(CandidateDirs) < 1:
        print('Error! Could not find any directories in ' + mainDirectory)
        sys.exit(1)

      #Check for valid directories
      DataDirs=[]
      MissDirs=[]
      for dirind in range(len(CandidateDirs)):
        if os.path.isfile(mainDirectory+'/'+CandidateDirs[dirind]+'/sfincsOutput.h5'):
          DataDirs.append(CandidateDirs[dirind])
        else:
          MissDirs.append(CandidateDirs[dirind])
      if len(DataDirs)<1:
        print('Could not find any files sfincsOutput.h5 in the directories!')
        sys.exit('Could not find any files sfincsOutput.h5 in the directories!')

      #sort directories after rN
      sortQuant                              = np.zeros((len(DataDirs)))
      for dirind in range(len(DataDirs)):
        file = h5py.File(mainDirectory + '/' + DataDirs[dirind] + '/sfincsOutput.h5','r')
        sortQuant[dirind]=file[sortafter][()]
        Nspecies=file['Nspecies'][()]
        file.close

      sortind=np.argsort(sortQuant)
      self.DataDirs=[None]*len(DataDirs)
      for ind in range(len(DataDirs)):
        self.DataDirs[ind]=DataDirs[sortind[ind]]

      self.sortafter=sortafter
      self.Nspecies=Nspecies
      self.Nruns = len(DataDirs)
      self.initiate_variables(self.Nruns,self.Nspecies)

      for ind in range(len(self.DataDirs)):
        fullDirectory = mainDirectory + "/" + self.DataDirs[ind]

        if verbose>0:
          print("*************************************************")
          print("Processing directory "+self.DataDirs[ind])
          print("*************************************************")

        
        file = h5py.File(fullDirectory + '/sfincsOutput.h5','r') 

        if file['RHSMode'][()]!=1:
          print('RHSMode was not = 1. Not implemented yet!')
          sys.exit('RHSMode was not = 1. Not implemented yet!')
        
        integerToRepresentTrue = file['integerToRepresentTrue'][()]
        if 'finished' in file:
          self.finished[ind] = (file['finished'][()]== integerToRepresentTrue)
        else:
          self.finished[ind] = False 
        #Quantities that should be the same for all runs
        NPeriods     = file['NPeriods'][()]
        psiAHat      = file['psiAHat'][()]
        Zs           = file['Zs'][()]
        mHats        = file['mHats'][()]
        includePhi1  = (file['includePhi1'][()]==integerToRepresentTrue)
        withAdiabatic= (file['withAdiabatic'][()]==integerToRepresentTrue)
        withNBIspec  = (file['withNBIspec'][()]==integerToRepresentTrue)
        if withAdiabatic:
          adiabaticZ =(file['adiabaticZ'][()]==integerToRepresentTrue)
          adiabaticMHat  = file['adiabaticMHat'][()]
        if withNBIspec:
          NBIspecZ   = file['NBIspecZ'][()]

        #consistency checks
        if self.NPeriods is None: #First run
          self.NPeriods     = NPeriods
          self.psiAHat      = psiAHat
          self.Zs           = Zs
          self.mHats        = mHats
          self.includePhi1  = includePhi1
          self.withAdiabatic= withAdiabatic
          self.withNBIspec  = withNBIspec
          if withAdiabatic:
            self.adiabaticZ    = adiabaticZ
            self.adiabaticMHat = adiabaticMHat
          if withNBIspec:
            self.NBIspecZ = NBIspecZ
        else:
          if self.NPeriods!=NPeriods:
            sys.exit('Different NPeriods for different runs!')
          if self.psiAHat!=psiAHat:
            sys.exit('Different psiAHat for different runs!')
          if np.any(self.Zs!=Zs):
            sys.exit('Different Zs for different runs!')
          if np.any(self.mHats!=mHats):
            sys.exit('Different mHats for different runs!')
          if self.includePhi1!=includePhi1:
            sys.exit('Different includePhi1 for different runs!')
          if self.withAdiabatic!=withAdiabatic:
            sys.exit('Different withAdiabatic for different runs!')
          if self.withNBIspec!=withNBIspec:
            sys.exit('Different withNBIspec for different runs!')
          if withAdiabatic and self.adiabaticZ!=adiabaticZ:
            sys.exit('Different adiabaticZ for different runs!')
          if withAdiabatic and self.adiabaticMHat!=self.adiabaticMHat:
            sys.exit('Different self.adiabaticMHat for different runs!')
          if withNBIspec and self.NBIspecZ!=NBIspecZ:
            sys.exit('Different NBIspecZ for different runs!')

        self.Ntheta[ind]    = file['Ntheta'][()]
        self.Nzeta[ind]     = file['Nzeta'][()]
        self.Nxi[ind]       = file['Nxi'][()]
        self.Nx[ind]        = file['Nx'][()]
        self.NL[ind]        = file['NL'][()]
        self.solverTolerance[ind] = file['solverTolerance'][()]
        self.theta[ind]     = file['theta'][()]
        self.zeta[ind]      = file['zeta'][()]
        self.psiHat[ind]    = file['psiHat'][()]
        self.rN[ind]        = file['rN'][()]  
        self.rHat[ind]      = file['rHat'][()]  
        self.GHat[ind]      = file['GHat'][()]
        self.IHat[ind]      = file['IHat'][()]
        self.B0OverBBar[ind]= file['B0OverBBar'][()] 
        self.iota[ind]      = file['iota'][()]
        self.VPrimeHat[ind] = file['VPrimeHat'][()]
        self.FSABHat2[ind]  = file['FSABHat2'][()]
        self.alpha[ind]     = file['alpha'][()]
        self.Delta[ind]     = file['Delta'][()]
        self.nu_n[ind]      = file['nu_n'][()]
        self.BHat[ind]                     = file['BHat'][()]
        self.dBHatdtheta[ind]              = file['dBHatdtheta'][()]
        self.dBHatdzeta[ind]               = file['dBHatdzeta'][()]
        self.BHat_sub_psi[ind]             = file['BHat_sub_psi'][()]
        self.BHat_sup_theta[ind]           = file['BHat_sup_theta'][()]
        self.BHat_sup_zeta[ind]            = file['BHat_sup_zeta'][()]
        self.dBHat_sub_psi_dtheta[ind]     = file['dBHat_sub_psi_dtheta'][()]
        self.dBHat_sub_psi_dzeta[ind]      = file['dBHat_sub_psi_dzeta'][()]
        self.dIHat_dpsiHat[ind]            = file['dBHat_sub_theta_dpsiHat'][()]
        self.dGHat_dpsiHat[ind]            = file['dBHat_sub_zeta_dpsiHat'][()]
        self.dBHat_sup_theta_dpsiHat[ind]  = file['dBHat_sup_theta_dpsiHat'][()]
        self.dBHat_sup_theta_dzeta[ind]    = file['dBHat_sup_theta_dzeta'][()]
        self.dBHat_sup_zeta_dpsiHat[ind]   = file['dBHat_sup_zeta_dpsiHat'][()]
        self.dBHat_sup_zeta_dtheta[ind]    = file['dBHat_sup_zeta_dtheta'][()]
        self.dBHatdpsiHat[ind]             = file['dBHatdpsiHat'][()]
        try:
          self.gpsiHatpsiHat[ind]=file['gpsiHatpsiHat'][()]
        except:
          if verbose>0:
            print('old version where gpsiHatpsiHat was not stored')

        self.nHats[ind]        = file['nHats'][()]
        self.THats[ind]        = file['THats'][()]
        self.dnHatdpsiN[ind]   = file['dnHatdpsiN'][()]
        self.dnHatdrN[ind]     = file['dnHatdrN'][()]
        self.dnHatdrHat[ind]   = file['dnHatdrHat'][()]
        self.dTHatdpsiN[ind]   = file['dTHatdpsiN'][()]
        self.dTHatdrN[ind]     = file['dTHatdrN'][()]
        self.dTHatdrHat[ind]   = file['dTHatdrHat'][()]
        self.dPhiHatdpsiN[ind] = file['dPhiHatdpsiN'][()]
        self.dPhiHatdrN[ind]   = file['dPhiHatdrN'][()]
        self.dPhiHatdrHat[ind] = file['dPhiHatdrHat'][()]
        self.EParallelHat[ind] = file['EParallelHat'][()]
        if withAdiabatic:
          self.adiabaticNHat[ind] = file['adiabaticNHat'][()]
          self.adiabaticTHat[ind] = file['adiabaticTHat'][()]
        if withNBIspec:
          self.NBIspecNHat[ind]   = file['NBIspecNHat'][()]
        if self.finished[ind] and 'FSABFlow' in file:
          self.FSABFlow[ind]            =file['FSABFlow'][()][:,-1]
          self.FSABjHat[ind]            =file['FSABjHat'][()][-1]
          self.NTV[ind]                 =file['NTV'][()][:,-1]
          if not(self.includePhi1):
            # print file['particleFlux_vm_rHat'][()][:,-1]
            #print self.particleFlux_vm_rHat[ind]
            self.particleFlux_vm_rHat[ind]=file['particleFlux_vm_rHat'][()][:,-1]
            self.particleFlux_vm_rN[ind]  =file['particleFlux_vm_rN'][()][:,-1]
            self.particleFlux_vm_psiN[ind]=file['particleFlux_vm_psiN'][()][:,-1]
            self.heatFlux_vm_rHat[ind]    =file['heatFlux_vm_rHat'][()][:,-1]
            self.heatFlux_vm_rN[ind]      =file['heatFlux_vm_rN'][()][:,-1]
            self.heatFlux_vm_psiN[ind]    =file['heatFlux_vm_psiN'][()][:,-1]
            self.momentumFlux_vm_rHat[ind]=file['momentumFlux_vm_rHat'][()][:,-1]
            self.momentumFlux_vm_rN[ind]  =file['momentumFlux_vm_rN'][()][:,-1]
            self.momentumFlux_vm_psiN[ind]=file['momentumFlux_vm_psiN'][()][:,-1]
            if 'classicalParticleFlux_rHat' in file:
              self.classicalParticleFlux_rHat[ind]=file['classicalParticleFlux_rHat'][()][:,-1]
              self.classicalParticleFlux_rN[ind]  =file['classicalParticleFlux_rN'][()][:,-1]
              self.classicalParticleFlux_psiN[ind]=file['classicalParticleFlux_psiN'][()][:,-1]
          else:#if self.includePhi1
            self.didNonlinearCalculationConverge[ind] = (file['didNonlinearCalculationConverge'][()]==integerToRepresentTrue)
            self.Phi1Hat[ind]             =file['Phi1Hat'][()]
            self.particleFlux_vd_rHat[ind]=file['particleFlux_vd_rHat'][()][:,-1]
            self.particleFlux_vd_rN[ind]  =file['particleFlux_vd_rN'][()][:,-1]
            self.particleFlux_vd_psiN[ind]=file['particleFlux_vd_psiN'][()][:,-1]
            self.heatFlux_vd_rHat[ind]    =file['heatFlux_vd_rHat'][()][:,-1]
            self.heatFlux_vd_rN[ind]      =file['heatFlux_vd_rN'][()][:,-1]
            self.heatFlux_vd_psiN[ind]    =file['heatFlux_vd_psiN'][()][:,-1]
            self.momentumFlux_vd_rHat[ind]=file['momentumFlux_vd_rHat'][()][:,-1]
            self.momentumFlux_vd_rN[ind]  =file['momentumFlux_vd_rN'][()][:,-1]
            self.momentumFlux_vd_psiN[ind]=file['momentumFlux_vd_psiN'][()][:,-1]
            if 'classicalParticleFlux_rHat' in file:
              self.classicalParticleFlux_rHat[ind]=file['classicalParticleFlux_rHat'][()][:,-1]
              self.classicalParticleFlux_rN[ind]  =file['classicalParticleFlux_rN'][()][:,-1]
              self.classicalParticleFlux_psiN[ind]=file['classicalParticleFlux_psiN'][()][:,-1]
          #end if self.includePhi1
        #end if self.finished[ind]
        file.close
      #end for ind in range(len(self.DataDirs))
    else: #collapseErScans=True 
      #print('abnormal case: collapseErScans=True')
      #print('verbose='+str(verbose))
      # This reduces the 2D parameter scan to 1D. It returns a sfincsScan object like
      # if a simple scan over radius had been made
      RadialAndErScan=sfincsRadialAndErScan(mainDirectory,verbose=verbose)
      Nradii=RadialAndErScan.Nradii
      Nspecies=RadialAndErScan.Erscans[0].Nspecies
      self.initiate_variables(Nradii,Nspecies)
      self.DataDirs=[None]*Nradii

      bestErind=-1*np.ones((Nradii),dtype=int) 
      for i in range(Nradii):
          bestErind[i]=RadialAndErScan.Erscans[i].choose_Erscan_run_with_best_Er()
     
      #print('bestErind='+str(bestErind))
      attrs=dir(RadialAndErScan.Erscans[0])
      for attr in attrs:
        if attr[0]!='_' and not(hasattr(getattr(RadialAndErScan.Erscans[0],attr), '__call__')):
          val=copy.deepcopy(getattr(RadialAndErScan.Erscans[0],attr))
          if not(isinstance(val,(np.ndarray,list))):
            setattr(self,attr,val)
          elif attr=='Zs' or attr=='mHats':
            setattr(self,attr,val)
          else: #if len(val.shape[0])==2: or 1
            tmp=getattr(self,attr) 
            for i in range(Nradii):
              if bestErind[i]==-1:
                tmp[i]=np.nan
              else:
                tmp[i]=copy.deepcopy(getattr(RadialAndErScan.Erscans[i],attr)[bestErind[i]])
            setattr(self,attr,tmp)

      self.sortafter='rN'
      self.Nruns=Nradii
      self.Nspecies=Nspecies
      self.mainDir=mainDirectory
      self.reduced=True


  def disp(self,radialCoord='rHat'):
    print('--------------------------------------------------------------------------------------')
    print('Data from radial scan in ' + self.mainDir)
    print('--------------------------------------------------------------------------------------')
    print('Zs       = '+str(self.Zs))
    print('mHats    = '+str(self.mHats))
    if self.includePhi1:
      print('includePhi1 = True')
      fluxTerm='_vd_'+radialCoord
    else:
      fluxTerm='_vm_'+radialCoord
    if self.withAdiabatic:
      print('withAdiabatic = True, adiabaticZ='+str(self.adiabaticZ)+', adiabaticMHat='+str(self.adiabaticMHat))
    if self.withNBIspec:
      print('withNBIspec =True, NBIspecZ='+str(self.NBIspecZ))

    print(radialCoord+' = ')
    print(getattr(self,radialCoord))
    print('particleFlux'+fluxTerm+' = ')
    print(getattr(self,'particleFlux'+fluxTerm))
    print('heatFlux'+fluxTerm+' = ')
    print(getattr(self,'heatFlux'+fluxTerm))
    print('momentumFlux'+fluxTerm+' = ')
    print(getattr(self,'momentumFlux'+fluxTerm))
    print('FSABFlow = ')
    print(self.FSABFlow)


  def choose_Erscan_run_with_best_Er(self,ErQuantity='dPhiHatdrN'):
    # returns the index of the run with the lowest radial current
    if not(self.sortafter=='dPhiHatdrN' or self.sortafter=='dPhiHatdrHat' or self.sortafter=='dPhiHatdpsiN'):
      sys.exit("Error using Ersearch: The sfincsScan data was not sorted after radial electric field. "+
               "Please use the option sortafter='dPhiHatdrN' in the initialisation.")
    Er=getattr(self,ErQuantity)
    if np.any(np.diff(Er))==0:
      sys.exit('Er values are not unique. This is probably not an Er scan!')
    if self.includePhi1:
      particleFlux=self.particleFlux_vd_rHat
    else:
      particleFlux=self.particleFlux_vm_rHat
    absJr=np.abs(np.sum(particleFlux*self.Zs,axis=1))
    # Unbeliavable that the current is zero to machine precision! This is instead because there was no result stored!
    goodindxs=np.where(absJr!=0.0)[0]
    return goodindxs[np.argmin(absJr[goodindxs])]

  def Ersearch(self,ErQuantity='dPhiHatdrN',verbose=0,launch='no',launchCommand='sbatch',interptype='quad',
               jobfilefrom='nearest'):
    #launch can be 'yes', 'no' or 'ask'
    if not(self.sortafter=='dPhiHatdrN' or self.sortafter=='dPhiHatdrHat' or self.sortafter=='dPhiHatdpsiN'):
      sys.exit("Error using Ersearch: The sfincsScan data was not sorted after radial electric field. "+
               "Please use the option sortafter='dPhiHatdrN' in the initialisation.")
    if self.Nruns<2:
      sys.exit('Error using Ersearch: Number of runs is less than 2!')

    Er=getattr(self,ErQuantity)
    if self.includePhi1:
      particleFlux=self.particleFlux_vd_rHat
    else:
      particleFlux=self.particleFlux_vm_rHat
    Jr=np.sum(particleFlux*self.Zs,axis=1)
    if verbose>0 or launch=='ask': 
      #print '---------'
      #print particleFlux
      #print np.abs(particleFlux*self.Zs)
      #print Jr
      #print np.max(np.abs(particleFlux*self.Zs),axis=1)
      rel_error=np.abs(Jr)/np.max(np.abs(particleFlux*self.Zs),axis=1)
      #rel_error_min=np.minimum(np.abs(Jr/np.max(np.abs(particleFlux*self.Zs),axis=0)))
      print(ErQuantity+',   Jr (a.u),   relative error')
      for ind in range(self.Nruns):
        if Jr[ind]!=0.0:
          print('{:10.8e} {: 8.6e} {:8.6e}'.format(Er[ind],Jr[ind],rel_error[ind]))
        else:
          print('{:10.8e}        nan       nan'.format(Er[ind]))

    if np.any(Jr==0.0) or np.any(np.isnan(Jr)):
      #print('removing nans')
      # Unbeliavable that the current is zero to machine precision! This is instead beacuse there was no result stored!
      goodindxs=np.where(np.logical_and(Jr!=0.0,np.logical_not(np.isnan(Jr))))[0]
      #print(goodindxs)
      Er = Er[goodindxs]
      Jr = Jr[goodindxs]
      rel_error = rel_error[goodindxs]
    
    #print 'Jr='+str(Jr)
    #print np.sign(Jr)
    #print np.diff(np.sign(Jr))==0.0
    
    if len(Jr)<2:
      print('Less than two (finished) runs. Skipping this radius.')
    else:
      extrapNeeded='No'
      if np.all(np.diff(np.sign(Jr))==0.0):
        extrapNeeded='Right'
        if abs(Jr[0])<abs(Jr[-1]):
          extrapNeeded='Left'

      if np.any(np.diff(Er)==0.0):
        sys.exit('Error using Ersearch: The radial electric fields in this scan '+self.mainDir+' are not unique! Perhaps not an Er scan.')
      if self.Nruns==2 or extrapNeeded=='Left':
        newEr = Er[0]+(Er[1]-Er[0])/(Jr[1]-Jr[0])*(0.0-Jr[0])
      elif extrapNeeded=='Right':
        newEr = Er[-1]+(Er[-2]-Er[-1])/(Jr[-2]-Jr[-1])*(0.0-Jr[-1])
      elif interptype=='lin' or len(Jr)==2:
        Lind=np.where(np.diff(np.sign(Jr))!=0)[0][0]
        newEr=Er[Lind]+(Er[Lind+1]-Er[Lind])/(Jr[Lind+1]-Jr[Lind])*(0.0-Jr[Lind])
      else:
        Lind=np.where(np.diff(np.sign(Jr))!=0)[0][0]
        Hind=Lind+1
        if self.Nruns==3 or Lind==0:
          #define c by Er=Er[base]+c[0]*(Jr-Jr[base])+c[1]*(Jr-Jr[base])^2
          base=0
          DJ1=Jr[base+1]-Jr[base]
          DJ2=Jr[base+2]-Jr[base]
          A=np.array([[DJ2**2,-DJ1**2],[-DJ2,DJ1]])
          c=1/(DJ1*DJ2*(DJ2-DJ1))*np.matmul(A,Er[base+1:base+3]-Er[base])
          newEr = Er[base]+c[0]*(0.0-Jr[base])+c[1]*(0.0-Jr[base])**2
        elif Lind==len(Jr)-2:
          base=len(Jr)-3
          DJ1=Jr[base+1]-Jr[base]
          DJ2=Jr[base+2]-Jr[base]
          A=np.array([[DJ2**2,-DJ1**2],[-DJ2,DJ1]])
          c=1/(DJ1*DJ2*(DJ2-DJ1))*np.matmul(A,Er[base+1:base+3]-Er[base])
          newEr = Er[base]+c[0]*(0.0-Jr[base])+c[1]*(0.0-Jr[base])**2
        else:
          #Make two quadratic approximations
          bas1=Lind-1
          bas2=Lind

          DJ1=Jr[bas1+1]-Jr[bas1]
          DJ2=Jr[bas1+2]-Jr[bas1]
          A=np.array([[DJ2**2,-DJ1**2],[-DJ2,DJ1]])
          c=1/(DJ1*DJ2*(DJ2-DJ1))*np.matmul(A,Er[bas1+1:bas1+3]-Er[bas1])
          Er1=Er[bas1]+c[0]*(0.0-Jr[bas1])+c[1]*(0.0-Jr[bas1])**2

          DJ1=Jr[bas2+1]-Jr[bas2]
          DJ2=Jr[bas2+2]-Jr[bas2]
          A=np.array([[DJ2**2,-DJ1**2],[-DJ2,DJ1]])
          c=1/(DJ1*DJ2*(DJ2-DJ1))*np.matmul(A,Er[bas2+1:bas2+3]-Er[bas2])
          Er2=Er[bas2]+c[0]*(0.0-Jr[bas2])+c[1]*(0.0-Jr[bas2])**2
          
          x=(0.0-Jr[Lind])/(Jr[Hind]-Jr[Lind])
          newEr = Er1*(1-x)+Er2*x
        #end if self.Nruns==3 or Lind==0
        if (newEr-Er[Lind])*(newEr-Er[Lind+1])>0: #newEr is not in the interval => use linear approximation
          newEr=Er[Lind]+(Er[Lind+1]-Er[Lind])/(Jr[Lind+1]-Jr[Lind])*(0.0-Jr[Lind])
      #end if self.Nruns==2
      closestind=np.argmin(np.abs(newEr-Er))

      if launch=='ask':
        answer=inp('Launch the calculation at the new value: '+ErQuantity + ' = {:8.6f} ? (ret=yes,n=no):'.format(newEr))
        if len(answer)==0:
          launch='yes'

      if launch=='yes':
        newDataDir=self.mainDir + '/' + ErQuantity + '{:8.6f}'.format(newEr)
        print(newDataDir)
        launchindeed=False
        if not(os.path.isdir(newDataDir)):
          os.mkdir(newDataDir) 
          launchindeed=True
        else:
          print('The directory already exists. It might contain a failed calculation.')
          answer=inp('Launch the calculation again? (ret=yes,n=no):')        
          if len(answer)==0:
            launchindeed=True
        if launchindeed:
          #copy jobfile
          print('jobfilefrom='+jobfilefrom)
          print('maindir='+self.mainDir)
          print('closest datadir='+self.DataDirs[closestind])
          if jobfilefrom=='nearest':
             jobfilesource=self.mainDir + '/' + self.DataDirs[closestind]+'/job.sfincsScan'
          else:
             jobfilesource=jobfilefrom
          with open(jobfilesource) as f:
            oldjob_file=f.readlines()
          newjob_fid=open(newDataDir+'/job.sfincsScan','w')
          for line in oldjob_file:
            newjob_fid.write(line)
          newjob_fid.close()

          #copy input.namelist
          with open(self.mainDir + '/' + self.DataDirs[closestind]+'/input.namelist') as f:
            oldnamelist_file=f.readlines()
          newnamelist_fid=open(newDataDir+'/input.namelist','w')
          for line in oldnamelist_file:
            startind=line.find(ErQuantity)
            if startind<0:
              newnamelist_fid.write(line)
            else:
              newnamelist_fid.write(line[:startind]+ErQuantity+' = '+'{:8.6f}\n'.format(newEr))
          newnamelist_fid.close()

          env = dict(os.environ)
          stat=subprocess.call([launchCommand,'job.sfincsScan'],cwd=newDataDir,env=env)
          if stat > 0:
            print('Error submitting the file '+newDataDir+'/job.sfincsScan with '+launchCommand+' !')
            sys.exit(stat)

    return newEr

  def plot(self,xvarName,yvarNames):
    print(xvarName)
    xvar=getattr(self,xvarName)
    print(xvar)
    if not(isinstance(yvarNames,(list,tuple))):
      yvarNames=[yvarNames]
    #print(len(yvarNames))
    fig,ax=plt.subplots(len(yvarNames),1,sharex=True)
    axv=ax
    if len(yvarNames)==1:
     print('hej')
     axv=[ax]

    for iy in range(len(yvarNames)):
      print(iy)
      print(yvarNames[iy])
      yvar=getattr(self,yvarNames[iy])
      print(yvar)
      axv[iy].plot(xvar,yvar)
      axv[iy].set_ylabel(yvarNames[iy])
 
    axv[-1].set_xlabel(xvarName)
    return fig,ax    

  def plotConvScan(self,convParams=None,quantities=None,yscale='relative',radialCoord='rHat'):
    #find out which was the baseCase
    if not 'baseCase' in self.DataDirs:
      sys.exit('This is not a convergence scan!')
    baseCaseind=self.DataDirs.index('baseCase')
    Nspec=len(self.nHats[0])

    #check which variables were scanned
    if convParams is None:
      convParams=[]
      if np.any(np.diff(self.Ntheta)!=0.0):
        convParams.append('Ntheta')
      if np.any(np.diff(self.Nzeta)!=0.0):
        convParams.append('Nzeta')
      if np.any(np.diff(self.Nxi)!=0.0):
        convParams.append('Nxi')
      if np.any(np.diff(self.Nx)!=0.0):
        convParams.append('Nx')
      if np.any(np.diff(self.NL)!=0.0):
        convParams.append('NL')
      if np.any(np.diff(self.solverTolerance)!=0.0):
        convParams.append('solverTolerance')
    NParams=len(convParams)
    paraminds=[None]*NParams
    baseParams=[None]*NParams
    for ind in range(NParams):
      baseParams[ind]=getattr(self,convParams[ind])[baseCaseind]
      paraminds[ind]=np.where(getattr(self,convParams[ind])!=baseParams[ind])
    
    if quantities is None:
      quantities=['particleFlux','heatFlux','FSABFlow']
    Nquant=len(quantities)
    if self.includePhi1:
      vmvdstr='_vd_'
    else:
      vmvdstr='_vm_'
    for ind in range(Nquant):
      if quantities[ind] in ['particleFlux','heatFlux','momentumFlux']:
        quantities[ind]=quantities[ind]+vmvdstr+radialCoord
    
    figs=[None]*Nspec
    axes=[None]*Nspec
    for iSpec in range(Nspec):
      figs[iSpec],axes[iSpec]=plt.subplots(Nquant,NParams,figsize=(1+Nquant*2,1+NParams*2))
      for qi in range(Nquant):
        thisQuant=getattr(self,quantities[qi])[:,iSpec]
        if yscale=='rel' or yscale=='relative':
          thisQuant=thisQuant/thisQuant[baseCaseind]
        for pi in range(NParams):
          thisParam=getattr(self,convParams[pi])
          vari=paraminds[pi]
          ax=axes[iSpec][qi,pi]
          ax.plot(thisParam[vari],thisQuant[vari],'b+')
          ax.plot(thisParam[baseCaseind],thisQuant[baseCaseind],'r+')
          ax.set_xlabel(convParams[pi])
          if pi==0:
            ax.set_ylabel(quantities[qi])
    return figs,axes
    

  def save(self,filename,varNames):
    if not(filename.endswith('.dat')):
      sys.exit('Only saving to text files (.dat) possible!')
    if not(isinstance(varNames,(list,tuple))):
      varNames=[varNames]
    varNameList=[]
    varList=[]
    attrlist=dir(self)
    for vn in varNames:
      if vn not in attrlist:
        print('The variable '+vn+' does not exist. Not saving it!')
      else: 
        var=getattr(self,vn)
        if not(isinstance(var,(list,np.ndarray))):
          varList.append(var*np.ones((self.Nruns)))
          varNameList.append(vn)
        elif isinstance(var,np.ndarray) and len(var.shape)==2:
          for iSpec in range(self.Nspecies):
            varList.append(var[:,iSpec])
            varNameList.append(vn+'_'+str(iSpec+1))  
        else: #1D array, length=Nruns
          varList.append(np.array(var)) 
          varNameList.append(vn)

    f = open(filename, 'w')
    #write header line
    f.write(' '.join(varNameList))
    f.write('\n')  
    for i in range(self.Nruns):
        for j in range(len(varList)):
            f.write('{:11.7e} '.format(varList[j][i]))
        f.write('\n')
    f.close()
      

#################################################################################################################
#################################################################################################################
class sfincsRadialAndErScan:
#################################################################################################################
#################################################################################################################
  def __init__(self,headDirectory,verbose=0):
    #print('In sfincsRadialAndErScan')
    #print('verbose='+str(verbose))
    if headDirectory is None:
      sys.exit('Missing input headDirectory!')
    if headDirectory[-1]=='/':
       headDirectory=headDirectory[:-1]
    self.headDir=headDirectory
    if verbose>0:
      print("Extracting data from radial scan of Er scans in " + headDirectory)
    # Get a list of the subdirectories:
    dirList=os.listdir(headDirectory)
    CandidateDirs=[]
    for ind in range(len(dirList)):
      if os.path.isdir(headDirectory+'/'+dirList[ind]):
        CandidateDirs.append(dirList[ind])

    CandidateDirs = sorted(CandidateDirs)
    if len(CandidateDirs) < 1:
      print('Error! Could not find any directories in ' + headDirectory)
      sys.exit(1)

    unsrtErscans=[]
    for ind in range(len(CandidateDirs)):
      if verbose>1:
        print('Loading: '+headDirectory+'/'+CandidateDirs[ind])
        tmp=sfincsScan(headDirectory+'/'+CandidateDirs[ind],sortafter='dPhiHatdrN',verbose=verbose)
        unsrtErscans.append(tmp)
      else:
        try:
          tmp=sfincsScan(headDirectory+'/'+CandidateDirs[ind],sortafter='dPhiHatdrN',verbose=verbose)
          unsrtErscans.append(tmp)
        except:
          print('Could not load '+headDirectory+'/'+CandidateDirs[ind])

    self.Nradii=len(unsrtErscans)
    unsrt_rN=np.zeros((self.Nradii))
    for ind in range(self.Nradii):
      if np.any(unsrtErscans[ind].rN!=unsrtErscans[ind].rN[0]):
        sys.exit('In the directory '+unsrtErscans[ind].mainDir+', not all rN are equal. Is it probably not an Er scan!')
      unsrt_rN[ind]=unsrtErscans[ind].rN[0]

    sorti=np.argsort(unsrt_rN)
    self.rN=unsrt_rN[sorti]
    self.Erscans=[None]*self.Nradii
    for ind in range(self.Nradii):
       self.Erscans[ind]=unsrtErscans[sorti[ind]]
    
  def disp(self):
    print('-----------------------------------------------------------------------------------')
    print(' Radial and Er scan over '+str(self.Nradii)+' radii')
    print('-----------------------------------------------------------------------------------')
    for ind in range(self.Nradii):
      print ('radius ind {:2d}'.format(ind)+
             ', rN={:6.4f}'.format(self.rN[ind])+
             ', mainDirectory: '+self.Erscans[ind].mainDir)

  def Ersearch(self,ErQuantity='dPhiHatdrN',verbose=0,launch='no',interptype='quad',jobfilefrom='nearest'):
    #launch can be 'yes', 'no' or 'ask'
    newErQ=np.nan*np.zeros((self.Nradii))
    for ind in range(self.Nradii):
      print('------------------')
      print('radius index '+str(ind))
      print('------------------')
      if jobfilefrom=='parentdir':
         jobfilefrom=self.headDir+'/job.sfincsScan'
      if launch=='no':
        newErQ[ind]=self.Erscans[ind].Ersearch(ErQuantity=ErQuantity,verbose=verbose,
                                               interptype=interptype,jobfilefrom=jobfilefrom)
      else:
        newErQ[ind]=self.Erscans[ind].Ersearch(ErQuantity=ErQuantity,verbose=1,launch=launch,
                                               interptype=interptype,jobfilefrom=jobfilefrom)
        
        
    return newErQ
      

    

    
