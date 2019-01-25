from __future__ import division
import numpy as np
import f90nml
import h5py

import subprocess
from warnings import warn

from shutil import copytree
import os

from geomlib import bcgeom
from fluxcoorddiscr import fluxcoorddiscr

from calculate_classical_transport import calculate_classical_transport



class Normalization(object):
    def __init__(self,BBar,RBar,TBar,mBar,nBar,eBar,ePhiBar,units="SI"):
        if units == "SI":
            self.units=units
        else:
            print "units '" + units + "' are not supported. Supported units are: 'SI'"
        
        self.BBar=BBar
        self.RBar=RBar
        self.TBar=TBar
        self.mBar=mBar
        self.nBar=nBar
        self.eBar=eBar
        self.ePhiBar=ePhiBar
        self.vBar=np.sqrt(2*TBar/mBar)
        self.alpha = ePhiBar/TBar
        self.Delta = mBar*self.vBar/(eBar*BBar*RBar)

    def __str__(self):
        return str(self.RBar) + "," +str(self.BBar) + "," +str(self.TBar) + "," +str(self.mBar) + "," +str(self.nBar) + "," +str(self.eBar) + "," +str(self.ePhiBar)

    

def create_normalization(filename="norms.namelist"):
    """ creates a Normalization object from the filename of a namelist of normalizing quantities"""
    filename="." + "/" + filename
    file = f90nml.read(filename)
    units=file["unitSpecification"]["units"]
    norm_group_name="normalizationParameters"
    BBar=file[norm_group_name]["BBar"]
    RBar=file[norm_group_name]["RBar"]
    TBar=file[norm_group_name]["TBar"]
    mBar=file[norm_group_name]["mBar"]
    nBar=file[norm_group_name]["nBar"]
    eBar=file[norm_group_name]["eBar"]
    ePhiBar=file[norm_group_name]["ePhiBar"]
    return Normalization(BBar,RBar,TBar,mBar,nBar,eBar,ePhiBar,units)

##############################
                                                                                                        


class Species(object):
    def __init__(self,Zs,mHats,names):
        self.names = names
        self.Zs = Zs
        self.mHats = mHats
    def __str__(self):
        return str(self.names) + "," + str(self.Zs) + "," + str(self.mHats)

def create_species(normalization,filename="species",database_filename = os.path.abspath(__file__).rsplit("/",1)[0] + "/species_database.namelist"):
    """ creates a Species object from a .csv file containing species, and a Normalization object"""

    filename="." + "/" + filename
    names=[x.strip() for x in open(filename,'r').read().split('\n')[0].split(',')]
    database = f90nml.read(database_filename)
    Zs = np.array([database["speciesCharge"][x] for x in names])/normalization.eBar
    mHats = np.array([database["speciesMass"][x] for x in names])/normalization.mBar
    return Species(Zs,mHats,names)

##############################
             

class Sfincs_aux_inputs(object):
    
    def __init__(self,norm_filename='norms.namelist',species_filename='species'):
        self.norm_filename = norm_filename
        self.species_filename = species_filename
        
        
        self.normalization = create_normalization(norm_filename)
        self.species = create_species(species_filename)
        
        self.eBar = self.normalization.eBar
        self.RBar = self.normalization.RBar
        self.nBar = self.normalization.nBar
        self.Tbar = self.normalization.TBar
        self.mBar = self.normalization.mBar
        self.vBar = self.normalization.vBar
        self.BBar = self.normalization.BBar
        self.ePhiBar = self.normalization.ePhiBar
        self.Zs = self.species.Zs
        self.mHats = self.species.mHats
        self.species_name = self.species.names
        self.Delta = self.normalization.Delta
        self.alpha = self.normalization.alpha
        
    def __str__(self):
        return str(self.normalization) + "\n" + str(self.species)
        
##############################

class Sfincs_input(object):

    defaults = {'general' : {'rhsmode' : 1,
                             'outputfilename' : 'sfincsOutput.h5',
                             'savematlaboutput' : False,
                             'solvesystem' : True,
                             },
                'geometryparameters' : {'geometryscheme' : 1,
                                        'inputradialcoordinate' : 3,
                                        'inputradialcoordinateforgradients' : 4,
                                        'equilibriumfile' : "",
                                        'vmecradialoption' : 1,
                                        'vmec_nyquist_option': 1,
                                        'min_bmn_to_load' : 0.0,
                                    },
                'physicsparameters' : {'delta' : 4.5694e-3,
                                       'alpha': 1.0,
                                       'nu_n' : 8.330e-3,
                                       'er' : 0.0,
                                       'dphihatdpsihat' : 0.0,
                                       'dphihatdpsin' : 0.0,
                                       'dphihatdrhat' : 0.0,
                                       'dphihatdrn' : 0.0,
                                   },
                'speciesparameters' : {'zs' : [1.0],
                                       'mhats' : [1.0],
                                       'nhats' : [1.0],
                                       'thats' : [1.0],
                                       'dnHatdpsiHats' : [0.0],
                                       'dnhatdpsins' : [0.0],
                                       'dnhatdrhats' : [0.0],
                                       'dnhatdrns' : [0.0],
                                       'dthatdpsihats' : [0.0],
                                       'dthatdpsins' : [0.0],
                                       'dthatdrhats' : [0.0],
                                       'dthatdrns' : [0.0],
                                   },
                'resolutionparameters' : {'forceoddnthetaandnzeta' : True,
                                          'xmax': 5.0,
                                      },
            }


    def changevar(self,group,var,value):
        # Warning: this command will fail silently if the pattern is not found. Sorry about that.
        # Warning: case insensitive
        if type(value) == str:
            #strings must be enclosed in "" in namelists
            #may be wise to see if the string contains citation marks...
            if (value.find("'") != -1) or (value.find('"') != -1):
                print "Warning! String to changevar contains a ' or \" character." 
            value = '"' + value + '"'
        if (type(value) == list) or (type(value) == np.ndarray):
            # arrays are space seperated
            delimiter=' '
            value_temp = '' 
            for val in value:
                value_temp =  value_temp + str(val) + delimiter
            value = value_temp.rsplit(delimiter,1)[0]
            
        subprocess.call("sed -i -e '/\&"+group+"/I,/\&/{ s/^  "+var+" =.*/  "+var+" = "+str(value)+"/I } ' "+self.input_name, shell=True)
    
    def get_value_from_input_or_defaults(self,groupname,varname):
        varname = varname.lower()
        groupname = groupname.lower()
        inputs = f90nml.read(self.input_name)
        if not varname in inputs[groupname].keys():
            return Sfincs_input.defaults[groupname][varname]
        else:
            return inputs[groupname][varname]
    
    def __init__(self,input_name):
        self.input_name = input_name

    @property
    def output_filename(self):
        return self.get_value_from_input_or_defaults("general","outputfilename")

    @property
    def Delta(self):
        return self.get_value_from_input_or_defaults("physicsparameters","delta")

    @property
    def alpha(self):
        return self.get_value_from_input_or_defaults("physicsparameters","alpha")
    
    @property
    def nu_n(self):
        return self.get_value_from_input_or_defaults("physicsparameters","nu_n")

    @property
    def Zs(self):
        return self.get_value_from_input_or_defaults("speciesparameters","zs")

    @property
    def mHats(self):
        return self.get_value_from_input_or_defaults("speciesparameters","mhats")

    @property
    def nHats(self):
        return self.get_value_from_input_or_defaults("speciesparameters","nhats")

    @property
    def THats(self):
        return self.get_value_from_input_or_defaults("speciesparameters","thats")

    @property
    def Ntheta(self):
        return self.get_value_from_input_or_defaults("resolutionparameters","Ntheta")

    @property
    def Nzeta(self):
        return self.get_value_from_input_or_defaults("resolutionparameters","Nzeta")

    @property
    def Nxi(self):
        return self.get_value_from_input_or_defaults("resolutionparameters","Nxi")

    @property
    def Nx(self):
        return self.get_value_from_input_or_defaults("resolutionparameters","Nx")

    @property
    def NL(self): 
        return self.get_value_from_input_or_defaults("resolutionparameters","NL")

    @property
    def forceOddNthetaAndNzeta(self):
        return self.get_value_from_input_or_defaults("resolutionparameters","forceOddNthetaAndNzeta")
        

    @property
    def geometryScheme(self):
        return self.get_value_from_input_or_defaults("geometryparameters","geometryscheme")

    @property
    def inputRadialCoordinate(self):
        return self.get_value_from_input_or_defaults("geometryparameters","inputradialcoordinate")

    @property
    def inputRadialCoordinateForGradients(self):
        return self.get_value_from_input_or_defaults("geometryparameters","inputradialcoordinateforgradients")

    @property
    def equilibrium_name(self):
        return self.get_value_from_input_or_defaults("geometryparameters","equilibriumFile")

    @property
    def coordinate_wish(self):
        if self.inputRadialCoordinate == 0:
            return self.get_value_from_input_or_defaults("geometryparameters","psiHat_wish")
        elif self.inputRadialCoordinate == 1:
            return self.get_value_from_input_or_defaults("geometryparameters","psiN_wish")
        elif (self.inputRadialCoordinate == 2) or (self.inputRadialCoordinate == 4):
            return self.get_value_from_input_or_defaults("geometryparameters","rHat_wish")
        elif self.inputRadialCoordinate == 3:
            return self.get_value_from_input_or_defaults("geometryparameters","rN_wish")
        else:
            raise ValueError("inputRadialCoordinate should be 0,1,2,3,4; it is" + str(inputRadialCoordinate))

    @property
    def dTHatdss(self):
        if self.inputRadialCoordinateForGradients == 0:
            ret= self.get_value_from_input_or_defaults("speciesparameters","dTHatdpsiHats")
        elif self.inputRadialCoordinateForGradients == 1:
            ret= self.get_value_from_input_or_defaults("speciesparameters","dTHatdpsiNs")
        elif (self.inputRadialCoordinateForGradients == 2) or (self.inputRadialCoordinateForGradients == 4):
            ret= self.get_value_from_input_or_defaults("speciesparameters","dTHatdrHats")
        elif self.inputRadialCoordinateForGradients == 3:
            ret= self.get_value_from_input_or_defaults("speciesparameters","dTHatdrNs")
        else:
            raise ValueError("inputRadialCoordinateForGradients should be 0,1,2,3,4; it is" + str(inputRadialCoordinate))
        return np.array(ret)

    @property
    def dnHatdss(self):
        if self.inputRadialCoordinateForGradients == 0:
            ret= self.get_value_from_input_or_defaults("speciesparameters","dnHatdpsiHats")
        elif self.inputRadialCoordinateForGradients == 1:
            ret= self.get_value_from_input_or_defaults("speciesparameters","dnHatdpsiNs")
        elif (self.inputRadialCoordinateForGradients == 2) or (self.inputRadialCoordinateForGradients == 4):
            ret= self.get_value_from_input_or_defaults("speciesparameters","dnHatdrHats")
        elif self.inputRadialCoordinateForGradients == 3:
            ret= self.get_value_from_input_or_defaults("speciesparameters","dnHatdrNs")
        else:
            raise ValueError("inputRadialCoordinateForGradients should be 0,1,2,3,4; it is" + str(inputRadialCoordinate))
        return np.array(ret)

    @property
    def dPhiHatds(self):
        if self.inputRadialCoordinateForGradients == 0:
            ret= self.get_value_from_input_or_defaults("physicsparameters","dPhiHatdpsiHat")
        elif self.inputRadialCoordinateForGradients == 1:
            ret= self.get_value_from_input_or_defaults("physicsparameters","dPhiHatdpsiN")
        elif self.inputRadialCoordinateForGradients == 2:
            ret= self.get_value_from_input_or_defaults("physicsparameters","dPhiHatdrHat")
        elif self.inputRadialCoordinateForGradients == 3:
            ret= self.get_value_from_input_or_defaults("physicsparameters","dPhiHatdrN")
        elif self.inputRadialCoordinateForGradients == 4:
            ret= -self.get_value_from_input_or_defaults("physicsparameters","Er")
            
        else:
            raise ValueError("inputRadialCoordinateForGradients should be 0,1,2,3,4; it is" + str(inputRadialCoordinate))
        return np.array(ret)

    @property
    def A1(self):
        ret = self.dnHatdss/self.nHats + self.dTHatdss/self.THats + self.Zs * self.dPhiHatds/self.THats
        return np.array(ret)

    @property
    def A2(self):
        ret = self.dTHatdss/self.THats
        return np.array(ret)
    
    def __str__(self):
        return str(self.input_name)

##############################

class Sfincs_simulation(object):

    def absolute_path(self,path):
        if path[0]=='/':
            return path
        else:
            return self.dirname + "/" + path

    def clean(self):
        pass

    def copy(self,new_dirname):
        copytree(self.dirname,new_dirname)
        return Sfincs_simulation(new_dirname,input_name=self.input_name,norm_name=self.norm_name,species_name=self.species_name)
        
    def __init__(self,dirname,input_name="input.namelist",norm_name="norm.namelist",species_name="species",override_geometry_name=None):
        #description of simulation, for usage as legend in plots, etc.
        self.description=""

        self.dirname = dirname
        self.input_name = input_name
        self.norm_name = norm_name
        self.species_name = species_name
        
        
        self.input = Sfincs_input(self.absolute_path(self.input_name))
        self.normalization = create_normalization(self.absolute_path(self.norm_name))
        self.species = create_species(self.normalization,self.absolute_path(self.species_name))
        try:
            self.outputs=h5py.File(self.absolute_path(self.input.output_filename),'r')
        except IOError:
            self.hasOutput=False
        else:
            self.hasOutput=True

        if override_geometry_name is not None:
            self.equilibrium_name = override_geometry_name
        else:
            self.equilibrium_name = self.absolute_path(self.input.equilibrium_name)
            
        if self.input.geometryScheme == 11:
            # known problem: if self.input.equilibrium_name is changed
            # these quanitities will no longer correspond to those indicated by the input.namelist
            # in some sense, these belong to the input_file and not the simulation, but the input file doesn't right now know about it's directory and hence can't translate to absolute path.
            self.symmetry='StelSym'
            self.min_Bmn=0
            self.max_m=float("inf")
            self.maxabs_n=float("inf")
            self.signcorr=1
            verbose = 0

            self.zeroout_Deltaiota = 0.000

            try:
                self.geom = bcgeom(self.equilibrium_name,self.min_Bmn,self.max_m,self.maxabs_n,self.symmetry,self.signcorr,verbose)
            except IOError as e:
                raise IOError("Error loading geometry: " + str(e))
            
            self.Nperiods = self.geom.Nperiods
            self.psiAHat = self.geom.psi_a/(self.normalization.BBar*self.normalization.RBar**2)
            self.aHat = self.geom.minorradiusVMEC/self.normalization.RBar
            
            self.rNs = self.geom.rnorm
            self.rind=np.argmin(np.fabs(self.rNs - self.rN_wish))
            self.iota = self.geom.iota[self.rind]
            
            self.Booz = fluxcoorddiscr(self.geom,self.rind,self.input.Ntheta,self.input.Nzeta,u_zeroout_Deltaiota=self.zeroout_Deltaiota,name='Boozer')
            self.GHat=self.Booz.G/(self.normalization.BBar*self.normalization.RBar)
            self.IHat=self.Booz.I/(self.normalization.BBar*self.normalization.RBar)
            self.B00Hat = self.Booz.B00/(self.normalization.BBar)

    @property
    def gpsipsiHat(self):
        if self.input.geometryScheme == 11:
            return np.transpose(self.Booz.gpsipsi)/((self.normalization.BBar*self.normalization.RBar)**2)
        else:
            try:
                return self.outputs["gpsiHatpsiHat"][()]
            except KeyError:
                raise NotImplementedError("gpsipsiHat cannot be calculated for this simulation since it does not use a .bc file and does not contain gpsipsi in the output.")
            
    @property
    def BHat(self):
        if self.input.geometryScheme == 11:
            return np.transpose(self.Booz.B)/self.normalization.BBar
        else:
            try:
                return self.outputs["BHat"][()]
            except KeyError:
                raise NotImplementedError("BHat cannot be calculated for this simulation since it does not use a .bc file and does not contain BHat in the output.")

    @property
    def u(self):
        if self.input.geometryScheme == 11:
            return np.transpose(self.Booz.u_psi)
        else:
            raise NotImplementedError("u cannot be calculated for this simulation since it does not use a .bc file.")
        

    def FSA(self,X):
        if self.input.geometryScheme == 11:
            BHat = self.BHat
            # order of axis: zeta, theta, species, iteration
            # thus the below sum over zeta and theta
            return np.sum(X/BHat**2,axis=(0,1))/np.sum(1/BHat**2)
        else:
            raise NotImplementedError("Flux-surface average cannot be calculated for this simulation since it does not use a .bc file.")
            
    @property
    def mixed_col_NC_C_ratio(self):
        BHat = self.BHat
        u = self.u
        gpsipsi = self.gpsipsiHat
        FSABHat2 = self.FSA(BHat**2)
        ret = self.FSA(u**2 * BHat**2)*FSABHat2 - self.FSA(u*BHat**2)**2
        ret = ret/(self.FSA(gpsipsi/BHat**2) * self.FSA(BHat**2))
        return ret
        
    @property
    def rN_wish(self):
        if self.input.inputRadialCoordinate == 0:
            return np.sqrt(self.input.coordinate_wish/self.psiAHat)
        elif self.input.inputRadialCoordinate == 1:
            return np.sqrt(self.input.coordinate_wish)
        elif (self.input.inputRadialCoordinate == 2) or (self.input.inputRadialCoordinate == 4):
            return self.input.coordinate_wish/self.aHat 
        elif self.input.inputRadialCoordinate == 3:
            return self.input.coordinate_wish
        else:
            raise ValueError("inputRadialCoordinate should be 0,1,2,3,4; it is" + str(self.input.inputRadialCoordinate))

    @property
    def psiN(self):
        return self.geom.s[self.rind]

    @property
    def rN(self):
        return np.sqrt(self.psiN)

    @property
    def psiHat(self):
        return self.psiN * self.psiAHat

    @property
    def rHat(self):
        return self.rN * self.aHat

    

    @property
    def dnHatdpsiHats(self):
        if self.input.inputRadialCoordinateForGradients == 0:
            conversion_factor = 1.0
        elif self.input.inputRadialCoordinate == 1:
            conversion_factor = 1/self.psiAHat
        elif (self.input.inputRadialCoordinate == 2) or (self.input.inputRadialCoordinate == 4):
            conversion_factor = self.aHat/(2*self.psiAHat*np.sqrt(self.psiN))
        elif self.input.inputRadialCoordinate == 3:
            conversion_factor = 1/(2*self.psiAHat*np.sqrt(self.psiN))
        else:
            raise ValueError("inputRadialCoordinate should be 0,1,2,3,4; it is" + str(self.input.inputRadialCoordinate))
        return conversion_factor * self.input.dnHatdss

    @property
    def dTHatdpsiHats(self):
        if self.input.inputRadialCoordinateForGradients == 0:
            conversion_factor = 1.0
        elif self.input.inputRadialCoordinate == 1:
            conversion_factor = 1/self.psiAHat
        elif (self.input.inputRadialCoordinate == 2) or (self.input.inputRadialCoordinate == 4):
            conversion_factor = self.aHat/(2*self.psiAHat*np.sqrt(self.psiN))
        elif self.input.inputRadialCoordinate == 3:
            conversion_factor = 1/(2*self.psiAHat*np.sqrt(self.psiN))
        else:
            raise ValueError("inputRadialCoordinate should be 0,1,2,3,4; it is" + str(self.input.inputRadialCoordinate))
        return self.input.dTHatdss * conversion_factor
        
            
    @property
    def integerToRepresentTrue(self):
        return self.outputs["integerToRepresentTrue"]

    @property
    def includePhi1(self):
        return (self.outputs["includePhi1"][()] == self.integerToRepresentTrue)
    
    @property
    def GammaHat(self):
        if self.includePhi1:
            return self.outputs["particleFlux_vd_psiHat"][:,-1]
        else:
            return self.outputs["particleFlux_vm_psiHat"][:,-1]

    @property
    def GammaHat_C(self):
        if self.includePhi1:
            try:
                return self.outputs["classicalParticleFlux_psiHat"][:,-1]
            except KeyError:
                return calculate_classical_transport(self.input.Zs,self.input.mHats,self.input.nHats,self.input.THats,self.dnHatdpsiHats,self.dTHatdpsiHats,self.input.Delta,self.input.alpha,self.input.nu_n,self.gpsipsiHat,self.BHat,self.Phi1Hat)[0]
        else:
            try:
                return self.outputs["classicalParticleFluxNoPhi1_psiHat"][()]
            except KeyError:
                return calculate_classical_transport(self.input.Zs,self.input.mHats,self.input.nHats,self.input.THats,self.dnHatdpsiHats,self.dTHatdpsiHats,self.input.Delta,self.input.alpha,self.input.nu_n,self.gpsipsiHat,self.BHat)[0]

    @property
    def QHat(self):
        if self.includePhi1:
            return self.outputs["heatFlux_vd_psiHat"][:,-1]
        else:
            return self.outputs["heatFlux_vm_psiHat"][:,-1]

    @property
    def Phi1Hat(self):
        if self.includePhi1:
            return self.outputs["Phi1Hat"][:,:,-1]
        else:
            None
        
    @property
    def QHat_C(self):
        if self.includePhi1:
            try:
                return self.outputs["classicalHeatFlux_psiHat"][:,-1]
            except KeyError:
                return calculate_classical_transport(self.input.Zs,self.input.mHats,self.input.nHats,self.input.THats,self.dnHatdpsiHats,self.dTHatdpsiHats,self.input.Delta,self.input.alpha,self.input.nu_n,self.gpsipsiHat,self.BHat,self.Phi1Hat)[1]
        else:
            try:
                return self.outputs["classicalHeatFluxNoPhi1_psiHat"][()]
            except KeyError:
                return calculate_classical_transport(self.input.Zs,self.input.mHats,self.input.nHats,self.input.THats,self.dnHatdpsiHats,self.dTHatdpsiHats,self.input.Delta,self.input.alpha,self.input.nu_n,self.gpsipsiHat,self.BHat)[1]

    @property
    def FSABFlow(self):
        return self.outputs["FSABFlow"][:,-1]

    @property
    def collisionality(self):
        if self.input.geometryScheme == 11:
            B00Hat = self.B00Hat
            IHat = self.IHat
            GHat = self.GHat
            iota = self.iota
        RHat = np.fabs((GHat + iota*IHat)/B00Hat)
        nu_n = self.input.nu_n
        Zs = self.input.Zs
        nHats = self.input.nHats
        THats = self.input.THats

        Nspecies = len(Zs)
        nuHats = np.zeros(Nspecies)
        for a in range(Nspecies):
            for b in range(Nspecies):
                nuHats[a] =nuHats[a] +  (Zs[b]**2*nHats[b])
            nuHats[a] = nuHats[a] * (Zs[a]/THats[a])**2 
        nuHats = nuHats * nu_n *RHat
        return nuHats


    
    @property
    def collisionality_aa(self):
        if self.input.geometryScheme == 11:
            B00Hat = self.B00Hat
            IHat = self.IHat
            GHat = self.GHat
            iota = self.iota
        RHat = np.fabs((GHat + iota*IHat)/B00Hat)
        nu_n = self.input.nu_n
        Zs = self.input.Zs
        nHats = self.input.nHats
        THats = self.input.THats
        
        Nspecies = len(Zs)
        nuHats = np.zeros(Nspecies)
        for a in range(Nspecies):
            nuHats[a] =nuHats[a] +  Zs[a]**4*nHats[a]/THats[a]**2
        nuHats = nuHats * nu_n *RHat
        return nuHats

            
if __name__=="__main__":

    norm_filename = "norm.namelist"
    species_filename = "species"
    simul = Sfincs_simulation('.')

    print "GammaHat:"
    print simul.GammaHat
    print "GammaHat_C:"
    print simul.GammaHat_C
    print "QHat:"
    print simul.QHat
    print "QHat_C:"
    print simul.QHat_C
    print "FSABFlow:"
    print simul.FSABFlow

    print "Gamma_C/Gamma:"
    print simul.GammaHat_C/simul.GammaHat
    print "Q_C/Q:"
    print simul.QHat_C/simul.QHat

