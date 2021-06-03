from __future__ import division, print_function
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

elecharge = 1.60217662e-19

class Normalization(object):
    def __init__(self,BBar,RBar,TBar,mBar,nBar,eBar,ePhiBar,units="SI"):
        if units == "SI":
            self.units=units
        else:
            print("units '" + units + "' are not supported. Supported units are: 'SI'")
        
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
        self.PhiBar = self.ePhiBar/self.eBar
        
    def __str__(self):
        return str(self.RBar) + "," +str(self.BBar) + "," +str(self.TBar) + "," +str(self.mBar) + "," +str(self.nBar) + "," +str(self.eBar) + "," +str(self.ePhiBar)

    

def create_normalization(filename="norms.namelist"):
    """ creates a Normalization object from the filename of a namelist of normalizing quantities"""
    filename=filename
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

    filename=filename
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
        if type(value) == bool:
            if value == True:
                value = ".true."
            else:
                value = ".false."
        elif type(value) == str:
            #strings must be enclosed in "" in namelists
            #may be wise to see if the string contains citation marks...
            if (value.find("'") != -1) or (value.find('"') != -1):
                print("Warning! String to changevar contains a ' or \" character.")
            value = '"' + value + '"'
        elif (type(value) == list) or (type(value) == np.ndarray):
            # arrays are space seperated
            delimiter=' '
            value_temp = '' 
            for val in value:
                value_temp =  value_temp + str(val) + delimiter
            value = value_temp.rsplit(delimiter,1)[0]
        else:
            pass    
        subprocess.call("sed -i -e '/\&"+group+"/I,/\&/{ s/^  "+var+" =.*/  "+var+" = "+str(value)+"/I } ' "+self.input_name, shell=True)

    def changessvar(self,var,value):
        # Warning: this command will fail silently if the pattern is not found. Sorry about that.
        # Warning: case insensitive
        if type(value) == bool:
            if value == True:
                value = ".true."
            else:
                value = ".false."
        elif type(value) == str:
            #strings must be enclosed in "" in namelists
            #may be wise to see if the string contains citation marks...
            if (value.find("'") != -1) or (value.find('"') != -1):
                print("Warning! String to changevar contains a ' or \" character.")
            value = '"' + value + '"'
        elif (type(value) == list) or (type(value) == np.ndarray):
            # arrays are space seperated
            delimiter=' '
            value_temp = '' 
            for val in value:
                value_temp =  value_temp + str(val) + delimiter
            value = value_temp.rsplit(delimiter,1)[0]
        else:
            pass    
        subprocess.call("sed -i -e 's/^\!ss "+var+" =.*/\!ss "+var+" = "+str(value)+"/' "+self.input_name, shell=True)
    
    def get_value_from_input_or_defaults(self,groupname,varname):
        varname = varname.lower()
        groupname = groupname.lower()
        inputs = f90nml.read(self.input_name)
        if not varname in inputs[groupname].keys():
            return Sfincs_input.defaults[groupname][varname]
        else:
            return inputs[groupname][varname]#.split('!',1)[0].strip()
    
    def __init__(self,input_name):
        self.input_name = input_name

    def copy(self,new_dirname):
        # copy the directory and everything in it, not only
        # the input_namelist
        dirname = ("./" + self.input_name).rsplit("/",1)[0] +"/"
        copytree(dirname,new_dirname)
        return Sfincs_input(new_dirname)
        
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
        nu_n =  self.get_value_from_input_or_defaults("physicsparameters","nu_n")
        if nu_n<0:
            # TODO check that normalizations are correct
            lnLambda = (25.3e+0) - (1.15e+0)*np.log10(self.nHats[0]*(1e14)) + (2.30e+0)*np.log10(self.THats[0]*1000)
            eC = 1.6022e-19
            epsilon0 = 8.8542e-12
            mproton = 1.6726e-27
            nu_n = np.sqrt(mproton/(2*1000*eC)) * 4*np.sqrt(2*np.pi)*(1e20)*(eC**4)*lnLambda / (3*((4*np.pi*epsilon0)**2)*np.sqrt(mproton)*((1000*eC)**(1.5e+0)))
        return nu_n

    @property
    def Zs(self):
        return np.array(self.get_value_from_input_or_defaults("speciesparameters","zs"))

    @property
    def Nspecies(self):
        return len(self.Zs)

    @property
    def mHats(self):
        return np.array(self.get_value_from_input_or_defaults("speciesparameters","mhats"))

    @property
    def nHats(self):
        return np.array(self.get_value_from_input_or_defaults("speciesparameters","nhats"))

    @property
    def dnHatdrHats(self):
        return np.array(self.get_value_from_input_or_defaults("speciesparameters","dnhatdrhats"))

    @property
    def dTHatdrHats(self):
        return np.array(self.get_value_from_input_or_defaults("speciesparameters","dthatdrhats"))

    @property
    def dPhiHatdrHat(self):
        return np.array(self.get_value_from_input_or_defaults("physicsparameters","dphihatdrhat"))


    @property
    def Zeff(self):
        # exclude electrons from sum if electrons are included in the simulation
        noe = np.where(Zs!=-1)
        return np.sum(self.Zs[noe]**2*self.nHats[noe])/np.sum(self.Zs[noe]*self.nHats[noe])
    
    @property
    def THats(self):
        return np.array(self.get_value_from_input_or_defaults("speciesparameters","thats"))

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
        return Sfincs_simulation(new_dirname,input_name=self.input_name,norm_name=self.norm_name,species_name=self.species_name,override_geometry_name=self.equilibrium_name, load_geometry=self.load_geometry)
        
    def __init__(self,dirname,input_name="input.namelist",norm_name="norm.namelist",species_name="species",override_geometry_name=None, load_geometry=True, Booz=None, geom=None):
        #description of simulation, for usage as legend in plots, etc.
        self.description=""

        self.dirname = dirname
        self.input_name = input_name
        self.norm_name = norm_name
        self.species_name = species_name
        
        
        self.input = Sfincs_input(self.absolute_path(self.input_name))
        try:
            self.normalization = create_normalization(self.absolute_path(self.norm_name))
        except FileNotFoundError:
            self.normalization = Normalization(1,1,1.602176565e-16,1.672621898e-27,1e20,1.602176565e-19,1.602176565e-16,units="SI")
        try:
            self.species = create_species(self.normalization,self.absolute_path(self.species_name))
        except FileNotFoundError:
            pass
        try:
            self.outputs=h5py.File(self.absolute_path(self.input.output_filename),'r')
        except IOError:
            self.hasOutput=False
            print("no output: " + str(self.absolute_path(self.input_name)))
        else:
            self.hasOutput=True

        if override_geometry_name is not None:
            self.equilibrium_name = override_geometry_name
        else:
            self.equilibrium_name = self.absolute_path(self.input.equilibrium_name)
        self.geometry_loaded = False
        # halfloaded is set halfway during self.load_geometry() 
        # because some data has to be read to 
        # read in the correct radius in the Boozer file
        self.geometry_halfloaded = False
        if (self.input.geometryScheme == 11) and load_geometry:
            self.load_geometry(Booz,geom)
            
        
    def load_geometry(self,Booz,geom):
        # known problem: if self.input.equilibrium_name is changed
        # these quanitities will no longer correspond to those indicated by the input.namelist
        # in some sense, these belong to the input_file and not the simulation, but the input file doesn't right now know about it's directory and hence can't translate to absolute path.
        if geom is None:
            self.symmetry='StelSym'
            self.min_Bmn=0
            self.max_m=float("inf")
            self.maxabs_n=float("inf")
            self.signcorr=1
            verbose = 0

            self.zeroout_Deltaiota = 0.000

            try:
                print("Loading geometry " +self.equilibrium_name +" ...")
                self.geom = bcgeom(self.equilibrium_name,self.min_Bmn,self.max_m,self.maxabs_n,self.symmetry,self.signcorr,verbose)
            except IOError as e:
                raise IOError("Error loading geometry: " + str(e))
        else:
            self.geom = geom
        self._Nperiods = self.geom.Nperiods
        self._psiAHat = self.geom.psi_a/(self.normalization.BBar*self.normalization.RBar**2)
        self._aHat = self.geom.minorradiusW7AS/self.normalization.RBar

        self.rNs = self.geom.rnorm
        
        self.geometry_halfloaded = True
        self.rind=np.argmin(np.fabs(self.rNs - self.rN_wish))
        self._iota = self.geom.iota[self.rind]
        if Booz is None:
            self.Booz = fluxcoorddiscr(self.geom,self.rind,self.input.Ntheta,self.input.Nzeta,u_zeroout_Deltaiota=self.zeroout_Deltaiota,name='Boozer')
        else:
            self.Booz=Booz
        self._GHat=self.Booz.G/(self.normalization.BBar*self.normalization.RBar)
        self._IHat=self.Booz.I/(self.normalization.BBar*self.normalization.RBar)
        self._B00Hat = self.Booz.B00/(self.normalization.BBar)
        self.geometry_loaded = True
        
    
    @property
    def Nperiods(self):
        if self.geometry_loaded:
            return self._Nperiods
        else:
            try:
                return self.outputs["NPeriods"][()]
            except KeyError:
                raise NotImplementedError("Nperiods cannot be calculated for this simulation since it does not use a .bc file and does not contain NPeriods in the output.")

    @property
    def Nspecies(self):
        return self.input.Nspecies

    @property
    def Ntheta(self):
        return self.input.Ntheta

    @property
    def Nzeta(self):
        return self.input.Nzeta
    
    @property
    def Nx(self):
        return self.input.Nx

    @property
    def Nxi(self):
        return self.input.Nxi
    
    @property
    def NL(self):
        return self.input.NL
    
    @property
    def psiAHat(self):
        if self.geometry_halfloaded:
            return self._psiAHat
        else:
            try:
                return self.outputs["psiAHat"][()]
            except KeyError:
                raise NotImplementedError("psiAHat cannot be calculated for this simulation since it does not use a .bc file and does not contain psiAHat in the output.")
    
    @property
    def aHat(self):
        if self.geometry_halfloaded:
            return self._aHat
        else:
            try:
                return self.outputs["aHat"][()]
            except KeyError:
                raise NotImplementedError("aHat cannot be calculated for this simulation since it does not use a .bc file and does not contain aHat in the output.")

    @property
    def iota(self):
        if self.geometry_loaded:
            return self._iota
        else:
            try:
                return self.outputs["iota"][()]
            except KeyError:
                raise NotImplementedError("iota cannot be calculated for this simulation since it does not use a .bc file and does not contain iota in the output.")

    @property
    def B00Hat(self):
        if self.geometry_loaded:
            return self._B00Hat
        else:
            try:
                return self.outputs["B0OverBBar"][()]
            except KeyError:
                raise NotImplementedError("B00Hat cannot be calculated for this simulation since it does not use a .bc file and does not contain B00OverBBar in the output.")


    @property
    def GHat(self):
        if self.geometry_loaded:
            return self._GHat
        else:
            try:
                return self.outputs["GHat"][()]
            except KeyError:
                raise NotImplementedError("GHat cannot be calculated for this simulation since it does not use a .bc file and does not contain GHat in the output.")

    @property
    def IHat(self):
        if self.geometry_loaded:
            return self._IHat
        else:
            try:
                return self.outputs["IHat"][()]
            except KeyError:
                raise NotImplementedError("IHat cannot be calculated for this simulation since it does not use a .bc file and does not contain IHat in the output.")

    @property
    def gpsipsiHat(self):
        if self.geometry_loaded:
            return np.transpose(self.Booz.gpsipsi)/((self.normalization.BBar*self.normalization.RBar)**2)
        else:
            try:
                return self.outputs["gpsiHatpsiHat"][()]
            except KeyError:
                raise NotImplementedError("gpsipsiHat cannot be calculated for this simulation since it does not use a .bc file and does not contain gpsipsi in the output.")
            
    @property
    def BHat(self):
        if self.geometry_loaded:
            return np.transpose(self.Booz.B)/self.normalization.BBar
        else:
            try:
                return self.outputs["BHat"][()]
            except KeyError:
                raise NotImplementedError("BHat cannot be calculated for this simulation since it does not use a .bc file and does not contain BHat in the output.")

    @property
    def u(self):
        if self.geometry_loaded:
            return np.transpose(self.Booz.u_psi)
        else:
            raise NotImplementedError("u cannot be calculated for this simulation since it does not use a .bc file.")
        

    def FSA(self,X):
        BHat = self.BHat
        # order of axis: zeta, theta, species, iteration
        # thus the below sum over zeta and theta

        rank = len(X.shape)
        if rank == 2:
            return np.sum(X/BHat**2,axis=(0,1))/np.sum(1/BHat**2)
        elif rank == 3:
            return np.sum(X/BHat[:,:,np.newaxis]**2,axis=(0,1))/np.sum(1/BHat**2)
        elif rank == 4:
            X = X[:,:,:,-1]
            return np.sum(X/BHat[:,:,np.newaxis]**2,axis=(0,1))/np.sum(1/BHat**2)
        else:
            print("Unsupported rank for FSA:" + str(rank))
            raise ValueError 
            
    @property
    def mixed_col_NC_C_ratio(self):
        BHat = self.BHat
        u = self.u
        gpsipsi = self.gpsipsiHat
        FSABHat2 = self.FSA(BHat**2)
        ret = self.FSA(u**2 * BHat**2)*FSABHat2 - self.FSA(u*BHat**2)**2
        ret = ret/(self.FSA(gpsipsi/BHat**2) * self.FSA(BHat**2))
        return ret

    def mixed_col_DCni(self,ion_index = 0,impurity_index = -1):
        """The bulk ion gradient drive for the classical radial impurity flux
        Defined by Gamma_z^C = -n_z D^C_{n_i} d ln(ni)/dr"""
        impurity_index = -1
        BHat = self.BHat
        gpsipsi = self.gpsipsiHat
        return -self.FSA(gpsipsi/BHat**2)*self.ms[ion_index]*self.ns[ion_index]*self.Ts[ion_index]/(self.charges[impurity_index]*self.tau_ab[ion_index,impurity_index] * self.ns[impurity_index]*elecharge)

    def mixed_col_DCTi(self,ion_index = 0,impurity_index = -1):
        return -self.mixed_col_DCni(ion_index,impurity_index)/2

    def mixed_col_DCnz(self,ion_index = 0,impurity_index = -1):
        return -self.mixed_col_DCni(ion_index,impurity_index)/self.input.Zs[impurity_index]
    
    def mixed_col_DCni_r(self,ion_index = 0,impurity_index = -1):
        return self.mixed_col_DCni(ion_index,impurity_index)*(self.drHat_dpsiHat*self.normalization.RBar)**2
        
    def mixed_col_DCTi_r(self,ion_index = 0,impurity_index = -1):
        return self.mixed_col_DCTi(ion_index,impurity_index)*(self.drHat_dpsiHat*self.normalization.RBar)**2
   
    def mixed_col_DCnz_r(self,ion_index = 0,impurity_index = -1):
        return self.mixed_col_DCnz(ion_index,impurity_index)*(self.drHat_dpsiHat*self.normalization.RBar)**2
    
    def mixed_col_DNCni(self,ion_index = 0,impurity_index = -1):
        """The bulk ion gradient drive for the neoclassical radial impurity flux
        Defined by Gamma_z^C = -n_z D^C_{n_i} d ln(ni)/dr"""
        impurity_index = -1
        BHat = self.BHat
        u = self.u
        FSABHat2 = self.FSA(BHat**2)
        geom = self.FSA(u**2 * BHat**2) - self.FSA(u*BHat**2)**2/FSABHat2
        return -geom*self.ns[ion_index]*self.ms[ion_index]*self.Ts[ion_index]/(self.charges[impurity_index]*self.tau_ab[ion_index,impurity_index] * self.ns[impurity_index]*elecharge)

    def mixed_col_DNCTi(self,ion_index = 0,impurity_index = -1):
        return -self.mixed_col_DNCni(ion_index,impurity_index)/2

    def mixed_col_DNCnz(self,ion_index = 0,impurity_index = -1):
        return -self.mixed_col_DNCni(ion_index,impurity_index)/self.input.Zs[impurity_index]

    def mixed_col_DNCni_r(self,ion_index = 0,impurity_index = -1):
        return self.mixed_col_DNCni(ion_index,impurity_index)*(self.drHat_dpsiHat*self.normalization.RBar)**2
        
    def mixed_col_DNCTi_r(self,ion_index = 0,impurity_index = -1):
        return self.mixed_col_DNCTi(ion_index,impurity_index)*(self.drHat_dpsiHat*self.normalization.RBar)**2
        
    def mixed_col_DNCnz_r(self,ion_index = 0,impurity_index = -1):
        return self.mixed_col_DNCnz(ion_index,impurity_index)*(self.drHat_dpsiHat*self.normalization.RBar)**2
    
    
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
        if self.geometry_loaded:
            return self.geom.s[self.rind]
        else:
            try:
                return self.outputs["psiN"][()]
            except KeyError:
                raise NotImplementedError("psiN cannot be calculated for this simulation since it does not use a .bc file and does not contain psiN in the output.")

    @property
    def theta(self):
        return self.outputs["theta"][()]

    
    @property
    def zeta(self):
        return self.outputs["zeta"][()]
    
    
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
    def drHat_dpsiHat(self):
        if not np.isnan(self.aHat):
            return self.aHat/(2*self.psiAHat*np.sqrt(self.psiN))
        elif self.outputs["dTHatdpsiHat"][()][0] != 0:
            return self.outputs["dTHatdpsiHat"][()][0]/self.outputs["dTHatdrHat"][()][0]
        elif self.outputs["dnHatdpsiHat"][()][0] != 0:
            return self.outputs["dnHatdpsiHat"][()][0]/self.outputs["dnHatdrHat"][()][0]
        elif self.outputs["dPhiHatdpsiHat"][()][0] != 0:
            return self.outputs["dPhiHatdpsiHat"][()][0]/self.outputs["dPhiHatdrHat"][()][0]
        else:
            raise ValueError("Cannot calculate dr^/dpsi^ since all gradients are zero!")
        
    @property
    def dnHatdrHats(self):
        return  self.dnHatdpsiHats/self.drHat_dpsiHat

    @property
    def dTHatdrHats(self):
        return  self.dTHatdpsiHats/self.drHat_dpsiHat
    
    @property
    def dPhiHatdrHat(self):
        return  self.dPhiHatdpsiHat/self.drHat_dpsiHat
    
    @property
    def dnHatdpsiHats(self):
        if self.input.inputRadialCoordinateForGradients == 0:
            conversion_factor = 1.0
        elif self.input.inputRadialCoordinateForGradients == 1:
            conversion_factor = 1/self.psiAHat
        elif (self.input.inputRadialCoordinateForGradients == 2) or (self.input.inputRadialCoordinateForGradients == 4):
            conversion_factor = self.aHat/(2*self.psiAHat*np.sqrt(self.psiN))
        elif self.input.inputRadialCoordinateForGradients == 3:
            conversion_factor = 1/(2*self.psiAHat*np.sqrt(self.psiN))
        else:
            raise ValueError("inputRadialCoordinate should be 0,1,2,3,4; it is" + str(self.input.inputRadialCoordinate))
        return conversion_factor * self.input.dnHatdss

    @property
    def dTHatdpsiHats(self):
        if self.input.inputRadialCoordinateForGradients == 0:
            conversion_factor = 1.0
        elif self.input.inputRadialCoordinateForGradients == 1:
            conversion_factor = 1/self.psiAHat
        elif (self.input.inputRadialCoordinateForGradients == 2) or (self.input.inputRadialCoordinateForGradients == 4):
            conversion_factor = self.aHat/(2*self.psiAHat*np.sqrt(self.psiN))
        elif self.input.inputRadialCoordinateForGradients == 3:
            conversion_factor = 1/(2*self.psiAHat*np.sqrt(self.psiN))
        else:
            raise ValueError("inputRadialCoordinate should be 0,1,2,3,4; it is" + str(self.input.inputRadialCoordinate))
        return self.input.dTHatdss * conversion_factor

    @property
    def dPhiHatdpsiHat(self):
        if self.input.inputRadialCoordinateForGradients == 0:
            conversion_factor = 1.0
        elif self.input.inputRadialCoordinateForGradients == 1:
            conversion_factor = 1/self.psiAHat
        elif (self.input.inputRadialCoordinateForGradients == 2) or (self.input.inputRadialCoordinateForGradients == 4):
            conversion_factor = self.aHat/(2*self.psiAHat*np.sqrt(self.psiN))
        elif self.input.inputRadialCoordinateForGradients == 3:
            conversion_factor = 1/(2*self.psiAHat*np.sqrt(self.psiN))
        else:
            raise ValueError("inputRadialCoordinate should be 0,1,2,3,4; it is" + str(self.input.inputRadialCoordinate))
        return self.input.dPhiHatds * conversion_factor

    
    @property
    def ErHat(self):
        if self.input.inputRadialCoordinateForGradients == 0:
            conversion_factor = (2*self.psiAHat*np.sqrt(self.psiN))/self.aHat
        elif self.input.inputRadialCoordinateForGradients == 1:
            conversion_factor = (2*np.sqrt(self.psiN))/self.aHat
        elif (self.input.inputRadialCoordinateForGradients == 2) or (self.input.inputRadialCoordinateForGradients == 4):
            conversion_factor = 1.0
            #conversion_factor = self.aHat/(2*self.psiAHat*np.sqrt(self.psiN))
        elif self.input.inputRadialCoordinateForGradients == 3:
            conversion_factor = 1/self.aHat
        else:
            raise ValueError("inputRadialCoordinate should be 0,1,2,3,4; it is" + str(self.input.inputRadialCoordinate))
        return -self.input.dPhiHatds * conversion_factor

    
    @property
    def integerToRepresentTrue(self):
        return self.outputs["integerToRepresentTrue"]

    @property
    def includePhi1(self):
        return (self.outputs["includePhi1"][()] == self.integerToRepresentTrue)

    @property
    def Phi1Hat(self):
        if self.includePhi1:
            return self.outputs["Phi1Hat"][:,:,-1]
        else:
            prefactor=self.input.alpha * np.sum((self.Zs**2)*(self.nHats/self.THats))
            charge_perturbation = np.sum(self.Zs*self.n1Hat,axis=2)
            return charge_perturbation/prefactor

    @property
    def n1Hat(self):
        # zeta,theta,species
        return self.outputs["densityPerturbation"][:,:,:,-1]

    @property
    def n1(self):
        # zeta,theta,species
        return self.n1Hat * self.normalization.nBar

    
    @property
    def rmsN1Hat(self):
        n1 = self.total_nHat
        FSAn1 = self.FSA(n1)
        return np.sqrt(self.FSA((n1 - FSAn1)**2))/FSAn1

    @property
    def maxminN1Hat(self):
        n1 = self.total_nHat
        nmax = np.max(n1,axis=(0,1))
        nmin = np.min(n1,axis=(0,1))
        FSAn1 = self.FSA(n1)
        return (nmax - nmin)/FSAn1
    
    @property
    def total_nHat(self):
        # zeta,theta,species
        if self.includePhi1:
            return self.outputs["totalDensity"][:,:,:,-1]
        else:
            return self.outputs["totalDensity"][:,:,:,-1] + self.n1Hat_adiabatic


    @property
    def total_n(self):
        # zeta,theta,species
        return self.total_nHat * self.normalization.nBar


    @property
    def n1_adiabatic(self):
        return self.n1Hat_adiabatic * self.normalization.nBar

        
            
    @property
    def n1Hat2(self):
        # rms n1Hat
        # close to zero for small density assymmetry
        n1 = self.n1Hat
        FSAn1 = self.FSA(n1)
        return np.sqrt(self.FSA((n1 - FSAn1)**2))

    @property
    def total_Zeff(self):
        # exclude electrons from sum if electrons are included in the simulation
        noe = np.where(self.Zs!=-1)[0]
        n = self.total_nHat # total
        ne = np.sum(self.Zs[noe] * self.nHats[noe]) # zeroth order
        ret = np.sum(self.Zs[noe]**2 * n[:,:,noe],axis=2)/ne
        return ret
        
    @property
    def rms_total_Zeff(self):
        # rms Zeff calculated from n1Hat
        n1 = self.total_Zeff
        FSAn1 = self.FSA(n1)
        return np.sqrt(self.FSA((n1 - FSAn1)**2))

    
    @property
    def n1Hat_adiabatic(self):
        if self.includePhi1:
            return self.nHats * (np.exp(-self.Zs * self.input.alpha * self.Phi1Hat[:,:,np.newaxis]/self.THats) - 1)
        else:
            return -self.nHats * self.Zs * self.input.alpha * self.Phi1Hat[:,:,np.newaxis]/self.THats

    @property
    def boltzmann_factor(self):
        return self.Zs * self.input.alpha * self.Phi1Hat[:,:,np.newaxis]/self.THats

    
    @property
    def externalNHat(self):
        try:
            return self.outputs["externalN"][:,:,:]
        except KeyError:
            return 0.0

    
    @property
    def externalFlow(self):
        try:
            return self.outputs["externalFlow"][:,:,:]
        except KeyError:
            return 0.0

    @property
    def externalCurrent(self):
        if self.externalZs is None:
            print(self.dirname)
            externalZs = 0.0
        else:
            externalZs = self.externalZs
        
        return np.sum(externalZs * self.externalFlow,axis=2)

        
    @property
    def externalFlow_ms(self):
        return self.externalFlow * self.normalization.nBar * self.normalization.vBar

    @property
    def externalCurrent_Am2(self):
        return self.externalCurrent * self.normalization.nBar * self.normalization.vBar * self.normalization.eBar

    
    @property
    def externalCurrent_kAm2(self):
        return self.externalCurrent_Am2/1e3
     
    @property
    def externalFSABFlow(self):
        try:
            return self.outputs["externalFSABFlow"][:]
        except KeyError:
            return 0.0

    @property
    def externalZs(self):
        try:
            return self.outputs["externalZs"][:]
        except KeyError:
            return None
        
    @property
    def externalFSABjHat(self):
        if self.externalZs is None:
            print(self.dirname)
            externalZs = 0.0
        else:
            externalZs = self.externalZs
        return np.sum(externalZs * self.externalFSABFlow)

    @property
    def externalFSABjHatOverB0(self):
        return self.externalFSABjHat/self.B00Hat

    @property
    def externalFSABjOverB0_Am2(self):
        return self.externalFSABjHatOverB0 * self.normalization.eBar * self.normalization.nBar * self.normalization.vBar

    @property
    def externalFSABjOverB0_kAm2(self):
        return self.externalFSABjOverB0_Am2/1e3
    
    @property
    def FSAB2(self):
        return self.FSA(self.BHat**2)

    @property
    def B0Hat(self):
        return 2*np.fabs(self.psiAHat)/self.aHat**2
        
    @property
    def FSABExternalFlowOverRootFSAB2(self):
        return self.externalFSABFlow/(np.sqrt(self.FSAB2))#*self.FSA(self.externalNHat))

    @property
    def FSABExternalFlowOverRootFSAB2_sm2(self):
        return self.FSABExternalFlowUsingFSADensityOverRootFSAB2 * self.normalization.vBar * self.normalization.nBar

    @property
    def FSABExternalFlow_arturo(self):
        return self.externalFSABFlow * self.B0Hat/self.FSAB2

    @property
    def FSABFlow_arturo(self):
        return (self.FSABFlow * self.B0Hat/self.FSAB2)

    @property
    def FSABExternalFlow_arturo_sm2(self):
        return  self.FSABExternalFlow_arturo * self.normalization.vBar * self.normalization.nBar

    @property
    def FSABFlow_arturo_sm2(self):
        return  self.FSABFlow_arturo * self.normalization.vBar * self.normalization.nBar
    
    @property
    def externalNHat2(self):
        n1 = self.externalNHat
        FSAn1 = self.FSA(n1)
        print(FSAn1)
        return np.sqrt(self.FSA((n1 - FSAn1)**2))

         
    @property
    def VPrimeHat(self):
        return self.outputs["VPrimeHat"][()]

    @property
    def GammaHat(self):
        if self.includePhi1:
            ret=self.outputs["particleFlux_vd_psiHat"][:,-1]
        else:
            ret=self.outputs["particleFlux_vm_psiHat"][:,-1]
        return ret

    @property
    def GammaHat_rHat(self):
        if self.includePhi1:
            ret=self.outputs["particleFlux_vd_rHat"][:,-1]
        else:
            ret=self.outputs["particleFlux_vm_rHat"][:,-1]
        return ret

    # @property
    # def drHat_dpsiHat(self):
    #     tol = 1e-10
    #     if np.fabs(self.input.dTHatdss[0]) > tol:
    #         try:
    #             return self.outputs["dTHatdpsiHat"]/self.outputs["dTHatdrHat"]
    #         except KeyError:
    #             pass
    #     if np.fabs(self.input.dnHatdss[0]) > tol:
    #         try:
    #             return self.outputs["dnHatdpsiHat"]/self.outputs["dnHatdrHat"]
    #         except KeyError:
    #             pass
    #     if np.fabs(self.input.dPhiHatds) > tol:
    #         try:
    #             return self.outputs["dPhiHatdpsiHat"]/self.outputs["dPhiHatdrHat"]
    #         except KeyError:
    #             pass
        
      
            
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
    def GammaHat_C_rHat(self):
        return self.GammaHat_C * self.drHat_dpsiHat
        
    @property
    def QHat(self):
        if self.includePhi1:
            return self.outputs["heatFlux_vd_psiHat"][:,-1]
        else:
            return self.outputs["heatFlux_vm_psiHat"][:,-1]

    @property
    def QHat_rHat(self):
        if self.includePhi1:
            return self.outputs["heatFlux_vd_rHat"][:,-1]
        else:
            return self.outputs["heatFlux_vm_rHat"][:,-1]

    @property
    def ms(self):
        return self.input.mHats * self.normalization.mBar

    @property
    def Zs(self):
        return self.input.Zs
    
    
    @property
    def charges(self):
        return self.input.Zs * self.normalization.eBar
        
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
    def QHat_C_rHat(self):
        return self.QHat_C * self.drHat_dpsiHat
            
    
    @property
    def r_eff(self):
        return self.rN * self.aHat * self.normalization.RBar

        
    @property
    def dndrs(self):
        return self.dnHatdrHats  * self.normalization.nBar/self.normalization.RBar
    
    @property
    def dTdrs(self):
        return  self.dTHatdrHats  * self.normalization.nBar/self.normalization.RBar
    
    @property
    def dVdr(self):
        return  (self.VPrimeHat/self.drHat_dpsiHat) * self.normalization.RBar**2

    @property
    def Er(self):
        return self.ErHat * self.normalization.PhiBar/self.normalization.RBar

    @property
    def dPhidr(self):
        return self.dPhiHatdrHat * self.normalization.PhiBar/self.normalization.RBar

    
    @property
    def Gamma(self):
        return self.GammaHat * self.normalization.nBar*self.normalization.vBar/self.normalization.RBar

    @property
    def Gamma_s(self):
        """Gamma in units of 1/s"""
        return self.VPrimeHat * self.GammaHat * self.normalization.nBar*self.normalization.vBar*self.normalization.RBar**2

    
    @property
    def Gamma_sm2(self):
        """Gamma in units of 1/(sm^2)"""
        return self.Gamma_s/self.dVdr
    
    @property
    def Gamma_C(self):
        return self.GammaHat_C * self.normalization.nBar*self.normalization.vBar/self.normalization.RBar

    @property
    def Gamma_C_s(self):
        """Gamma_C in units of 1/s"""
        return self.VPrimeHat * self.GammaHat_C * self.normalization.nBar*self.normalization.vBar*self.normalization.RBar**2

    
    @property
    def Gamma_C_sm2(self):
        """Gamma_C in units of 1/(sm^2)"""
        return self.Gamma_C_s/self.dVdr

    @property
    def Gamma_r(self):
        return self.GammaHat_rHat * self.normalization.nBar*self.normalization.vBar/self.normalization.RBar

    @property
    def Gamma_C_r(self):
        return self.GammaHat_C_rHat * self.normalization.nBar*self.normalization.vBar/self.normalization.RBar

    @property
    def Q(self):
        """Fixed normalization"""
        return self.QHat *self.normalization.nBar*self.normalization.vBar**3*self.normalization.mBar/self.normalization.RBar

    
    @property
    def Q_W(self):
        """Q in units of W"""
        return self.VPrimeHat * self.QHat *self.normalization.nBar*self.normalization.vBar**3*self.normalization.mBar*self.normalization.RBar**2


    @property
    def Q_MW(self):
        """Q in units of MW"""
        return self.Q_W/1e6


    @property
    def Q_C(self):
        """Fixed normalization"""
        return self.QHat_C *self.normalization.nBar*self.normalization.vBar**3*self.normalization.mBar/self.normalization.RBar

    
    @property
    def Q_C_W(self):
        """Q_C in units of W"""
        return self.VPrimeHat * self.QHat_C *self.normalization.nBar*self.normalization.vBar**3*self.normalization.mBar*self.normalization.RBar**2

    
    @property
    def Q_C_MW(self):
        return self.Q_C_W/1e6
    @property
    def Q_r(self):
        """Fixed normalization"""
        return self.QHat_rHat *self.normalization.nBar*self.normalization.vBar**3*self.normalization.mBar

    @property
    def Q_C_r(self):
        """Fixed normalization"""
        return self.QHat_C_rHat *self.normalization.nBar*self.normalization.vBar**3*self.normalization.mBar

    @property
    def FSABVelocityUsingFSADensityOverRootFSAB2_ms(self):
        #print(self.FSABVelocityUsingFSADensityOverRootFSAB2)
        return self.FSABVelocityUsingFSADensityOverRootFSAB2 * self.normalization.vBar

        
    @property
    def FSABExternalVelocityUsingFSADensityOverRootFSAB2_ms(self):
        return self.FSABExternalVelocityUsingFSADensityOverRootFSAB2 * self.normalization.vBar
    
    @property
    def FSABFlow(self):
        return self.outputs["FSABFlow"][:,-1]

    @property
    def flow(self):
        return self.outputs["flow"][:,:,:,-1]

    @property
    def flow_ms(self):
        return self.flow * self.normalization.nBar * self.normalization.vBar

    @property
    def FSABVelocityUsingFSADensityOverRootFSAB2(self):
        return self.outputs["FSABVelocityUsingFSADensityOverRootFSAB2"][:,-1]
    
    @property
    def nuBar(self):
        # Inserting the definition of nu_n, this is
        # sqrt(2) * nBar * e**4 * lnLambda/(12*pi**(3/2) * sqrt(mBar)*eps0**2 * TBar**(3/2))
        # multiplying this with Z[a]**2 * Z[b]**2 * nHat[b]/(sqrt(mHat[a]) *THat[a]**(3/2))
        # gives 1/tau_{ab} as defined in my classical paper and my JPP paper.
        return self.input.nu_n*self.normalization.vBar/self.normalization.RBar

    @property
    def tau_ab(self):
        one_over_tau = np.zeros((self.input.Nspecies,self.input.Nspecies))
        for a in range(self.input.Nspecies):
            for b in range(self.input.Nspecies):
                one_over_tau[a,b] = self.nuBar*self.input.Zs[a]**2 * self.input.Zs[b]**2 * self.input.nHats[b]/(np.sqrt(self.input.mHats[a]) *self.input.THats[a]**(3/2))
        return 1/one_over_tau

    @property
    def Ts(self):
        return self.input.THats * self.normalization.TBar

    @property
    def ns(self):
        return self.input.nHats * self.normalization.nBar

    @property
    def dPhiHatdpsiHat(self):
        return self.outputs["dPhiHatdpsiHat"][()]

    
    @property
    def dPhiHatdrHat(self):
        return self.outputs["dPhiHatdrHat"][()]

    
    @property
    def dTHatdpsiHat(self):
        return self.outputs["dTHatdpsiHat"][()]

    @property
    def dnHatdpsiHat(self):
        return self.outputs["dnHatdpsiHat"][()]

    @property
    def nHats(self):
        return self.outputs["nHats"][()]

    @property
    def THats(self):
        return self.outputs["THats"][()]

    @property
    def FSABjHatOverB0(self):
        return self.outputs["FSABjHatOverB0"][()][-1]

    @property
    def FSABjOverB0_Am2(self):
        return self.FSABjHatOverB0 * self.normalization.eBar * self.normalization.nBar * self.normalization.vBar

    @property
    def FSABjOverB0_kAm2(self):
        return self.FSABjOverB0_Am2/1000
    
    
    @property
    def A1(self):
        ret = self.dnHatdpsiHats/self.nHats + self.dTHatdpsiHats/self.THats + self.input.alpha * self.input.Zs * self.dPhiHatdpsiHat/self.THats
        return ret

    @property
    def A2(self):
        ret = self.dTHatdpsiHats/self.THats
        return ret

    
    @property
    def A1_r(self):
        print(self.drHat_dpsiHat)
        return self.A1/(self.drHat_dpsiHat*self.normalization.RBar)
        
    @property
    def A2_r(self):
        return self.A2/(self.drHat_dpsiHat*self.normalization.RBar)

    @property
    def deltaT(self):
        if self.input.geometryScheme == 11:
            B00Hat = self.B00Hat
            IHat = self.IHat
            GHat = self.GHat
            iota = self.iota
        RHat = np.fabs((GHat + iota*IHat)/B00Hat)
        LT =self.THats/self.dTHatdrHats
        return  RHat/LT

    @property
    def deltaN(self):
        if self.input.geometryScheme == 11:
            B00Hat = self.B00Hat
            IHat = self.IHat
            GHat = self.GHat
            iota = self.iota
        RHat = np.fabs((GHat + iota*IHat)/B00Hat)
        LN =self.nHats/self.dnHatdrHats
        return  RHat/LN
        
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

    @property
    def collisionality_ab(self):
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
        nuHats = np.zeros((Nspecies,Nspecies))
        for a in range(Nspecies):
            for b in range(Nspecies):
                nuHats[a,b] =(Zs[b]**2*nHats[b]) * (Zs[a]/THats[a])**2 
        nuHats = nuHats * nu_n *RHat
        return nuHats

    @property
    def jrHat(self):
        return np.sum(self.input.Zs*self.GammaHat)
            
if __name__=="__main__":

    norm_filename = "norm.namelist"
    species_filename = "species"
    simul = Sfincs_simulation('.')

    print("GammaHat:")
    print(simul.GammaHat)
    print("GammaHat_C:")
    print(simul.GammaHat_C)
    print("QHat:")
    print(simul.QHat)
    print("QHat_C:")
    print(simul.QHat_C)
    print("FSABFlow:")
    print(simul.FSABFlow)

    print("Gamma_C/Gamma:")
    print(simul.GammaHat_C/simul.GammaHat)
    print("Q_C/Q:")
    print(simul.QHat_C/simul.QHat)

