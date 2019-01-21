from __future__ import division
import numpy as np
import f90nml
import h5py

import subprocess
from warnings import warn

from shutil import copytree

from geomlib import bcgeom
from fluxcoorddiscr import fluxcoorddiscr



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

def create_species(normalization,filename="species",database_filename = "/home/bstefan/documents/svn/stefan/scripts/python/PERFECT_utils/species_database.namelist"):
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
                                   },
                'speciesparameters' : {'zs' : [1.0],
                                       'mhats' : [1.0],
                                       'nhats' : [1.0],
                                       'thats' : [1.0],
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
    def nu_n(self):
        return self.get_value_from_input_or_defaults("physicsparameters","nu_n")

    @property
    def Zs(self):
        return self.get_value_from_input_or_defaults("speciesparameters","zs")

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
            
    def __str__(self):
        return str(self.input_name)

##############################

class Sfincs_simulation(object):

    def absolute_path(self,relative_path):
        return self.dirname + "/" + relative_path

    def copy(self,new_dirname):
        copytree(self.dirname,new_dirname)
        return Sfincs_simulation(new_dirname,input_name=self.input_name,norm_name=self.norm_name,species_name=self.species_name)
        
    def __init__(self,dirname,input_name="input.namelist",norm_name="norm.namelist",species_name="species"):
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
                self.geom = bcgeom(self.absolute_path(self.input.equilibrium_name),self.min_Bmn,self.max_m,self.maxabs_n,self.symmetry,self.signcorr,verbose)
            except IOError as e:
                raise IOError("Error loading geometry: " + str(e))
            
            self.Nperiods = self.geom.Nperiods
            self.psiA = self.geom.psi_a/(self.normalization.BBar*self.normalization.RBar**2)
            self.aHat = self.geom.minorradiusVMEC/self.normalization.RBar
            
            
            if self.input.inputRadialCoordinate == 0:
                self.rN_wish = np.sqrt(self.input.coordinate_wish/self.psiAHat)
            elif self.input.inputRadialCoordinate == 1:
                self.rN_wish = np.sqrt(self.input.coordinate_wish)
            elif (self.input.inputRadialCoordinate == 2) or (self.input.inputRadialCoordinate == 4):
                self.rN_wish = self.input.coordinate_wish/self.aHat 
            elif self.input.inputRadialCoordinate == 3:
                self.rN_wish = self.input.coordinate_wish
            else:
                raise ValueError("inputRadialCoordinate should be 0,1,2,3,4; it is" + str(self.input.inputRadialCoordinate))
            self.rN = self.geom.rnorm
            self.rind=np.argmin(np.fabs(self.rN - self.rN_wish))
            self.iota = self.geom.iota[self.rind]
            
            self.Booz = fluxcoorddiscr(self.geom,self.rind,self.input.Ntheta,self.input.Nzeta,u_zeroout_Deltaiota=self.zeroout_Deltaiota,name='Boozer')
            self.GHat=self.Booz.G/(self.normalization.BBar*self.normalization.RBar)
            self.IHat=self.Booz.I/(self.normalization.BBar*self.normalization.RBar)
            self.B00Hat = self.Booz.B00/(self.normalization.BBar)
            
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
            return self.outputs["particleFlux_vd_psiHat"][()]

    @property
    def QHat(self):
        if self.includePhi1:
            return self.outputs["heatFlux_vd_psiHat"][:,-1]
        else:
            return self.outputs["heatFlux_vd_psiHat"][()]

    @property
    def FSABFlow(self):
        if self.includePhi1:
            return self.outputs["FSABFlow"][:,-1]
        else:
            return self.outputs["FSABFlow"][()]

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

            
            
if __name__=="__main__":

    norm_filename = "norm.namelist"
    species_filename = "species"
    simul = Sfincs_simulation('.')

    print "GammaHat:"
    print simul.GammaHat
    print "QHat:"
    print simul.QHat
    print "FSABFlow:"
    print simul.FSABFlow


