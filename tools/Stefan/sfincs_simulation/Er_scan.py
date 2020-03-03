from __future__ import division
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

import os


from sfincs_simulation import Sfincs_simulation

class Er_scan(object):

    def __getattr__(self, name):
        """Get attributes from the underlying simulations by interpolation"""
        # REMEMBER: we cannot acccess attributes of this class with dot notation
        # use self.__dict__["variable name"] instead.
        y = [getattr(s,name) for s in self.__dict__["simuls"]]
        x = self.__dict__["Ers"]
        return np.array([self.__dict__["interpolator"](x,y)(r) for r in self.__dict__["roots"]])

    def __init__(self,dirname, interpolator =PchipInterpolator ,input_name="input.namelist",norm_name="norm.namelist",species_name="species",override_geometry_name=None, load_geometry=False):
        """ dirname: The directory where the Er scan is located in. 
                     This directory will be scanned for subdirs. """

        # TODO: may not work on absolute paths. gotta test
        self.subdirs = [name for name in os.listdir(dirname) if os.path.isdir(os.path.join(dirname, name))]
        
        # sort simulations based on Ers
        simuls = [Sfincs_simulation("./" + dirname + "/" + subdir, input_name,norm_name,species_name,override_geometry_name, load_geometry) for subdir in self.__dict__["subdirs"]]
        #inputRadialCoordinateForGradients = np.array([s.input.inputRadialCoordinateForGradients for s in simuls])
        # NOTE: this assumes all the simulations use the same inputRadialCoordinateForGradients. TODO: remove this assumption.
        Ers = np.array([s.Er for s in simuls])
        tmp = sorted(zip(Ers,simuls))
        self.simuls = [b for a,b in tmp]
        self.Ers = [a for a,b in tmp]
        self.interpolator = interpolator
        
        jrs = np.array([s.jrHat for s in self.__dict__["simuls"]])
        maxjr = np.max(jrs)
        minjr = np.min(jrs)
        if not ((maxjr > 0) and (minjr < 0)):
            raise ValueError("Ambipolar radial electric field cannot be found by interpolation. maxjr = " + str(maxjr) + ", minjr = " + str(minjr))
            # TODO: give estimate of ambipolar Er value based on extrapolation
        else:
            self.roots = self.solve_for_ambipolar_Er(self.Ers,jrs)

        self.root_types = self.classify_roots(self.__dict__["roots"])


    def solve_for_ambipolar_Er(self,Ers,jrs,NEr_fine = 500):
        # code taken from sfincsScanPlot_2
        # assumes Ers is sorted
        Ermin = Ers[0]
        Ermax = Ers[-1]
        interpolator = self.__dict__["interpolator"](Ers,jrs)
        Er_fine = np.linspace(Ermin, Ermax, num=NEr_fine)
        radialCurrent_fine = interpolator(Er_fine)
        # Approximately find points where the sign of radial current flips:
        positiveCurrent = (radialCurrent_fine>0)
        signFlips = (positiveCurrent[:-1] != positiveCurrent[1:])
        numRoots = sum(1 for x in signFlips if x)
        #print "Number of roots found for E_r: ",numRoots
        roots = []
        #for i in range(numRoots):
        for index,value in enumerate(signFlips):
            if value:
                roots.append(brentq(interpolator,Er_fine[index],Er_fine[index+1]))
                # Convert standard array to numpy array:
                roots = np.sort(np.array(roots))

        return roots


    def classify_roots(self,roots):
        # Classify roots as ion, unstable, or electron root
        # Code taken from sfincsScanPlot_2
        if len(roots)==1:
            if roots[0] >0:
                root_types=['electron']
            else:
                root_types=['ion']
        elif len(roots)==3:
                root_types=['ion','unstable','electron']
        else:
            # If the number of roots is not 1 or 3, then flag the case as complicated:
            root_types=['unknown']*len(roots)

        return root_types
