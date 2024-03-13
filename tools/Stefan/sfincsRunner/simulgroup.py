from __future__ import division, print_function

import subprocess
import os
from shutil import copy

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from sfincs_simulation import Sfincs_simulation

from sfincs_utils import moveScan, scanInDir, modifyScan1, modifyScan2, changeVar
from simulation import Simulation

def subdirs(dir):
    subdirs = [os.path.join(dir, f) for f in os.listdir(dir) if os.path.isdir(os.path.join(dir, f))]
    return subdirs

class Simulgroup(object):
    """Object for handling a dynamic group of SFINCS simulations required to arrive at a converged ambipolar electric field.

    Creates simulation in a directory to find parameters for numerical convergence and then the ambipolar electric field. 

    States:
    None: look for previous simulations and determine next state.  [default state]
    0: initial resolution scan
    1: Er to find more ambipolar solution
    2: final Er scan
    3: more ambipolar resolution scan
    9: Done [terminating state]"""
    
    params = ["Ntheta","Nzeta","Nxi","Nx"]
    NErScan2Simuls = 3
    archivedir1 = "OLD_scan1" #used when going from 0 to 1
    archivedir2 = "scan1" # used when going from [0 or 3] to 2
    archivedir3 = "OLD_scan2" #used when going from 1 to 3
    def __init__(self,dirname, state = None):
        self.upper = 0.0
        self.lower = 0.3
        self.dirname = dirname
        self.simullist = []
        self.scan = None
        self.auto = True
        self.setup = False # turned True after first 'run' call

        if state is None:
            # DEFAULT: 
            # determine initial state of simulation by looking at
            # the old simulations present in the directory
            self.determine_initial_state()
        else:
            # NON-DEFAULT:
            # manual override of initial state
            self.state = state
        #print(self.dirname + " : initial state:" + str(self.state))
        
    def determine_initial_state(self):
        # look at existing archive dirs to retrace steps.
        # 1st: look for OLD_scan1
        # 2nd: look for OLD_scan2
        # 3rd: look for scan1
        # 4th (TODO?): maybe look at current state in main directory?
        if os.path.isdir(os.path.join(self.dirname, type(self).archivedir1)):
            self.simullist = self.simullist +  [Simulation.fromDir(d) for d in subdirs(self.dirname + "/" + type(self).archivedir1)]
            if os.path.isdir(os.path.join(self.dirname, type(self).archivedir3)):
                state = 3
                self.simullist = self.simullist +  [Simulation.fromDir(d) for d in subdirs(self.dirname + "/" + type(self).archivedir3)]
            else:
                state = 1
        if os.path.isdir(os.path.join(self.dirname, type(self).archivedir2)):
            self.simullist = self.simullist +  [Simulation.fromDir(d) for d in subdirs(self.dirname + "/" + type(self).archivedir2)]
            state = 2
        else:
            state = 0
        
        self.state = state
        return state
        
    def setup_initial_state(self):
        """Run once at the start of the first run call"""
        if self.setup == False:
            if (self.state  == 0) or (self.state == 3):
                modifyScan1(self.dirname,N=2,lower=self.lower,upper=self.upper)
                if self.state == 0:
                    self.params = type(self).params
            elif (self.state  == 1) or (self.state == 2):

                modifyScan2(self.dirname,N=type(self).NErScan2Simuls,lower=0.3,upper=0.3)
            existing_dirnames = [s.dirname for s in self.simullist]
            for d in subdirs(self.dirname):
                if d[-1] == "/":
                    d = d[:-1]
                if d.rsplit("/",1)[-1] not in ["OLD_scan1","OLD_scan2","scan1"]:
                    if d not in existing_dirnames:
                        self.simullist.append(Simulation.fromDir(d))
            self.startNewSimuls()
        self.setup = True

    def to_state0(self):
        # this will never be called
        # since we cannot transition back to state0
        # and initial state is handled by self.setup_initial_state()
        pass
    
    def to_state1(self):
        newdir = self.dirname + "/" + type(self).archivedir1
        self.archiveScan(newdir)
        modifyScan2(self.dirname,N=type(self).NErScan2Simuls,lower=0.3,upper=0.3)
        # look at convergence in old scan
        # and change parameters
        tmpscan = Scan1(self.dirname + "/" + type(self).archivedir1)
        conv,convparams = tmpscan.check_convergence(self.params)
        params_to_scan = {}
        for ii,p in enumerate(self.params):
            if conv[ii]:
                changeVar(self.dirname,"resolutionParameters", p,convparams[ii])
        self.startNewSimuls()
        self.state = 1
        
    def to_state2(self):
        newdir = self.dirname + "/" + type(self).archivedir2
        self.archiveScan(newdir)
        # check scan1 to get converged values
        tmpscan = Scan1(self.dirname + "/" + type(self).archivedir2)
        conv,convparams = tmpscan.check_convergence(self.params)
        for ii,p in enumerate(self.params):
            changeVar(self.dirname,"resolutionParameters", p,convparams[ii])
        # setup Er scan
        modifyScan2(self.dirname,N=type(self).NErScan2Simuls,lower=0.3,upper=0.3)
        self.startNewSimuls()
        self.state = 2
   
    
    def to_state3(self):
        self.upper = 0.0
        self.lower = 0.3
        # this is only called when going from state 1 to state 3
        Er = self.scan.roots
        newdir = self.dirname + "/" + type(self).archivedir3
        self.archiveScan(newdir)
        # check old Scan1 to get previously converged parameters
        # to only scan previously unconverged parameters
        tmpscan = Scan1(self.dirname + "/" + type(self).archivedir1)
        conv,convparams = tmpscan.check_convergence(type(self).params)
        params_to_scan = {}
        for ii,p in enumerate(type(self).params):
            if not conv[ii]:
                params_to_scan[p] = 2
        self.params = list(params_to_scan.keys())
        changeVar(self.dirname,"physicsParameters", "dPhiHatdrHat",Er)
        modifyScan1(self.dirname,N=0,lower=self.lower,upper=self.upper,**params_to_scan)
        self.startNewSimuls()
        self.state = 3
        

    def run(self,force=False):
        """"Checks status of simulations and automatically launches new jobs as needed."""
        if self.setup == False:
            self.setup_initial_state()
        for job in self.scan.jobs:
            # relaunch individual jobs if needed (OOM or TIME)
            # unless force is True, 
            # in which case they are always restarted.
            if job.auto:
                job.run(force=force)
        if self.scan.eval_status() == "DONE":
            # all jobs in scan done, check if scan converged
            # exact checks performed will depend on our state
            if self.state == 0 or self.state == 3:
                conv,convparams = self.scan.check_convergence(self.params)
                if not all(conv):
                    # "if not converged, check ambipolar"
                    ambi = self.scan.check_ambipolarity(0.1)
                    if not ambi:
                        # "if not ambpolar, do Er scan (state 1)"
                        self.to_state1()
                    else:
                        # add extra points to scan
                        # will try to add more points
                        # until at least 1 new simulation is started
                        while True:
                            self.upper = self.upper + 0.3
                            params_to_scan = {}
                            for ii,p in enumerate(self.params):
                                if not conv[ii]:
                                    params_to_scan[p] = 2
                                    modifyScan1(self.dirname,N=0,lower=self.lower,upper=self.upper,**params_to_scan)
                            self.startNewSimuls()
                            if len(self.scan.jobs) > 0:
                                break
                else:
                    self.to_state2()
            elif self.state == 1 or self.state == 2:
                Er = self.scan.roots
                if Er is None:
                    alpha = 0.1
                    newEr = self.scan.extrapolate_for_root()
                    if newEr is not None:
                        if newEr >= self.scan.maxEr:
                            l = newEr - alpha * (newEr - self.scan.maxEr)
                            u = newEr + alpha * (newEr - self.scan.maxEr)
                        elif newEr < self.scan.minEr:
                            l = newEr - alpha * (self.scan.minEr - newEr)
                            u = newEr + alpha * (self.scan.minEr - newEr)
                        else:
                            # this should be impossible
                            raise ValueError("Impossible!")
                        modifyScan2(self.dirname,N=type(self).NErScan2Simuls,lower=l,upper=u,absolute=True)
                        self.startNewSimuls()
                    else:
                        "Extrapolation failed?"
                        "write to log of problematic simulations and wait"
                elif self.state == 1:
                    # TODO is this good enough?
                    # perform a new resolution scan
                    self.to_state3()
                else:
                    # check if Er has datapoints sufficiently close on both sides
                    # if not, make a new scan with points close to Er added
                    Ers  = self.scan._Ers
                    if any((Ers<Er) * (Ers>Er*0.9)):
                        if any((Ers>Er) * (Ers<Er*1.1)):
                            # points close enough exist. We are done!
                            self.state = 9
                        else:
                            u = Er * 1.1
                            modifyScan2(self.dirname,N=2,lower=Er,upper=u,absolute=True)
                            self.startNewSimuls()
                            # no Er close above
                            # Er close underneath
                    else:
                        if any((Ers>Er) * (Ers<Er*1.1)):
                            l = Er*0.9
                            modifyScan2(self.dirname,N=2,lower=l,upper=Er,absolute=True)
                            self.startNewSimuls()
                            # no Er close underneath
                            # Er close above
                        else:
                            # No Er close above or underneath
                            l = Er * 0.9
                            u = Er * 1.1
                            modifyScan2(self.dirname,N=3,lower=l,upper=u,absolute=True)
                            self.startNewSimuls()


    def archiveScan(self,archivedir):
        # update the dirname of the simulations in simullist
        for s in self.simullist:
            if s.dirname in self.scan.subdirs:
                s.dirname = archivedir +"/"+ s.dirname.rsplit("/",1)[-1]
        # move the actual directories to correspond to the new names
        moveScan(self.dirname,archivedir)
        
        

    def startNewSimuls(self):
        stype, dirnames, jobIDs = scanInDir(self.dirname)
        if stype == 1:
            self.scan = Scan1(self.dirname,[Simulation(d,ID) for d,ID in zip(dirnames,jobIDs)])
        elif stype == 2:
            self.scan = Scan2(self.dirname,[Simulation(d,ID) for d,ID in zip(dirnames,jobIDs)])
        else:
            # this should never happen
            raise ValueError("Incorrect type of scan '"  + str(stype) + "' in state " + str(self.state) + ".")

        # do not add duplicates
        # TODO: test that dirnames are on the same format
        existing_dirnames = [s.dirname for s in self.simullist]
        for d,ID in zip(dirnames,jobIDs):
            if d not in existing_dirnames:
                self.simullist.append(Simulation(d,ID))
                

    @property
    def statestring(self):
        state = self.state
        if state == 0:
            ret = "conv"
        elif state == 1:
            ret = "Er 1"
        elif state == 2:
            ret = "Er 2"
        elif state == 3:
            ret = "cnv2"
        elif state == 9:
            ret = "DONE"
        else:
            ret = "????"
        return ret

    @property
    def short_dirname(self):
        if self.dirname[-1] == "/":
            dirname = self.dirname[:-1]
        else:
            dirname = self.dirname
        drs = dirname.rsplit("/",1)
        return drs[1]

    @property
    def info(self):
        return "simulgroup info"

    @property
    def typestring(self):
        return "group"


class Scan(object):
    def __init__(self,joblist,status=None):
        self.jobs = joblist
        if status is not None:
            self.status = status
        else:
            self.eval_status()
    
    def eval_status(self):
        if all([(j.status == "DONE") for j in self.jobs]):
            self.status="DONE"
        else:
            self.status="WAIT"
        return self.status

    def __getitem__(self,i):
        return self.jobs[i]

    def __len__(self):
        return len(self.jobs)

    def append(self,item):
        self.jobs.append(item)

    def remove(self,item):
        self.jobs.remove(item)

    def __str__(self):
        return str([j.dirname for j in self.jobs])
        


class Scan1(Scan):
    """Object for representing a SfincsScan of type 1 (numerical convergence).

    dirname -- directory of the scan.
    joblist -- A list of Simulation objects for the jobs in the scan.
    
    Scan1 is aware of both 'active' jobs (in the joblist) and finished simulations (present in the directory). Convergence tests are performed with both active and finished simulations."""

    def __init__(self,dirname,joblist=[],status=None):
        self.dirname = dirname
        super().__init__(joblist,status)

    @property
    def subdirs(self):
        return [os.path.join(self.dirname, f) for f in os.listdir(self.dirname) if (os.path.isdir(os.path.join(self.dirname, f)) and f not in ["OLD_scan1","OLD_scan2","scan1"])]

    def check_convergence(self, params = ["Ntheta","Nzeta","Nxi","Nx"],tol = 1e-6,convTol = 1e-1):
        Nparams = len(params)
        conv = [False] * Nparams
        convparams = [None] * Nparams
        dirs =  [[] for x in range(0,Nparams)]
        simuls =  [[] for x in range(0,Nparams)]
        # we don't use the dirs purely in the scan, but all subdirs
        # thus, older simultations from previous scans are also accounted for
        for i,p in enumerate(params):
            for d in self.subdirs:
                if d.rsplit("/",1)[1] == "baseCase":
                    cond = True
                elif p == "Nx":
                    cond = (p in d and not "Nxi" in d)
                else:
                    cond = (p in d)
                if cond:
                    dirs[i].append(d)

                    s = Sfincs_simulation(d,load_geometry=False)
                    simul = {}
                    simul["Ntheta"] = s.Ntheta
                    simul["Nzeta"] = s.Nzeta
                    simul["Nxi"] = s.Nxi
                    simul["Nx"] = s.Nx
                    simul["d"] = d
                    try:
                        simul["Gamma"] = s.GammaHat
                        simul["FSABFlow"] = s.FSABFlow
                        simul["Q"] = s.QHat
                    except KeyError:
                        # skip this simulation
                        print("Simulation '" + d + "' has no output")
                        continue
                    except AttributeError:
                        # skip this simulation
                        print("Simulation '" + d + "' has no output")
                        continue
                    simuls[i].append(simul)
        
        for i,p in enumerate(params):
            ss = simuls[i]
            # sort by the value of the param
            Nps  = [s[p] for s in ss]
            tmp = sorted(zip(Nps,ss), key = lambda x: x[0])
            ss = [b for a,b in tmp]
            Nps = [a for a,b in tmp]
            if not len(ss)>1:
                # need more than 1 simul to check convergence
                print("Need more than 1 simulation to assess convergence. Affected Parameter: '" + p + "'. Dir: '" + self.dirname + "'")
                continue
            dX = None
            D = Nps[-1]/Nps[-2] - 1.0
            if np.fabs(D) > tol:
                dX = D
            dGamma = ss[-1]["Gamma"]/ss[-2]["Gamma"] - 1.0
            dF = ss[-1]["FSABFlow"]/ss[-2]["FSABFlow"] - 1.0
            dQ = ss[-1]["Q"]/ss[-2]["Q"] - 1.0
            
            if dX is not None:
                dGdX = dGamma/dX
                dFdX = dF/dX
                dQdX = dQ/dX
            if all([(np.fabs(dGdX)< convTol).all(),(np.fabs(dFdX)< convTol).all(),(np.fabs(dQdX) < convTol).all()]):
                conv[i] = True
                convparams[i] = Nps[-2]
                
        return (conv,convparams)
        
    def check_ambipolarity(self,tol=0.2):
        for d in self.subdirs:
            if d == "baseCase":
                s = Sfincs_simulation(d,load_geometry=False)
                break
        else:
            for d in self.subdirs: 
                if "N" in d:
                    s = Sfincs_simulation(d,load_geometry=False)
                    break
            scale = max(np.fabs(s.GammaHat))
        return (s.jrHat/scale < tol)

class Scan2(Scan):    
    """Object for representing a SfincsScan of type 2 (Er scan).

    dirname -- directory of the scan.
    joblist -- A list of Job objects for the jobs in the scan.
    
    Scan2 is aware of both 'active' jobs (in the joblist) and finished simulations (present in the dirname directory). Convergence tests are performed with both active and finished simulations."""

    interpolator = PchipInterpolator
    def __init__(self,dirname,joblist=[],status=None):
        self.dirname = dirname
        self.subdirs = [os.path.join(self.dirname, f) for f in os.listdir(self.dirname) if (os.path.isdir(os.path.join(self.dirname, f)) and f not in ["OLD_scan1","OLD_scan2","scan1"])]
        super().__init__(joblist,status)

    def get_converged_simuls(self):
        simuls = []
        for subdir in self.subdirs:
            try:
                s = Sfincs_simulation(subdir,load_geometry=False)
            except:
                pass
            else:
                simuls.append(s)
        # remove unconverged simulations
        simuls = [s for s in simuls if s.converged]
        return simuls
        
    def extrapolate_for_root(self):
        simuls = self.get_converged_simuls()
        Ers = np.array([s.dPhiHatdrHat for s in simuls])
        tmp = sorted(zip(Ers,simuls), key = lambda x: x[0])
        simuls = [b for a,b in tmp]
        Ers = np.array([a for a,b in tmp])
        jrs = np.array([s.jrHat for s in simuls])
        if jrs[-1] * (Ers[-1] - Ers[-2]) < 0:
            # extrapolate up
            newEr = Ers[-1] - jrs[-1] * (Ers[-1] - Ers[-2])/(jrs[-1] - jrs[-2])
        elif jrs[0] *(Ers[1] - Ers[0]) > 0:
            # extrapolate down
            newEr = Ers[0] - jrs[0] * (Ers[1] - Ers[0])/(jrs[1] - jrs[0])
        else:
            # this is awful
            newEr = None
        return newEr


    @property
    def roots(self):
        simuls = self.get_converged_simuls()
        # sort simulations based on Er
        # NOTE: the sorting assumes all the simulations use the same inputRadialCoordinateForGradients. This should be true in all normal situations.
        
        Ers = np.array([s.dPhiHatdrHat for s in simuls])
        tmp = sorted(zip(Ers,simuls), key = lambda x: x[0])
        simuls = [b for a,b in tmp]
        Ers = [a for a,b in tmp]
        jrs = np.array([s.jrHat for s in simuls])
        maxjr = np.max(jrs)
        minjr = np.min(jrs)
        if not ((maxjr > 0) and (minjr < 0)):
            roots = None
            return roots
        else:
            roots = self.solve_for_ambipolar_Er(Ers,jrs)
            root_types = self.classify_roots(roots)
            if len(roots) == 1:
                return roots[0]
            else:
                raise ValueError("multiple roots!")

    def solve_for_ambipolar_Er(self,Ers,jrs,NEr_fine = 500):
        # code taken from sfincsScanPlot_2
        # assumes Ers is sorted
        Ermin = Ers[0]
        Ermax = Ers[-1]
        Er_fine = np.linspace(Ermin, Ermax, num=NEr_fine)
        interpolator = type(self).interpolator(Ers,jrs)
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

    @property
    def _Ers(self):
        simuls = self.get_converged_simuls()
        Ers = np.array([s.dPhiHatdrHat for s in simuls])
        return sorted(Ers)
    

    @property
    def maxEr(self):
        return self._Ers[-1]
        
    @property
    def minEr(self):
        return self._Ers[0]
        

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

