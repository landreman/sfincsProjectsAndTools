from __future__ import division, print_function

import subprocess
import os
from shutil import copy

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from sfincs_simulation import Sfincs_simulation


class OutputNotFound(Exception):
    def __init__(self, dirname):
        self.foo = dirname

def copyInput(dirname,newdirname):
    if not os.path.exists(newdirname):
        os.makedirs(newdirname)
        copy(dirname + "/input.namelist", newdirname)
        copy(dirname + "/job.sfincsScan", newdirname)


class Scan(object):
    def __init__(self,joblist,status=None):
        self.jobs = joblist
        if status is not None:
            self.status = status
        else:
            self.eval_status()
    
    def eval_status(self):
        if any([(j.status == "WAITING" or j.status == "TIME") for j in self.jobs]):
            self.status="WAITING"
        elif all([(j.status == "DONE") for j in self.jobs]):
            self.status="DONE"
    
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
    def check_convergence(self, params = ["Ntheta","Nzeta","Nxi","Nx"],tol = 1e-6,convTol = 1e-1):
        Nparams = len(params)
        conv = [False] * Nparams
        convparams = [None] * Nparams
        dirs =  [[] for x in range(0,Nparams)]
        simuls =  [[] for x in range(0,Nparams)]
        # we don't use the dirs purely in the scan, but all subdirs
        # thus, older simultations from previous scans are also accounted for
        subdirs = [os.path.join(self.dirname, f) for f in os.listdir(self.dirname) if os.path.isdir(os.path.join(self.dirname, f))]
        for i,p in enumerate(params):
            for d in subdirs:
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
                print(self.jobs)
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
        subdirs = [f for f in os.listdir(self.dirname) if os.path.isdir(os.path.join(self.dirname, f))]
        for d in subdirs:
            if d == "baseCase":
                s = Sfincs_simulation(self.dirname + "/" + d,load_geometry=False)
                break
        else:
            for d in subdirs: 
                if "N" in d:
                    s = Sfincs_simulation(self.dirname + "/" + d,load_geometry=False)
                    break
        scale = max(s.GammaHat)
        return (s.jrHat/scale < tol)

class Scan2(Scan):
    interpolator = PchipInterpolator


    def extrapolate_for_root(self):
        subdirs = [os.path.join(self.dirname, f) for f in os.listdir(self.dirname) if os.path.isdir(os.path.join(self.dirname, f))]
        simuls = [Sfincs_simulation(subdir,load_geometry=False) for subdir in subdirs]
        
        # remove unconverged simulations
        simuls = [s for s in simuls if s.converged]
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
        subdirs = [f for f in os.listdir(self.dirname) if os.path.isdir(os.path.join(self.dirname, f))]
        simuls = [Sfincs_simulation(self.dirname + "/" + subdir,load_geometry=False) for subdir in subdirs]
        
        # remove unconverged simulations
        simuls = [s for s in simuls if s.converged]
    
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
        subdirs = [f for f in os.listdir(self.dirname) if os.path.isdir(os.path.join(self.dirname, f))]
        simuls = [Sfincs_simulation(self.dirname + "/" + subdir,load_geometry=False) for subdir in subdirs]
        # remove unconverged simulations
        simuls = [s for s in simuls if s.converged]
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


class Job(object):
    donestring = "Done with the main solve."
    oomstring = "oom-kill"
    maybeoomstring = "killed"
    timestring = "time limit"

    @staticmethod
    def getLatestJobID(dirname):
        completed = subprocess.run("ls " + dirname + " | grep .out.", shell=True, capture_output=True, text=True)
        filenames = completed.stdout.split()
        if len(filenames) > 0:
            ids = [int(fn.rsplit(".",1)[-1]) for fn in filenames]
            ids.sort()
            jobid = str(ids[-1])
        else:
            jobid = None
        return jobid

    def __init__(self,dirname,jobid):
        self.dirname = dirname
        self.jobid = jobid
        if jobid is not None:
            self.stdout = dirname + "/job.out." + jobid
            self.stderr = dirname + "/job.err." + jobid
        self.jobfile = dirname + "/job.sfincsScan"

    @classmethod
    def fromDir(cls,dirname):
        jobid = Job.getLatestJobID(dirname)
        return cls(dirname,jobid)
        
    def retryJobID(self):
        self.jobid =  Job.getLatestJobID(self.dirname)

    def check_squeue(self):
        completed = subprocess.run(["squeue","-j",self.jobid], capture_output=True, text=True)
        if len(completed.stdout) > 0:
            status = "WAITING"
        else:
            status = "???"
        return status


    def getTimeLimit(self):
        with open(self.jobfile,'r') as f:
            t = f.read().rsplit("#SBATCH --time=",1)[1].split(None,1)[0]
        h,m,s = t.split(":")
        h = int(h)
        m = int(m)
        s = int(s)
        return h*3600 + m*60 + s

    def setTimeLimit(self,s):
        h = s//3600
        m = (s%3600)//60
        s = s%60 
        t = str(h).zfill(2) +  ":" + str(m).zfill(2) + ":" + str(s).zfill(2)
        new_filecontent = []
        with open(self.jobfile,'r') as f:
            for l in f.readlines():
                if "#SBATCH --time=" in l:
                    l = "#SBATCH --time=" + t + "\n"
                new_filecontent.append(l)
        with open(self.jobfile,'w') as f:
            for l in new_filecontent:
                f.write(l)
        
    @property
    def status(self):
        if self.jobid is None:
            return "MISSING JOB ID"

        if not os.path.isfile(self.stdout):
            status = self.check_squeue()
        else:
            with open(self.stdout,'r') as f:
                out = f.read()
            with open(self.stderr,'r') as f:
                err = f.read()

            if Job.donestring in out:
                status = "DONE"
            elif Job.oomstring in err.lower():
                status = "OOM"
            elif Job.timestring in err.lower():
                status = "TIME"
            elif Job.maybeoomstring in err.lower():
                # sometimes OOM results in a different error message
                status = "OOM"
            else:
                status = self.check_squeue()
        return status

def scanInDir(dirname):
    """ Launches a scan inside the named directory.
    dirname -- string with path to directory to scan in
    returns: a list of Jobs in the scan. Will also include old jobs not yet finished.
    """
    tmp = os.getcwd()
    os.chdir(dirname)
    # this breaks backwards compatibility, python3.7 and higher only!
    finished = subprocess.run("yes | sfincsScan",shell=True, capture_output=True, text=True)
    os.chdir(tmp)
    stdout = finished.stdout
    search = "Here are the directories that will be created:\n["
    i1 = stdout.find(search)
    i2 = i1 + stdout[i1:].find("]")
    subdirnames = [x.strip()[1:-1] for x in stdout[i1+len(search):i2].split(',') ] 
    dirnames = [dirname + "/" + sdn for sdn in subdirnames if len(sdn) > 0]
    jobids = []
    search = "Submitted batch job "
    lenS = len(search)
    for i in range(len(dirnames)):
        i2 = i2 + stdout[i2:].find(search) + lenS
        jobids.append(stdout[i2:].split(None,1)[0])

    stype = scanType(dirname)
    subdirs = [os.path.join(dirname, f) for f in os.listdir(dirname) if os.path.isdir(os.path.join(dirname, f))]
    # TODO: Check if this works
    olddirs = [(d,Job.getLatestJobID(d)) for d in subdirs if d not in (dirnames + ['prelim', 'Erscan'])]
    print("These old unfinished simulations were automatically added to the list of simulations:")
    print(olddirs)
    # add jobs that were previously run but ran out of time
    # or OOM
    # in the former case, nothing will actually be done
    # except that the list of unfinished simulations will be non-empty
    # so that the program waits for user intervention
    for (x,y) in olddirs:
        s = Job(x,y).status
        if s == "OOM" or s == "TIME":
            dirnames.append(x)
            jobids.append(y)

    if stype == 1:
        scan = Scan1([Job(x,y) for (x,y) in zip(dirnames,jobids)])
    elif stype == 2:
        scan = Scan2([Job(x,y) for (x,y) in zip(dirnames,jobids)])
    else:
        scan = Scan([Job(x,y) for (x,y) in zip(dirnames,jobids)])
    scan.dirname = dirname
    return scan

def scanType(dirname):
    # check which type of scan we launched
    filename = dirname + "/input.namelist"
    with open(filename) as f:
        for l in f.readlines():
            ls = l.split("=")
            if ls[0].strip() == "!ss scanType":
                return int(ls[1].strip())
        return None
        


def sbatchInDir(dirname):
    """ Launches a job inside the named directory.
    dirname -- string with path to directory to scan in
    returns: a job object, containing jobid and dirname and functions to check status.
    """
    tmp = os.getcwd()
    os.chdir(dirname)
    # this breaks backwards compatibility, python3.7 and higher only!
    finished = subprocess.run(["sbatch", "job.sfincsScan"], capture_output=True, text=True)
    stdout = finished.stdout
    search = "Submitted batch job "
    i1 = stdout.find(search) + len(search)
    #print("stdout")
    #print(stdout)
    #print(stdout[i1:])
    jobid = stdout[i1:].split(None,1)[0]
    #jobid="1"
    os.chdir(tmp)    
    return Job(dirname,jobid)

def addToN(dirname,Nplus = 2):
    filename = os.getcwd() + "/" + dirname + "/job.sfincsScan"
    modified_file = []
    with open(filename) as f:
        for l in f.readlines():
            if l.split("=")[0].strip() == "#SBATCH --nodes":
                N = float(l.split("=")[1])
                l = "#SBATCH --nodes=" + str(int(N+Nplus))
                scale = 1+Nplus/N
            elif "srun -n" in l:
                ls1 = l.split("-n",1)
                ls2 = ls1[1].strip().split(None,1)
                n = int(float(ls2[0]) * scale)
                l = ls1[0] + "-n " + str(n) + " " + ls2[1] 
            else:
                l = l[:-1] # remove trailing \n for standardization
            modified_file.append(l + "\n")
    with open(filename,'w') as f:
        for l in modified_file:
            f.write(l)
    return int(N+Nplus)

def modifyScan2(dirname, N=5,lower=0.3, upper=0.3, dPhidr = None):
    filename = dirname + "/input.namelist"
    modified_file = []
    with open(filename) as f:
        for l in f.readlines():
            if l.split("=")[0].strip() == "dPhiHatdrHat":
                if dPhidr is None:
                    dPhiHatdrHat = float(l.split("=")[1].split("!")[0])
                else:
                    dPhiHatdrHat = dPhidr
                

            if l.split("=")[0].strip() == "!ss scanType":
                l = "!ss scanType = 2"
            elif l.split("=")[0].strip() == "!ss dPhiHatdrHatMin":
                l = "!ss dPhiHatdrHatMin = " + str(dPhiHatdrHat * (1-lower))
            elif l.split("=")[0].strip() == "!ss dPhiHatdrHatMax":
                l = "!ss dPhiHatdrHatMax = " + str(dPhiHatdrHat * (1+upper))
            elif l.split("=")[0].strip() == "!ss NErs":
                l = "!ss NErs = " + str(N)
            else:
                l = l[:-1] # remove trailing \n for standardization

            modified_file.append(l + "\n")
    with open(filename,'w') as f:
        for l in modified_file:
            f.write(l)


def modifyScan1(dirname, N=2,lower=0.3, upper=0.0, Ntheta=None, Nzeta = None, Nxi = None, Nx = None):
    filename = dirname + "/input.namelist"
    # fall-back to default number of runs unless specific number specified
    if Ntheta is None:
        Ntheta = N
    if Nzeta is None:
        Nzeta = N
    if Nxi is None:
        Nxi = N
    if Nx is None:
        Nx = N

    modified_file = []
    with open(filename) as f:
        for l in f.readlines():
            if l.split("=")[0].strip() == "!ss scanType":
                l = "!ss scanType = 1"
            elif l.split("=")[0].strip() == "!ss NthetaMinFactor":
                l = "!ss NthetaMinFactor = " + str(1-lower)
            elif l.split("=")[0].strip() == "!ss NthetaMaxFactor":
                l = "!ss NthetaMaxFactor = " + str(1+upper)
            elif l.split("=")[0].strip() == "!ss NzetaMinFactor":
                l = "!ss NzetaMinFactor = " + str(1-lower)
            elif l.split("=")[0].strip() == "!ss NzetaMaxFactor":
                l = "!ss NzetaMaxFactor = " + str(1+upper)
            elif l.split("=")[0].strip() == "!ss NxMinFactor":
                l = "!ss NxMinFactor = " + str(1-lower)
            elif l.split("=")[0].strip() == "!ss NxMaxFactor":
                l = "!ss NxMaxFactor = " + str(1+upper)
            elif l.split("=")[0].strip() == "!ss NxiMinFactor":
                l = "!ss NxiMinFactor = " + str(1-lower)
            elif l.split("=")[0].strip() == "!ss NxiMaxFactor":
                l = "!ss NxiMaxFactor = " + str(1+upper)
            elif l.split("=")[0].strip() == "!ss solverToleranceNumRuns":
                l = "!ss solverToleranceNumRuns = 0"
            elif l.split("=")[0].strip() == "!ss NLNumRuns":
                l = "!ss NLNumRuns = 0"
            elif l.split("=")[0].strip() == "!ss NxNumRuns":
                l = "!ss NxNumRuns = " + str(Nx)
            elif l.split("=")[0].strip() == "!ss NxiNumRuns":
                l = "!ss NxiNumRuns = " + str(Nxi)
            elif l.split("=")[0].strip() == "!ss NzetaNumRuns":
                l = "!ss NzetaNumRuns = " + str(Nzeta)
            elif l.split("=")[0].strip() == "!ss NthetaNumRuns":
                l = "!ss NthetaNumRuns = " + str(Ntheta)

            else:
                l = l[:-1] # remove trailing \n for standardization

            modified_file.append(l + "\n")
    with open(filename,'w') as f:
        for l in modified_file:
            f.write(l)


def changeVar(dirname,group,var,value):
    # Warning: this command will fail silently if the pattern is not found. Sorry about that.
    # Warning: case insensitive
    filename = dirname + "/input.namelist"

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
        # escape slashes
        value = value.replace("/","\/")
    elif (type(value) == list) or (type(value) == np.ndarray):
        # arrays are space seperated
        delimiter=' '
        value_temp = '' 
        for val in value:
            value_temp =  value_temp + str(val) + delimiter
        value = value_temp.rsplit(delimiter,1)[0]
    else:
        pass
    subprocess.call("sed -i -e '/\&"+group+"/I,/\&/{ s/^  "+var+" =.*/  "+var+" = "+str(value)+"/I } ' "+filename, shell=True)



if __name__=="__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    if argc>1:
        dirname = argv[1]
    else:
        dirname = "."
    
    scan = scanInDir(dirname)
    print(scan.status)
    job = sbatchInDir(dirname+"/baseCase")
    print(job.status)
    job2 = Job.fromDir(dirname)
    print(job2.status)
    
    addToN(dirname)
