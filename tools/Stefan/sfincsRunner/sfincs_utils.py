from __future__ import division, print_function

import subprocess
import os
from shutil import copy, move

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from sfincs_simulation import Sfincs_simulation

def getLatestJobIDInDir(dirname):
    """Get the ID of the latest job that has been run in a directory.
    returns None if no jobs have been run.
    Warning: will return the ID of the latest job that was run, and not jobs that are merely queued. Use the instance method 'self.getLatestJobID()' if possible."""
    tmp = subprocess.run("ls " + dirname + " | grep .out.", shell=True, capture_output=True, text=True) # list of .out files.
    filenames = tmp.stdout.split()
    if len(filenames) > 0:
        ids = [int(fn.rsplit(".",1)[-1]) for fn in filenames]
        ids.sort()
        jobID = str(ids[-1])
    else:
        jobID = None
    return jobID


class ErrOut(object):
    """An object to parse stdout and stderr produced by a Sfincs run. The main purpose is to deduce the "status" of that run, which is used within the simulation object to set its own status.

    Usage:
    ErrOut(jobID).status()"""


    donestring = "Done with the main solve."
    oomstring = "oom-kill"
    maybeoomstring = "killed"
    timestring = "time limit"


    def __init__(self,dirname,jobID=None):
        self.dirname = dirname
        if jobID is None:
            #get latest
            self.jobID = getLatestJobIDInDir(dirname)
        else:
            self.jobID = jobID
        
        if self.jobID is not None:
            if not os.path.isfile(self.stdout):
                raise ValueError("Output does not exist")
    
    @property
    def jobID(self):
        return self._jobID

    @jobID.setter
    def jobID(self,jobID):
        self._jobID = jobID
        if jobID is not None:
            self.stdout = self.dirname + "/job.out." + str(self.jobID)
            self.stderr = self.dirname + "/job.err." + str(self.jobID)
        else:
            self.stdout = None
            self.stderr = None


    @property
    def status(self):
        # TODO check for KSP nonconvergence
        # could be from 9999 or some threshold (which would look like a TIME)
        if self.stdout is not None:
            with open(self.stdout,'r') as f:
                out = f.read()
            with open(self.stderr,'r') as f:
                err = f.read()
            if ErrOut.donestring in out:
                status = "DONE"
            elif (len(err) == 0) and (len(out) > 0):
                status= "RUNNING"
            elif ErrOut.oomstring in err.lower():
                status = "OOM"
            elif ErrOut.timestring in err.lower():
                status = "TIME"
            elif ErrOut.maybeoomstring in err.lower():
                # sometimes OOM results in a different error message
                status = "OOM"
            else:
                status = None
        else:
            status = None
        return status



def squeue(u="sbul",j=None):
    if j is not None:
        tmp = subprocess.run(["squeue","-j",j], capture_output=True, text=True)
    else:
        tmp = subprocess.run(["squeue","-u",u], capture_output=True, text=True)
    return tmp.stdout

def scancel(u="sbul",j=None):
    if j is not None:
        tmp = subprocess.run(["scancel","-j",j], capture_output=True, text=True)
    else:
        tmp = subprocess.run(["scancel","-u",u], capture_output=True, text=True)
    return tmp.stdout
    

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
    os.chdir(tmp)    
    return jobid


def copyInput(dirname,newdirname):
    if not os.path.exists(newdirname):
        os.makedirs(newdirname)
        copy(dirname + "/input.namelist", newdirname)
        copy(dirname + "/job.sfincsScan", newdirname)




def scanInDir(dirname):
    """ Launches a scan inside the named directory.
    dirname -- string with path to directory to scan in
    returns: a tuple with the scan type, a list of directories in which simulations were started, and the jobIDs of the jobs. Will also include old jobs that were not finished."""
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
    olddirs = [d for d in subdirs if d not in (dirnames + [os.path.join(dirname, f) for f in ['OLD_scan1', 'scan1', 'OLD_scan2']])]
    #print("These old simulations were automatically added to the list of simulations:")
    #print(olddirs)
    # add jobs that were previously run but ran out of time
    # or OOM
    # in the former case, nothing will actually be done
    # except that the list of unfinished simulations will be non-empty
    # so that the program waits for user intervention
    for d in olddirs:
        s = ErrOut(d).status
        if s == "OOM" or s == "TIME":
            dirnames.append(x)
            jobids.append(y)

    return stype, dirnames, jobids

def scanType(dirname):
    # check which type of scan we launched
    filename = dirname + "/input.namelist"
    with open(filename) as f:
        for l in f.readlines():
            ls = l.split("=")
            if ls[0].strip() == "!ss scanType":
                return int(ls[1].strip())
        return None
        

def modifyScan2(dirname, N=5,lower=0.3, upper=0.3, dPhidr = None, absolute=False):
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
                if not absolute:
                    l = "!ss dPhiHatdrHatMin = " + str(dPhiHatdrHat * (1-lower))
                else:
                    l = "!ss dPhiHatdrHatMin = " + str(lower)
            elif l.split("=")[0].strip() == "!ss dPhiHatdrHatMax":
                if not absolute:
                    l = "!ss dPhiHatdrHatMax = " + str(dPhiHatdrHat * (1+upper))
                else:
                    l = "!ss dPhiHatdrHatMax = " + str(upper)
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

def moveScan(olddir,newdir):
    filelist = os.listdir(olddir)
    os.mkdir(newdir)
    for _f in filelist:
        f = olddir + "/" + _f
        nf = newdir + "/" + _f
        if os.path.isdir(f) and (_f not in ["OLD_scan1","OLD_scan2","scan1"]):
            move(f,nf)
        elif _f in ["input.namelist", "job.sfincsScan"]:
            copy(f,nf)
        

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
