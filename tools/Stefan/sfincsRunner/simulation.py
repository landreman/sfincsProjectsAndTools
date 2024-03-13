from __future__ import division, print_function

import subprocess
import os
from shutil import copy

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from sfincs_simulation import Sfincs_simulation

from sfincs_utils import sbatchInDir, squeue, scancel, ErrOut, getLatestJobIDInDir
# for testing purposes:
# from simulMockup import sbatchInDir, squeue, scancel


class OutputNotFound(Exception):
    def __init__(self, dirname):
        self.info = dirname


        

class Simulation(object):
    """Object that represents a run of SFINCS with a given set of physics and resolution parameters. A 'simulation' can correspond to several jobs on the cluster, if the earlier jobs fail to converge; however, numerical resolution and physics parameters are never allowed to change, as this is considered a 'different' simulation.
    The simulation object is responsble for updating itself and relaunching, in the event that a job fails.

    Usage (if the jobID of a simulation running in dirname is known):
    > Simulation(dirname,jobID) 
    Usage (if only the directory is known; will try to deduce jobID from finished or queued/running jobs)
    > Simulation.fromDir(dirname): 

    Run input variables (that a simulation may update if a run fails):
    - time: the time limit
    - N nodes: number of nodes
    - icntl parameter to MUMPS
    - KSP parameters
    - Preconditioner settings (maybe not?)

    Run output variables:
    - status (also available if the run fails)
    - SFINCS outputs (only if the run is successful)
    
    Constants after initiation (create a new simulation object if these need to change)
    - dirname: simulation directory
    - input namelist
    - job file
    - Most variables in the input namelist (maybe only physics and resolution params?)    
    """
    
    @classmethod
    def fromDir(cls,dirname):
        jobID = getLatestJobIDInDir(dirname)
        return cls(dirname,jobID)

    def __init__(self,dirname,jobID):
        self._dirname = dirname # needs to be set before jobID setter
        self.jobID = jobID # needs to be set before dirname
        self.dirname = dirname 
        self.auto = True

        
    def run(self, rerun = True, force=False, timelimit = 4*3600, Nodelimit = 16):
        """Runs the simulation.

        Inputs:
        rerun -- update time limit and number of nodes and re-run the simulation if it is in a failed state (TIME, OOM). Default: True.
        force -- Re-run even successful simulations. Cancels running simulations, if the jobID is known. Default: False.
        timelimit -- Maximum time a simulation will set. Default: 4h
        Nodelimit -- Maximum number of nodes a simulation will set. Default: 16"""
        if self.status == "QUEUED":
            if force:
                tmp = scancel(j=self.jobID)
                self.jobID = sbatchInDir(self.dirname)
        
        elif self.status == "DONE":
            if force:
                self.jobID = sbatchInDir(self.dirname)
        elif self.status == "RUNNING":
            if force:
                tmp = scancel(j=self.jobID)
                self.jobID = sbatchInDir(self.dirname)
        elif self.status == "TIME":
            # TODO maybe do KSP change if above a threshold number of KSP iter.
            if force or rerun:
                time = 2*self.time
                if time <= timelimit:
                    print("Set new time to " + str(time))
                    # here, the job will simply wait
                    # rerunning this piece of code whenever run() is called
                    self.time = time
                    self.jobID = sbatchInDir(self.dirname)

        elif self.status == "OOM":
            if force or rerun:
                N = self.N + 2
                if N <= Nodelimit:
                    print("Set new N to " + str(N))
                    # here, the job will simply wait
                    # rerunning this piece of code whenever run() is called
                    self.N = N
                    self.jobID = sbatchInDir(self.dirname)
        else:
            # unknown status
            if force:
                self.jobID = sbatchInDir(self.dirname)
        

    def getLatestJobID(self):
        if not self.inQueue():
            # if we are queued, we can trust the jobID we already have
            # otherwise, we get it from here
            self.jobID =  getLatestJobIDInDir(self.dirname)
        return self.jobID

    @property
    def statestring(self):
        if self.status == "QUEUED":
            ret = "QUEU"
        elif self.status == "OOM":
            ret = "OOM"
        elif self.status == "TIME":
            ret = "TIME"
        elif self.status == "DONE":
            ret = "DONE"
        elif self.status == "RUNNING":
            ret = "RUN"
        else:
            ret = "???"
        return ret

    @property
    def short_dirname(self):
        if self.dirname[-1] == "/":
            dirname = self.dirname[:-1]
        else:
            dirname = self.dirname
        drs = dirname.rsplit("/",2)
        return drs[1] +"/" + drs[2]

    @property
    def info(self):
        return "simulation info"

    @property
    def typestring(self):
        return "simul"


    @property
    def dirname(self):
        return self._dirname

    @dirname.setter
    def dirname(self,dn):
        self._dirname = dn
        self.jobfile = self._dirname + "/job.sfincsScan"
        if self.jobID is not None:
            self.stdout = self._dirname + "/job.out." + str(self.jobID)
            self.stderr = self._dirname + "/job.err." + str(self.jobID)
    
    
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
        if self.jobID is None:
            # retry jobID once
            self.getLatestJobID()
            if self.jobID is None:
                # no job has been run or is running.
                status = "???"

        elif not os.path.isfile(self.stdout):
            if self.inQueue():
                status = "QUEUED"
        else:
            status = ErrOut(self.dirname,self.jobID).status
        return status


    def inQueue(self):
        tmp = squeue(j=self.jobID)
        if len(tmp) > 0:
            ret = True
        else:
            ret = False
        return ret

     # TODO return default if not found in file???
     # TODO write time and N as specfic cases of a general function
    @property
    def time(self):
        with open(self.jobfile,'r') as f:
            t = f.read().rsplit("#SBATCH --time=",1)[1].strip().split(None,1)[0]
        h,m,s = t.split(":")
        h = int(h)
        m = int(m)
        s = int(s)
        return h*3600 + m*60 + s

    @time.setter
    def time(self,s):
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
    def N(self):
        with open(self.jobfile,'r') as f:
            N = int(f.read().rsplit("#SBATCH --nodes=",1)[1].strip().split(None,1)[0])
        return N

    @N.setter
    def N(self,N):
        new_filecontent = []
        with open(self.jobfile,'r') as f:
            for l in f.readlines():
                if "#SBATCH --nodes=" in l:
                    Nold = float(l.split("=")[1])
                    l = "#SBATCH --nodes=" + str(N) + "\n"
                    scale = N/Nold
                elif "srun -n" in l:
                    ls1 = l.split("-n",1)
                    ls2 = ls1[1].strip().split(None,1)
                    n = int(float(ls2[0]) * scale)
                    l = ls1[0] + "-n " + str(n) + " " + ls2[1] 
            
                new_filecontent.append(l)
        with open(self.jobfile,'w') as f:
            for l in new_filecontent:
                f.write(l)


    @property
    def gmres_r(self):
        # TODO return default if not found in file
        # TODO write all these as specfic cases of a general function?
        with open(self.jobfile,'r') as f:
            l = f.read()
        if "-ksp_gmres_restart" in l:
            r = int(l.rsplit("-ksp_gmres_restart",1)[1].strip().split(None,1)[0])
        else:
            # PETSc default
            r = 30
        return r


    @gmres_r.setter
    def gmres_r(self,r):
        new_filecontent = []
        with open(self.jobfile,'r') as f:
            for l in f.readlines():
                if "srun" in l:
                    if "-ksp_gmres_restart" in l:
                        ls = l.rsplit("-ksp_gmres_restart",1)
                        try:
                            ls1 = ls[1].strip().split(None,1)[1]
                        except IndexError:
                            ls1 = ""
                        l = ls[0] +"-ksp_gmres_restart " + str(r) + " " + ls1
                    else:
                        l = l[:-2] + " -ksp_gmres_restart " + str(r) + "\n"
                new_filecontent.append(l)
        with open(self.jobfile,'w') as f:
            for l in new_filecontent:
                f.write(l)


    @property
    def cntl1(self):
        with open(self.jobfile,'r') as f:
            l = f.read()
        if "-mat_mumps_cntl_1" in l:
            g = float(l.rsplit("-mat_mumps_cntl_1",1)[1].strip().split(None,1)[0])
        else:
            # MUMPS default for unsymmetric or general matrixes
            g =  0.01 
        return g

    @cntl1.setter
    def cntl1(self,r):
        new_filecontent = []
        with open(self.jobfile,'r') as f:
            for l in f.readlines():
                if "srun" in l:
                    if "-mat_mumps_cntl_1" in l:
                        ls = l.rsplit("-mat_mumps_cntl_1",1)
                        try:
                            ls1 = ls[1].strip().split(None,1)[1]
                        except IndexError:
                            ls1 = ""
                        l = ls[0] +"-mat_mumps_cntl_1 " + str(r) +" " + ls1
                    else:
                        l = l[:-2] + " -mat_mumps_cntl_1 " + str(r) + "\n"
                new_filecontent.append(l)
        with open(self.jobfile,'w') as f:
            for l in new_filecontent:
                f.write(l)


    @property
    def icntl23(self):
        with open(self.jobfile,'r') as f:
            l = f.read()
        if "-mat_mumps_icntl_23" in l:
            i = int(l.rsplit("-mat_mumps_icntl_23",1)[1].strip().split(None,1)[0])
        else:
            # MUMPS default, meaning that each processor gets a workspace
            # estimated during analysis
            i = 0
        return i


    @icntl23.setter
    def icntl23(self,r):
        new_filecontent = []
        with open(self.jobfile,'r') as f:
            for l in f.readlines():
                if "srun" in l:
                    if "-mat_mumps_icntl_23" in l:
                        ls = l.rsplit("-mat_mumps_icntl_23",1)
                        try:
                            ls1 = ls[1].strip().split(None,1)[1]
                        except IndexError:
                            ls1 = ""
                        l = ls[0] +" -mat_mumps_icntl_23 " + str(r) + " " + ls1
                    else:
                        l = l[:-2] + " -mat_mumps_icntl_23 " + str(r) + "\n"
                new_filecontent.append(l)
        with open(self.jobfile,'w') as f:
            for l in new_filecontent:
                f.write(l)


    # TODO: implement
    @property
    def precondX(self):
        pass

    @precondX.setter
    def precondX(self,r):
        pass

    @property
    def precondCol(self):
        pass

    @precondCol.setter
    def precondCol(self,r):
        pass

    def __str__(self):
        return self.dirname + ": " + self.status
        
