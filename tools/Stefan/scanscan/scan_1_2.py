#!/usr/bin/env python3

from event_loop import Abstact_event_loop
from sfincs_utils import scanInDir, modifyScan1, modifyScan2, copyInput, changeVar, addToN, sbatchInDir
from set_resolution import set_resolution
from shutil import copy, move
import os

# states:
#   0: scan1
#   1: Er to find more ambipolar solution
#   2: final Er
#   3: Done

# transition between states: move to different dir or change things
# for state 1 or 2, can be in different subdir
# since state 0 is default, it should be in topdir

class Scan_1_2_event_loop(Abstact_event_loop):
    params = ["Ntheta","Nzeta","Nxi","Nx"]
    NErScan2Simuls = 3
    
    def to_state0(self,i):
        # this is only called when going from state 1 to state 0
        # not when initially entering state 0 at the beginning
        # TODO?: restore "max" in the scan
        self.states[i] = 0
        d = self.topdirlist[i]
        nd = self.topdirlist[i] + "/prelim"
        filelist = os.listdir(d) 
        os.mkdir(nd)
        for _f in filelist:
            f = d + "/" + _f
            nf = nd + "/" + _f
            if _f not in ["input.namelist", "job.sfincsScan"]:
                move(f,nf)
            else:
                copy(f,nf)
        # TODO?:
        # only scan previously unconverged parameters
        params_to_scan = {}
        for ii,p in enumerate(type(self).params):
            if not self.conv[i][ii]:
                params_to_scan[p] = 2
        changeVar(d,"physicsParameters", "dPhiHatdrHat",Er)
        modifyScan1(d,N=0,upper=self.upper,**params_to_scan)
        self.jobs[i] = scanInDir(d)
        


    def to_state1(self,i):
        self.states[i] = 1
        newdir = self.topdirlist[i] + "/Erscan"
        copyInput(self.topdirlist[i],newdir)
        modifyScan2(newdir,N=type(self).NErScan2Simuls,lower=0.3,upper=0.3)
        self.jobs[i] = scanInDir(newdir)

    def to_state2(self,i):
        self.states[i] = 2
        newdir =self.topdirlist[i] + "/Erscan"
        copyInput(self.topdirlist[i],newdir)
        # convparams[i] contain converged resolution parameters
        # so use those in the new dir
        for ii,p in enumerate(type(self).params):
            changeVar(newdir,"resolutionParameters", p,self.convparams[i][ii])
        modifyScan2(newdir,N=type(self).NErScan2Simuls,lower=0.3,upper=0.3)
        self.jobs[i] = scanInDir(newdir)

    def start(self):
        self.upper = 0.0
        Nr = len([f for f in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), f))])
        self.topdirlist = [str(x).zfill(2) for x in range(Nr)]
        self.states = [0] * Nr
        self.jobs = [None] * Nr
        self.conv = [[False] * len(type(self).params)] * Nr
        self.convparams = [[None] * len(type(self).params)] * Nr
        for i,d in enumerate(self.topdirlist):
            # set baseline resolution
            if False:
                set_resolution(d)
            # set parameters for basic scan1
            modifyScan1(d,N=2,lower=0.3,upper=0.0)
            # start scan, put list of jobs in self.jobs
            self.jobs[i] = scanInDir(d)

    def f(self):
        if all([x == 3 for x in self.states]):
            self.print("All done. Goodbye!")
            exit(0)

        for i,jl in enumerate(self.jobs):
            # remove completed jobs
            self.print("====== " + str(i).zfill(2) + " ======")
            self.print(jl)
            self.print(self.states[i])
            jl = [j for j in jl if j.status != "DONE"]
            for j in jl:
                if j.status == "OOM":
                    # remove old job, add a new job with
                    # 2 more compute nodes
                    newN = addToN(j.dirname,Nplus=2)
                    if newN <= 16:
                        jl.remove(j)
                        jl.append(sbatchInDir(j.dirname))
                    else:
                        #Do not remove from jobs[i]"
                        #Print problem to log
                        self.print("**ERR OOM " + j.dirname + "\n")
                        # check if the user has manually started a newer job
                        # by getting the latest (highest) jobID
                        j.retryJobID()
                elif j.status == "TIME":
                    newTime = j.getTimeLimit() * 2
                    if newTime <= 14400: # 4h in seconds
                        jl.remove(j)
                        jl.append(sbatchInDir(j.dirname))
                    else:
                        #Do not remove from jobs[i]"
                        #Print problem to log
                        self.print("**ERR TIME " + j.dirname + "\n")
                        # check if the user has manually started a newer job
                        # by getting the latest (highest) jobID
                        j.retryJobID()
            self.jobs[i].jobs = jl
            if len(jl) == 0:
                # all jobs are finished
                # what to do next depends on our state
                if self.states[i] == 0:
                    self.conv[i],convparams = self.jobs[i].check_convergence(type(self).params)
                    # only update params that are not None
                    for iii,c in enumerate(convparams):
                        if c is not None:
                            self.convparams[i][iii] = c
                    if not all(self.conv[i]):
                        # "if not converged, check ambipolar"
                        ambi = self.jobs[i].check_ambipolarity(0.05)
                        if not ambi:
                            # "if not ambpolar, do Er scan (state 1)"
                            # TODO: could update resolution parameters
                            # based on convergence
                            self.to_state1(i)
                        else:
                            # "if ambipolar, add more points to convergence scan"
                    
                            self.upper = self.upper + 0.3
                            # add extra points to scan
                            params_to_scan = {}
                            for ii,p in enumerate(type(self).params):
                                if not self.conv[i][ii]:
                                    params_to_scan[p] = 2
                            modifyScan1(self.jobs[i].dirname,N=0,upper=self.upper,**params_to_scan)
                            self.jobs[i] = scanInDir(self.jobs[i].dirname)
                    else:
                        # do final Er scan (state 2)
                        self.to_state2(i)
                else: 
                    # state = 1 or state = 2
                    # check if Er has been solved for

                    # TODO: can in principle return many roots
                    # but this is not typical with how Sfincs is run 
                    # and not supported right now. Script will fail if it happens.
                    Er = self.jobs[i].roots
                    if Er is None:
                        # did not find Er, add more points to scan
                        # based on extrapolation for Er
                        alpha = 0.1
                        newEr = self.jobs[i].extrapolate_for_root()
                        if newEr is not None:
                            if newEr > self.jobs[i].maxEr:
                                l = alpha * (1 - self.jobs[i].maxEr/newEr)
                            elif newEr < self.jobs[i].minEr:
                                l = alpha * (self.jobs[i].minEr/newEr - 1)
                            else:
                                # this should be impossible
                                raise ValueError("Impossible!")
                            modifyScan2(self.jobs[i].dirname,N=type(self).NErScan2Simuls,lower=l,upper=l,dPhidr = newEr)
                            self.jobs[i] = scanInDir(self.jobs[i].dirname)
                        else:
                            "Extrapolation failed?"
                            "write to log of problematic simulations and wait"
                            
                    elif self.states[i] == 1:
                        # perform a new resolution scan
                        self.to_state0(i)
                    else:
                        # check if Er has datapoints sufficiently close on both sides
                        # if not, make a new scan with points close to Er added
                        Ers  = self.jobs[i]._Ers
                        print("!!!!!!!!!!!!!")
                        print(Er)
                        print(Ers)
                        if any((Ers<Er) * (Ers>Er*0.9)):
                            if any((Ers>Er) * (Ers<Er*1.1)):
                                # points close enough exist. We are done!
                                self.states[i] = 3
                            else:
                                modifyScan2(self.jobs[i].dirname,N=2,lower=0,upper=l,dPhidr = Er)
                                self.jobs[i] = scanInDir(self.jobs[i].dirname)
                                # no Er close above
                                # Er close underneath
                        else:
                            if any((Ers>Er) * (Ers<Er*1.1)):
                                modifyScan2(self.jobs[i].dirname,N=2,lower=0.1,upper=0,dPhidr = Er)
                                self.jobs[i] = scanInDir(self.jobs[i].dirname)
                                
                                # no Er close underneath
                                # Er close above
                            else:
                                # No Er close above or underneath
                                print("bad under over")
                                modifyScan2(self.jobs[i].dirname,N=3,lower=0.1,upper=0.1,dPhidr = Er)
                                print(self.jobs[i].dirname)
                                self.jobs[i] = scanInDir(self.jobs[i].dirname)
                                
if __name__ == "__main__":
    f = open("log.txt","a")
    el = Scan_1_2_event_loop(f)
    el.event_loop(10.0,{})
