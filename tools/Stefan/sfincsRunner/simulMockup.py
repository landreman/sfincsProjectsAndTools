import os
from shutil import copy

def squeue(j=1111111, u = "sbul"):
    """Mocks the output of squeue. For testing purposes without actually queueing up jobs."""
    j = str(j)
    return "             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)\n           " + j + "       rvs  rvs_cpu     " + u + "  R    5:08:37      1 ravc2720\n"
    
def scancel(j=1111111,u = "sbul"):
    """Mocks the output of scancel. Since no job is ever actually queued up, this does basically nothing except mimic the interface of the non-mock version."""
    return ""

def sbatchInDir(dirname,j=1111111, status = "DONE"):
    """Mocks the output of running sbatchInDir() by returing a fake jobID and copying the stdout and stderr of an old simulation.
 
    Inputs:
    dirname -- directory that sbatch is run in. (It's a fake run.)
    j -- the jobID of the faked job.
    status -- The status of the faked job. Different stdout and stderr will be created to mimic this status. Valid values: 'DONE', 'OOM', 'TIME'."""
    
    # directory of this file
    thisdir = os.path.abspath(__file__).rsplit("/",1)[0]

    j = str(j)
    if status == "DONE":
        sourcedir = thisdir + "/tests/ErrOut/DONE"
        sourceID = "4520731"
    elif status == "OOM":
        sourcedir = thisdir + "/tests/ErrOut/OOM"
        sourceID = "4520737"
    elif status == "TIME":
        sourcedir = thisdir + "/tests/ErrOut/TIME"
        sourceID = "4520742"
    else:
        raise ValueError("sbatchInDir (mock) called with invalid status.")

    copy(sourcedir + "/job.out." + sourceID, dirname +"/job.out." + j)
    copy(sourcedir + "/job.err." + sourceID, dirname +"/job.err." + j)

    return j


if __name__=="__main__":
    print(squeue())
    tmpdir = "tmp__simulMockup"
    try:
        os.mkdir(tmpdir)
    except FileExistsError:
        sbatchInDir(tmpdir,j=666,status="OOM")
