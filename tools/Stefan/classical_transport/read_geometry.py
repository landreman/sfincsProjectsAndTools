import numpy as np
import sys,os

from calculate_classical_transport import calculate_classical_transport 
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Hakan/pythonBoozerFilesAndGeom'))
from bcgeom import bcgeom
from fluxcoorddiscr import fluxcoorddiscr

def read_psiA(geom,signcorr=2):
    if signcorr == 1:
        psiA = geom.psi_a
    elif signcorr == 2:
        psiA = geom.psi_a
    else:
        raise ValueError("signcorr must be 1 or 2.")
    return psiA

def read_geometry(inputRadialCoordinate,rN_wish,Ntheta,Nzeta,equilibriumFile,min_Bmn_to_load,max_m=float("inf"),maxabs_n = float("inf"),symmetry="unknown",signcorr=2,BBar=1,RBar=1):

    # use the right coordinate:
    # in bcgeom file, rnorm is sqrt(normalized_toroidal_flux)
    # and s is normalized_toroidal_flux

    geom = bcgeom(equilibriumFile,min_Bmn_to_load,max_m,maxabs_n ,symmetry,signcorr)
    psiA = read_psiA(geom,signcorr)
    psiAHat = psiA/(BBar*RBar**2)

    a = geom.minorradiusVMEC
    aHat = a/RBar

    # find which flux surface the use wants to look at
    if inputRadialCoordinate == 0:
        # wish coordinate is psiHat
        psiN_wish = rN_wish/psiAHat
    elif inputRadialCoordinate == 1:
        # wish coordinate is psiN
        psiN_wish = rN_wish
    elif inputRadialCoordinate == 2:
        # wish coordinate is rHat
        psiN_wish = rN_wish**2/aHat**2
    elif inputRadialCoordinate == 3:
        # wish coordinate is rN
        psiN_wish = rN_wish**2

    rind = np.argmin(np.abs(geom.s-psiN_wish))
    psiN = geom.s[rind]


    #print psiN
    #print rind

    Booz = fluxcoorddiscr(geom,rind,Ntheta,Nzeta,name='Boozer')
    Booz.s = psiN
    Booz.r = np.sqrt(psiN)
    return (Booz,geom)
    
if __name__=="__main__":
    argv = sys.argv
    if len(argv) == 1:
        equilibriumFile = '/home/bstefan/dwnloads/sfincs/fortran/version3/examples/filteredW7XNetCDF_2species_magneticDrifts_withEr/../../../../equilibria/w7x_standardConfig.bc'
        signcorr = 1
    elif len(argv) == 2:
        equilibriumFile  = sys.argv[1]
        signcorr = 1
    elif len(argv) == 3:
        equilibriumFile  = sys.argv[1]
        signcorr = int(sys.argv[2])
    else:
        raise Warning("Too many input arguments to " + argv[0])
        namelist_filename = sys.argv[1]
        signcorr = int(sys.argv[2])

    
    min_Bmn_to_load = 0.001
    max_m=float("inf")
    maxabs_n = float("inf")
    symmetry="StelSym"

    BBar = 1
    RBar = 1
    Ntheta = 51
    Nzeta = 29
    inputRadialCoordinate = 3
    rN_wish = 0.6
    
    
    Booz,geom=read_geometry(inputRadialCoordinate,rN_wish,Ntheta,Nzeta,equilibriumFile,min_Bmn_to_load,max_m,maxabs_n ,symmetry,signcorr,BBar,RBar)
    print Booz.disp()
    print geom.disp()
