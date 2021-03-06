import numpy as np
import sys,os

from calculate_classical_transport import calculate_classical_transport 
#sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Hakan/pythonBoozerFilesAndGeom'))
from geomlib import bcgeom, vmecgeom
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

    if equilibriumFile.split('.')[-1] == "nc":
        raise Warning("vmec files are not supported!")
    else:
        geom = bcgeom(equilibriumFile,min_Bmn_to_load,max_m,maxabs_n ,symmetry,signcorr,verbose=0)
        
    psiA = read_psiA(geom,signcorr)
    psiAHat = psiA/(BBar*RBar**2)

    a = geom.minorradiusVMEC
    if a is np.nan:
        a = geom.minorradiusW7AS

    if a is np.nan:
        raise ValueError("minor radius from geom object is NaN!")
    
    aHat = a/RBar

    # find which flux surface the use wants to look at
    if inputRadialCoordinate == 0:
        # wish coordinate is psiHat
        rN_wish = np.sqrt(rN_wish/psiAHat)
    elif inputRadialCoordinate == 1:
        # wish coordinate is psiN
        rN_wish = np.sqrt(rN_wish)
    elif inputRadialCoordinate == 2:
        # wish coordinate is rHat
        rN_wish = rN_wish/aHat
    elif inputRadialCoordinate == 3:
        # wish coordinate is rN
        rN_wish = rN_wish
    elif inputRadialCoordinate == 4:
        # wish coordinate is rHat
        rN_wish = rN_wish/aHat
    else:
        raise ValueError("inputRadialCoordinate should be 0,1,2,3,4; it is" + str(inputRadialCoordinate))


    rind = np.argmin(np.abs(np.sqrt(geom.s)-rN_wish))
    psiN = geom.s[rind]
    
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
