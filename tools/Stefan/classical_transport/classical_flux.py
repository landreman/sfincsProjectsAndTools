#!/usr/bin/env python
import sys
from read_namelist import read_namelist
from read_geometry import read_geometry, read_psiA
from calculate_classical_transport import calculate_classical_transport
from numpy import sqrt, nan

def classical_flux(namelist_filename,max_m,maxabs_n,symmetry,signcorr,BBar=1,RBar=1):

    (Zs,mHats,nHats,THats,dNHatdrHats, dTHatdrHats, dPhiHatdrHat,Delta,alpha,nu_n,inputRadialCoordinate,inputRadialCoordinateForGradients,rN_wish,equilibriumFile,min_Bmn_to_load,Ntheta,Nzeta,Phi1Hat) = read_namelist(namelist_filename)
    
    Booz,geom=read_geometry(inputRadialCoordinate,rN_wish,Ntheta,Nzeta,equilibriumFile,min_Bmn_to_load,max_m,maxabs_n ,symmetry,signcorr,BBar,RBar)

    BHat = Booz.B/BBar
    nablaPsiHat2 = Booz.gpsipsi/((BBar*RBar)**2)
    # 2018-04-26:
    # confirmed that this is the same value as what is read into SFINCS, at least for
    # RBar = 1 m and BBar = 1 T

    
    #convert to proper gradients
    psiAHat = read_psiA(geom,signcorr)/(BBar*RBar**2)
    psiN = Booz.s
    psiHat = Booz.s * psiAHat
    
    
    a = geom.minorradiusVMEC
    if a is nan:
        a = geom.minorradiusW7AS

    if a is nan:
        raise ValueError("minor radius from geom object is NaN!")
        
    aHat = a/RBar

    if inputRadialCoordinateForGradients == 0:
        conversion_factor = 1.0
    elif inputRadialCoordinateForGradients == 1:
        conversion_factor = 1/psiAHat
    elif inputRadialCoordinateForGradients == 2 or inputRadialCoordinateForGradients == 4:
        conversion_factor = aHat/(2*psiAHat*sqrt(psiN))
    elif inputRadialCoordinateForGradients == 3:
        conversion_factor = 1/(2*psiAHat*sqrt(psiN))
    else:
        raise ValueError("inputRadialCoordinateForGradients must be 0,1,2,3 or 4")

    dNHatdrHats =[conversion_factor*e for e in dNHatdrHats]
    dTHatdrHats =[conversion_factor*e for e in dTHatdrHats]
    dPhiHatdrHat = conversion_factor*dPhiHatdrHat

    return calculate_classical_transport(Zs,mHats,nHats,THats,dNHatdrHats, dTHatdrHats,Delta,alpha,nu_n, nablaPsiHat2,BHat,Phi1Hat)

if __name__=="__main__":
    argv = sys.argv
    if len(argv) == 2:
        namelist_filename = sys.argv[1]
        signcorr = 1
    elif len(argv) == 3:
        namelist_filename = sys.argv[1]
        signcorr = int(sys.argv[2])
    else:
        raise Warning("Too many input arguments to " + argv[0])
        namelist_filename = sys.argv[1]
        signcorr = int(sys.argv[2])

    max_m = float("inf")
    maxabs_n = float("inf")
    symmetry = "StelSym"
    classical_particle_flux =  classical_flux(namelist_filename,max_m,maxabs_n,symmetry,signcorr)
    print str(classical_particle_flux)[1:-1]
