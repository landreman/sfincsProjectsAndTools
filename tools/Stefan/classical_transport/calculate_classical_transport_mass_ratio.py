from __future__ import division

from numpy import sum,sqrt,array,pi,zeros
from integrals import F,K
from FSA import FSA

def calculate_classical_transport(Zs,mHats,nHats,THats,ddpsiHat_nHats, ddpsiHat_THats, ddpsiHat_PhiHat,Delta,alpha,nu_n, nablaPsiHat2,BHat):
    """Calculates the classical transport given masses, charge and profiles of all the species, and a geometric factor |\nabla \psi|^2. Everything should be given in terms of the internal normalization in SFINCS, and the species dependent quantities should be passed as a list, tuple or numpy.array where Zs[i] and mHats[i] is charge-number and mass of species i, and so on."""

    # get number of species
    Nspecies = len(mHats)
    # check that the number of species is consistent
    assert(all([len(x) == Nspecies for x in [Zs,nHats,THats,ddpsiHat_nHats, ddpsiHat_THats]]))

    classical_fluxes = zeros(Nspecies)

    a = 1
    b = 0
    #y = THats[a] * mHats[b]/(THats[b] * mHats[a])
            
    classical_fluxes[a] += Zs[b]*FSA(nablaPsiHat2 * nHats[a] * nHats[b]/BHat**2,BHat) *( \
                                                                                         + Zs[a]/4* sqrt(mHats[b]/THats[b])   * ddpsiHat_nHats[b]/nHats[b]  \
                                                                                         - Zs[a]/8 * sqrt(mHats[b]/THats[b])  * ddpsiHat_THats[b]/THats[b] \

    )

    classical_fluxes = 2*Delta**2 * nu_n * classical_fluxes

    return classical_fluxes

if __name__ == "__main__":
    # debug against mass ratio expanded expression previously calculated
    import sys, os
    print __file__
    print os.path.dirname(__file__)

    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Hakan/pythonBoozerFilesAndGeom'))
    from bcgeom import bcgeom
    from fluxcoorddiscr import fluxcoorddiscr

    e =  1.6021766208e-19
    epsilon0 = 8.854187817620389e-12
    mD = 3.349854702e-27
    me = 9.10938356e-31
    mC = 1.9926467051999998e-26

    geom = bcgeom('w7x-sc1.bc',1e-16,float("inf"),float("inf"),'StelSym',1)

    rind=34
    Ntheta = 91
    Nzeta = 81
    Booz = fluxcoorddiscr(geom,rind,Ntheta,Nzeta,name='Boozer')
    BHat = Booz.B
    nablaPsiHat2 = Booz.gpsipsi

    mBar = mD 
    nBar = 1e20 #reference n in m^-3
    TBar = 1000*e # reference T of 1keV in J
    BBar = 1 # 1 Tesla
    RBar = 1 # 1 meter
    alpha = 1 # ePhiBar/TBar
    Delta = sqrt(2*mBar*TBar)/(e*BBar*RBar)

    # Inventend parameters for D-C-e plasma
    lnLambda = 17
    nu_n = (nBar * RBar * e**4 * lnLambda/(12 * pi**(3/2) * epsilon0**2 * TBar**2)) 
    mHats = array([mD/mBar, mC/mBar, me/mBar])
    Zs = array([1,6,-1])

    nHats = array([1, 0.001, 1.006])
    THats = array([1.0, 1.1, 1.2])
    ddpsiHat_nHats = -2 * nHats
    ddpsiHat_THats = -2 * THats
    ddpsiHat_PhiHat = 0

    print e
    print THats*TBar
    print nHats*nBar
    print lnLambda
    print mBar*mHats
    
    normalized_classical_transport = calculate_classical_transport(Zs,mHats,nHats,THats,ddpsiHat_nHats, ddpsiHat_THats, ddpsiHat_PhiHat,Delta,alpha,nu_n, nablaPsiHat2,BHat)
    classical_transport = (nBar*sqrt(2*TBar/mBar)/RBar) * normalized_classical_transport
    print classical_transport
    quasi_neutrality = sum(Zs*classical_transport)
    print quasi_neutrality
