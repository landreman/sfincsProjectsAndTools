from __future__ import division

from numpy import sum,sqrt,array,zeros,exp
from integrals import F,K
from FSA import FSA

def calculate_classical_transport(Zs,mHats,NHats,THats,ddpsiHat_NHats, ddpsiHat_THats, ddpsiHat_PhiHat,Delta,alpha,nu_n,nablaPsiHat2,BHat,Phi1Hat=None):
    """Calculates the classical transport given masses, charge and profiles of all the species, and a geometric factor |\nabla \psi|^2. Everything should be given in terms of the internal normalization in SFINCS, and the species dependent quantities should be passed as a list, tuple or numpy.array where Zs[i] and mHats[i] is charge-number and mass of species i, and so on."""

    # get number of species
    Nspecies = len(mHats)
    # check that the number of species is consistent
    assert(all([len(x) == Nspecies for x in [Zs,NHats,THats,ddpsiHat_NHats, ddpsiHat_THats]]))

    nHats = NHats[:]
    #include Phi1Hat in nHats if appropriate
    if Phi1Hat is None:
        includePhi1 = False
    else:
        includePhi1 = True
        for i in range(Nspecies):
            nHats[i] = nHats[i] * exp(-Zs[i]*alpha*Phi1Hat/THats[i])
    
    #print (Zs,mHats,NHats,THats,ddpsiHat_NHats, ddpsiHat_THats, ddpsiHat_PhiHat,Delta,alpha,nu_n,nablaPsiHat2,BHat,Phi1Hat)
    
    classical_fluxes = zeros(Nspecies)

    for a in range(0,Nspecies):
        for b in range(0,Nspecies):
            y = THats[a] * mHats[b]/(THats[b] * mHats[a])

            classical_fluxes[a] += Zs[b]*FSA(nablaPsiHat2 * nHats[a] * nHats[b]/BHat**2,BHat) *( \
                                                                                               + Zs[a] * sqrt(mHats[b]/THats[b]) * (1 + mHats[b]/mHats[a]) * F(y) * (ddpsiHat_NHats[b]/NHats[b] + Zs[b] * alpha * ddpsiHat_PhiHat/THats[b]) \
                                                                                               + Zs[a] * sqrt(mHats[b]/THats[b]) * (1 + mHats[b]/mHats[a]) * K(y) * ddpsiHat_THats[b]/THats[b] \
                                                                                               - Zs[b] * sqrt(mHats[a]/THats[a]) * (1 + mHats[a]/mHats[b]) * F(1/y) * (ddpsiHat_NHats[a]/NHats[a] + Zs[a] * alpha * ddpsiHat_PhiHat/THats[a]) \
                                                                                                 - Zs[b] * sqrt(mHats[a]/THats[a]) * (1 + mHats[a]/mHats[b]) * K(1/y) * ddpsiHat_THats[a]/THats[a] \
            )
            if includePhi1:
                classical_fluxes[a] += alpha*Zs[a] *Zs[b]**2 * FSA(nablaPsiHat2 * nHats[a] * nHats[b] *Phi1Hat/BHat**2,BHat) * (\
                                                                                                            sqrt(mHats[b]/THats[b]) * (1 + mHats[b]/mHats[a]) * F(y) * ddpsiHat_THats[b]/(THats[b]**2) \
                                                                                                            - sqrt(mHats[a]/THats[a]) * (1 + mHats[a]/mHats[b]) * F(1/y) * ddpsiHat_THats[a]/(THats[a]**2) \
                )
            
    classical_fluxes = 2*Delta**2 * nu_n * classical_fluxes

    return classical_fluxes

if __name__ == "__main__":
    # debug against mass ratio expanded expression previously calculated
    import sys, os
    import numpy
    print __file__
    print os.path.dirname(__file__)

    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Hakan/pythonBoozerFilesAndGeom'))
    from bcgeom import bcgeom
    from fluxcoorddiscr import fluxcoorddiscr

    # can artificially increase mass and Z to try to match Z>>1, mC/mD >> 1 asymptotic results
    mass_increase = 1
    Z_increase = 1
    
    e =  1.6021766208e-19
    epsilon0 = 8.854187817620389e-12
    mD = 3.349854702e-27
    me = 9.10938356e-31
    mC = mass_increase * 1.9926467051999998e-26
    ZD = 1
    Ze = -1
    ZC = Z_increase * 6
    
    pi = numpy.pi

    geom = bcgeom('w7x-sc1.bc',1e-16,float("inf"),float("inf"),'StelSym',1)

    rind=34
    Ntheta = 91
    Nzeta = 81
    Booz = fluxcoorddiscr(geom,rind,Ntheta,Nzeta,name='Boozer')
    BHat = Booz.B
    nablaPsiHat2 = Booz.gpsipsi

    mBar = mD 
    nBar = 1e20 #reference number density in m^-3
    TBar = 1000*e # reference T of 1keV in J
    BBar = 1 # 1 Tesla
    RBar = 1 # 1 meter
    alpha = 1 # ePhiBar/TBar
    Delta = sqrt(2*mBar*TBar)/(e*BBar*RBar)

    # Inventend parameters for D-C-e plasma
    lnLambda = 17
    nu_n = (nBar * RBar * e**4 * lnLambda/(12 * pi**(3/2) * epsilon0**2 * TBar**2))
    mHats = array([mD/mBar, mC/mBar, me/mBar])
    Zs = array([ZD,ZC,Ze])

    nHats = array([1, 0.001, ZD*1 + ZC * 0.001])
    THats = array([1, 1, 1])
    ddpsiHat_nHats = -2 * nHats
    ddpsiHat_THats = -2 * THats
    ddpsiHat_PhiHat = 0

    normalized_classical_transport = calculate_classical_transport(Zs,mHats,nHats,THats,ddpsiHat_nHats, ddpsiHat_THats, ddpsiHat_PhiHat,Delta,alpha,nu_n, nablaPsiHat2,BHat)
    classical_transport = (nBar*sqrt(2*TBar/mBar)/RBar) * normalized_classical_transport/Z_increase
    print classical_transport
    quasi_neutrality = sum(Zs*normalized_classical_transport)
    print quasi_neutrality
