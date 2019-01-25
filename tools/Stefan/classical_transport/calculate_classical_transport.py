from __future__ import division

from numpy import sum,sqrt,array,zeros,exp
from integrals import Fa as F, F2a as F2
from FSA import FSA

def calculate_classical_transport(Zs,mHats,NHats,THats,ddpsiHat_NHats, ddpsiHat_THats, Delta,alpha,nu_n,nablaPsiHat2,BHat,Phi1Hat=None,verbose=0):
    """Calculates the classical transport given masses, charge and profiles of all the species, and a geometric factor |\nabla \psi|^2. Everything should be given in terms of the internal normalization in SFINCS, and the species dependent quantities should be passed as a list, tuple or numpy.array where Zs[i] and mHats[i] is charge-number and mass of species i, and so on."""

    # get number of species
    Nspecies = len(mHats)
    # check that the number of species is consistent
    assert(all([len(x) == Nspecies for x in [Zs,NHats,THats,ddpsiHat_NHats, ddpsiHat_THats]]))

    #Phi1Hat=None
    nHats = NHats[:]
    #include Phi1Hat in nHats if appropriate
    if Phi1Hat is None:
        includePhi1 = False
        Phi1Hat = 0.0
    else:
        includePhi1 = True
        for i in range(Nspecies):
            nHats[i] = nHats[i] * exp(-Zs[i]*alpha*Phi1Hat/THats[i])

    #print (Zs,mHats,NHats,THats,ddpsiHat_NHats, ddpsiHat_THats, ddpsiHat_PhiHat,Delta,alpha,nu_n,nablaPsiHat2,BHat,Phi1Hat)

    classical_fluxes_OLD = zeros(Nspecies)
    classicalPF = zeros(Nspecies)
    classicalHF = zeros(Nspecies)
    
    for a in range(0,Nspecies):
        for b in range(0,Nspecies):

            xab2 = mHats[a]*THats[b]/(mHats[b]*THats[a])

            # Braginskii matrix elements
            Mab00 = - (1+mHats[a]/mHats[b]) * (1 + xab2)
            Mab01 = - 1.5 * (1+mHats[a]/mHats[b])
            Mab11 = -(13. + 16. * xab2 + 30. * xab2**2)/4.
            Nab11 = 27. * mHats[a]/(4.*mHats[b])
            Mab00 = Mab00/((1+xab2)**(2.5))
            Mab01 = Mab01/((1+xab2)**(2.5))
            Mab11 = Mab11/((1+xab2)**(2.5))
            Nab11 = Nab11/((1+xab2)**(2.5))

            #print Mab00,Mab01,Mab11,Nab11
            
            G1ab =  FSA(nablaPsiHat2 * nHats[a] * nHats[b]/BHat**2,BHat)
            G2ab =  FSA(nablaPsiHat2 * nHats[a] * nHats[b] *Phi1Hat/BHat**2,BHat)

            #print G1ab,G2ab

            classicalPF[a] = classicalPF[a] \
                             + Zs[b]**2  *( \
                                            G1ab * Mab00 * (THats[a] * ddpsiHat_NHats[a]/(NHats[a] * Zs[a]) - THats[b] * ddpsiHat_NHats[b]/(NHats[b] * Zs[b])) \
                                            + G2ab * alpha * Mab00 * ( ddpsiHat_THats[a]/THats[a] -  ddpsiHat_THats[b]/THats[b]) \
                                            + G1ab * ((Mab00 - Mab01) *  ddpsiHat_THats[a]/Zs[a] - (Mab00 - xab2*Mab01) *  ddpsiHat_THats[b]/Zs[b]) \
                                        )

            classicalHF[a] = classicalHF[a] \
                             + Zs[b]**2 * ( \
                                                      G1ab * Mab01 * (THats[a] * ddpsiHat_NHats[a]/(NHats[a] * Zs[a]) - THats[b] * ddpsiHat_NHats[b]/(NHats[b] * Zs[b])) \
                                                      + G2ab * alpha * Mab01 * ( ddpsiHat_THats[a]/THats[a] -  ddpsiHat_THats[b]/THats[b]) \
                                                      + G1ab * ((Mab01 - Mab11) *  ddpsiHat_THats[a]/Zs[a] - (Mab01 + Nab11) *  ddpsiHat_THats[b]/Zs[b]) \
                                                  )
            # print classicalPF[a],classicalHF[a]

            classical_fluxes_OLD[a] += Zs[b]*FSA(nablaPsiHat2 * nHats[a] * nHats[b]/BHat**2,BHat) *sqrt(mHats[b]/THats[b]) * (1 + mHats[b]/mHats[a]) *F(THats[a] * mHats[b]/(THats[b] * mHats[a])) *( \
                                                                                               +   (Zs[a] * (ddpsiHat_NHats[b]/NHats[b] - ddpsiHat_THats[b]/THats[b]/2)- Zs[b] * (THats[a]/THats[b]) *  (ddpsiHat_NHats[a]/NHats[a] - ddpsiHat_THats[a]/THats[a]/2)) \
                                                                                               + 1.5* (THats[a] * mHats[b])/(THats[b] * mHats[a]+THats[a] * mHats[b])*(Zs[a] * ddpsiHat_THats[b]/THats[b] - (mHats[a]/mHats[b])  * Zs[b] * ddpsiHat_THats[a]/THats[a]) \
            )
                
                
            if includePhi1:
                classical_fluxes_OLD[a] += alpha*Zs[a] *Zs[b]**2 * FSA(nablaPsiHat2 * nHats[a] * nHats[b] *Phi1Hat/BHat**2,BHat) * sqrt(mHats[b]/THats[b]) * (1 + mHats[b]/mHats[a]) * F(THats[a] * mHats[b]/(THats[b] * mHats[a]))/THats[b] * (ddpsiHat_THats[b]/THats[b] -  ddpsiHat_THats[a]/THats[a])


        classicalPF[a] = Zs[a] * Delta**2 * nu_n * sqrt(mHats[a]) * classicalPF[a]/(2*THats[a]**1.5)
        classicalHF[a] = -Zs[a] * Delta**2 * nu_n * sqrt(mHats[a]) * classicalHF[a]/(2*sqrt(THats[a]))
        classicalHF[a] = classicalHF[a] + 2.5 * THats[a] * classicalPF[a]
    classical_fluxes_OLD = 2*Delta**2 * nu_n * classical_fluxes_OLD
    
    return classicalPF,classicalHF

if __name__ == "__main__":
    # debug against mass ratio expanded expression previously calculated
    import sys, os
    import numpy
    #sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Hakan/pythonBoozerFilesAndGeom'))
    from geomlib import bcgeom
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

    geom = bcgeom('w7x-sc1.bc',1e-16,float("inf"),float("inf"),'StelSym',1,verbose=0)

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

    normalized_classical_transport = calculate_classical_transport(Zs,mHats,nHats,THats,ddpsiHat_nHats, ddpsiHat_THats ,Delta,alpha,nu_n, nablaPsiHat2,BHat)
    classical_transport = (nBar*sqrt(2*TBar/mBar)/RBar) * normalized_classical_transport/Z_increase
    print classical_transport
    quasi_neutrality = sum(Zs*normalized_classical_transport)
    print quasi_neutrality
