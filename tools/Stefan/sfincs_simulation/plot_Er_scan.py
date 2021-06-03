#!/usr/bin/env python3

from __future__ import division, print_function
import numpy as np
from Er_scan import Er_scan
import matplotlib.pyplot as plt

def extract_particle_fluxes(Er_scan):
    Gammas = []
    Ers = []
    for s in Er_scan.simuls:
        Gammas.append(s.Gamma_s)
        Ers.append(s.Er)
    Gammas = np.array(Gammas) # Gammas[simulation,species]
    Ers = np.array(Ers) # Ers[simulation]
    return Ers, Gammas

def plot_Ers_Gammas(Ers,Gammas):
    Nspecies = len(Gammas[0])
    for i in range(Nspecies):
        plt.plot(Ers,Gammas[:,i])
    
def plot_Er_scan(Er_scan):
    Ers, Gammas = extract_particle_fluxes(Er_scan)

    plot_Ers_Gammas(Ers, Gammas)
    
    plt.legend(Er_scan.simuls[0].species.names)
        
    Nroots = len(Er_scan.roots)
    if (Nroots == 3) and (Er_scan.use_roots == "i&e" or Er_scan.use_roots == "e&i"):
        Nroots = 2
    elif (Nroots == 3) and (Er_scan.use_roots == "i" or Er_scan.use_roots == "e"):
        Nroots = 2
        
    for iroot in range(Nroots):
        if Nroots>1:
            Gamma_s = Er_scan.Gamma_s[iroot]
            Er = Er_scan.Er[iroot]
        else:
            Gamma_s = Er_scan.Gamma_s
            Er = Er_scan.Er
        for i in range(Nspecies):
            plt.plot(Er,Gamma_s[i],'*')
    plt.xlabel(r"$E_r/[\rm{V}\rm{m}^{-1}]$")
    plt.ylabel(r"$\Gamma/[\rm{s}^{-1}]$")


if __name__ == "__main__":
    Es = Er_scan(".")
    plot_Er_scan(Es)
    plt.savefig("Erscan.pdf")

    print(Es.Er)
    
    
