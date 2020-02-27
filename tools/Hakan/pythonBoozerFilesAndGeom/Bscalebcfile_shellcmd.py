#!/usr/bin/env python
from __future__ import division
import sys, copy, os, argparse, inspect, math, subprocess
import numpy as np
import datetime

import mnFourierlib
from fluxcoorddiscr import fluxcoorddiscr
from geomlib import bcgeom
#import netCDF4
from netCDF4 import Dataset
#import timeit
import scipy.integrate as integrate

#
# Note that some programs cannot read the output from this,
# because m0b and n0b do not reflect the content anymore
#

###############################################################
def parse_args(argv):
    parser = argparse.ArgumentParser(description='bcfile scaling to Baxis_zeta0')
    parser.add_argument('infile', type=str)
    parser.add_argument('outfile', type=str)
    parser.add_argument('Baxis_zeta0', type=float)
    parser.add_argument('min_Bmn', type=float)
    return parser.parse_known_args(argv[1:])

###############################################################
def Bscalebcfile(infile,outfile,Baxis_zeta0,min_Bmn=0,max_m=np.inf,maxabs_n=np.inf):
  Geom=bcgeom(infile,min_Bmn=min_Bmn,max_m=max_m,maxabs_n=maxabs_n,symmetry='StelSym')
  Booz=fluxcoorddiscr(Geom,0,15,101,name='Boozer')
  Baxis_zeta0_in=np.mean(Booz.B[:,0])
  print('Booz.B='+str(Booz.B))
  print('Booz.B[:,0]='+str(Booz.B[:,0]))
  print('Geom.s[0]='+str(Geom.s[0]))
  scalefactor=Baxis_zeta0/Baxis_zeta0_in
  for rind in range(len(Geom.B)):
    Geom.B[rind]*=scalefactor
  Geom.B00       *=scalefactor
  Geom.FSAB2     *=scalefactor**2
  Geom.Btheta    *=scalefactor
  Geom.Bphi      *=scalefactor
  Geom.torfluxtot*=scalefactor
  Geom.psi_a     *=scalefactor
  Geom.headertext.maincomment=('CC The magnetic field in this file has been scaled '+
                               'by a factor {:.3f}'.format(scalefactor)+' from the\n'+
                               'CC original to get B={:.2f}'.format(Baxis_zeta0)+' T '+
                               'on axis at phi=0.\n'+Geom.headertext.maincomment)
  Geom.write(outfile)
  
###############################################################

arg, arg_remains = parse_args(sys.argv)
Bscalebcfile(arg.infile,arg.outfile,arg.Baxis_zeta0,arg.min_Bmn)
