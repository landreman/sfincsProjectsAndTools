#!/usr/bin/env python

import matplotlib
#import matplotlib.pyplot as plt
import h5py
import numpy
import os, sys, inspect
import warnings
import matplotlib.ticker as ticker

show_rN = True
#show_rN = False

#makePDF = True
makePDF = False
for arg in sys.argv:
   if arg.lower()=='pdf':
      makePDF = True

if makePDF:
   matplotlib.use('PDF')
else:
   #import matplotlib.pyplot as plt
   matplotlib.use('qt5agg')

import matplotlib.pyplot as plt

print ("This is "+ inspect.getfile(inspect.currentframe()))

#########
##INPUT##
#########

filename = 'sfincsOutput.h5'

FigSize = (42,24)

font = {'size':40}
matplotlib.rc('font', **font)
matplotlib.rc('lines',markeredgewidth=0,markersize=3,linewidth=2.5)
matplotlib.rc('axes',linewidth=1.5)

matplotlib.rcParams['mathtext.default'] = 'it'
matplotlib.rcParams['text.usetex'] = True

zFactor = 1.0
withTickLabels = False
##W7-X##
#xAxisTicks = [r'$0$', r'$\pi/10$', r'$2\pi/10$', r'$3\pi/10$', r'$4\pi/10$']
##LHD
xAxisTicks = [r'$0$', r'$\pi/20$', r'$2\pi/20$', r'$3\pi/20$', r'$4\pi/20$']

yAxisTicks = [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']

fig = plt.figure(figsize=FigSize)
fig.patch.set_facecolor('white')
numRows = 2
numCols = 3
#iteration = 0
numContours = 100
ContourLevels = [-6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
#ContourLevels = [] ##If empty list the contours are calculated automatically
numLevels = 5

species = 3 - 1 ##Species index starts from 0

includePhi1InCalculation = True

zMin = None # None for automatic
zMax = None # None for automatic

##OLD dIotadpsiN = 0.8777 ##Manual input

##PUT FACTOR TO 1 TO INCLUDE THE TERM AND TO 0 TO EXCLUDE THE TERM##
##DEFAULT IS ALL 1##
TermPhi1Factor = 1.0
TermdPhi1dalphaFactor = 1.0
TermPhi1ShearFactor = 1.0
TermPhi1Exponential = 1.0

TermdBdalphaFactor = 1.0
TermdBdpsifactor = 1.0
TermBShearFactor = 1.0

phiExponent = 50 ##Exponent defining peakedness of the perturbed potential, should be integer > 0.
phiCenterZeta = 0.0 #0.0
phiCenterTheta = 0.0 #numpy.pi/4.0

#phiMaxShow = True

InfoBox = True
InfoBoxXcoord = 0.50
InfoBoxYcoord = -3.0
InfoBoxFormat = '{:1.2f}'
InfoBoxLabelSize = 37
InfoBoxLineSpacing = 2

####################

#############
##END INPUT##
#############


def fmt_cbar(x, pos):
   if x == 0.0:
      return r'${}$'.format(x)
   a, b = '{:.1e}'.format(x).split('e')
   b = int(b)
   return r'${} \cdot 10^{{{}}}$'.format(a, b)

def fmt_xy_axis(x, pos):
   return r'${}$'.format(x)

def fmt_xy_axis1(x, pos):
   return r'${:.1f}$'.format(x)

def fmt_xy_axis2(x, pos):
   return r'${:.2f}$'.format(x)

#for i in range(6):
print ("Processing file ",filename)
f = h5py.File(filename,'r')

theta = f["theta"][()]
zeta = f["zeta"][()]
Ntheta = f["Ntheta"][()]
Nzeta = f["Nzeta"][()]

zeta2D, theta2D = numpy.meshgrid(zeta, theta)

zeta2D = zeta2D.transpose()
theta2D = theta2D.transpose()

integerToRepresentTrue = f["integerToRepresentTrue"][()]

if includePhi1InCalculation:
   includePhi1 = f["includePhi1"][()] ##Check that Phi1 is included in SFINCS calculation

   if includePhi1 == integerToRepresentTrue:
      print("INCLUDING Phi1 IN CALCULATION.")
      Phi1Hat = f["Phi1Hat"][()] ## "zeta", "theta", "iteration"
      dPhi1Hatdtheta = f["dPhi1Hatdtheta"][()] ## "zeta", "theta", "iteration"
      dPhi1Hatdzeta = f["dPhi1Hatdzeta"][()] ## "zeta", "theta", "iteration"
   else:
      print("WARNING! YOU WANT TO PLOT FLUXES INCLUDING Phi1, BUT Phi1 IS NOT INCLUDED IN THE SFINCS OUTPUT FILE. MAKING CALCULATION WITHOUT Phi1.")

DHat = f["DHat"][()] ## "zeta", "theta" ##Inverse Jacobian (grad psi dot grad theta cross grad zeta)
VPrimeHat = f["VPrimeHat"][()] ##Scalar Int(1/DHat)
NPeriods = f["NPeriods"][()] ##Scalar Number of identical toroidal periods (e.g. 5 for W7-X, 10 for LHD, 4 for HSX)

BHat = f["BHat"][()] ## "zeta", "theta"
BHat_sub_psi = f["BHat_sub_psi"][()] ## "zeta", "theta"
dBHatdpsiHat = f["dBHatdpsiHat"][()] ## "zeta", "theta"
dBHatdtheta = f["dBHatdtheta"][()] ## "zeta", "theta"
dBHatdzeta = f["dBHatdzeta"][()] ## "zeta", "theta"
psiAHat = f["psiAHat"][()]

##OLD dIotadpsiHat = dIotadpsiN / psiAHat
dIotadpsiHat = f["diotadpsiHat"][()]
dIotadpsiN = dIotadpsiHat * psiAHat

dnHatdpsiHat = f["dnHatdpsiHat"][()] ## "species"
dTHatdpsiHat = f["dTHatdpsiHat"][()] ## "species"
Zs = f["Zs"][()] ## "species"
mHats = f["mHats"][()] ## "species"
THats = f["THats"][()] ## "species"
nHats = f["nHats"][()] ## "species"
B0OverBBar = f["B0OverBBar"][()]
GHat = f["GHat"][()]
IHat = f["IHat"][()]
iota = f["iota"][()]
q = 1.0/iota
nu_n = f["nu_n"][()]
alpha = f["alpha"][()]
iteration = f["NIterations"][()] - 1 #Results from last iteration
rN = f["rN"][()]

dqdpsiHat = - dIotadpsiN / (psiAHat * iota**2)

f.close()

deltaTheta = (numpy.amax(theta) - numpy.amin(theta))/(Ntheta - 1.0)
deltaZeta = (numpy.amax(zeta) - numpy.amin(zeta))/(Nzeta - 1.0)

def flux_surface_average(arg):
   # Add both the parameters and return them."
   fluxSurfAvg = (NPeriods*deltaTheta*deltaZeta*numpy.sum(arg/DHat)) / VPrimeHat
   return fluxSurfAvg;

print ("DHat min: " + str(numpy.amin(DHat)))
print ("DHat max: " + str(numpy.amax(DHat)))
print ("VPrimeHat: " + str(VPrimeHat))
#print ("zeta: " + str(zeta))
#print ("theta: " + str(theta))
#print ("Nzeta: " + str(Nzeta))
#print ("Ntheta: " + str(Ntheta))
print ("deltaZeta: " + str(deltaZeta))
print ("deltaTheta: " + str(deltaTheta))

#print ("VPrimeHat: " + str(VPrimeHat))
#print ("VPrimeHat test: " + str(NPeriods*deltaTheta*deltaZeta*numpy.sum(1.0/DHat)))
#print ("Flux-surface average test: " + str(flux_surface_average(1.0)))
#print ("Flux-surface average test 2: " + str(flux_surface_average(VPrimeHat)))
#print ("<B^2>: " + str(flux_surface_average(BHat**2)))

#sys.exit(0)

##SHIFT BOOZER ANGLES TO [-PI/NPeriods, PI/NPeriods] INTERVAL##
zetaMod = zeta.copy()
zetaMod[zetaMod > numpy.pi/NPeriods] += - 2.0*numpy.pi/NPeriods
thetaMod = theta.copy()
thetaMod[thetaMod > numpy.pi] += - 2.0*numpy.pi

zeta2DMod = zeta2D.copy()
zeta2DMod[zeta2DMod > numpy.pi/NPeriods] += - 2.0*numpy.pi/NPeriods
theta2DMod = theta2D.copy()
theta2DMod[theta2DMod > numpy.pi] += - 2.0*numpy.pi
###############################################################

##PERTURBED POTENTIAL##
#phiTMP = numpy.sqrt((zeta2D*NPeriods - numpy.pi)**2 + (theta2D - numpy.pi)**2)/(numpy.sqrt(2)*numpy.pi)
#phiTMP = numpy.sqrt((numpy.pi - numpy.absolute(zeta2D*NPeriods - phiCenterZeta*NPeriods))**2 + (numpy.pi - numpy.absolute(theta2D - phiCenterTheta))**2)
#phiTMP = (numpy.pi - numpy.absolute(zeta2DMod*NPeriods - phiCenterZeta*NPeriods)) + (numpy.pi - numpy.absolute(theta2DMod - phiCenterTheta))
#phiTMP = numpy.cos((zeta2D*NPeriods - phiCenterZeta*NPeriods)) + numpy.cos((theta2D - phiCenterTheta))
phiTMP = numpy.cos((zeta2DMod*NPeriods - phiCenterZeta*NPeriods)/2.0)**2 + numpy.cos((theta2DMod - phiCenterTheta)/2.0)**2
#phiTMP = numpy.amax(numpy.absolute(phiTMP)) - phiTMP 
phiTMP = phiTMP/numpy.amax(numpy.absolute(phiTMP))
phi2 = phiTMP**phiExponent
#phi2 = numpy.cos(phiTMP)**phiExponent
#phi = numpy.cos(zeta2D*NPeriods + theta2D)**4
#######################


DensityToNuFactor = numpy.absolute((GHat + iota*IHat)*nu_n*Zs**4 / (B0OverBBar*THats**2))

eta = nHats*dTHatdpsiHat / (THats*dnHatdpsiHat)

QuantityForQuasilinearFluxes = 1.0

QuantityForQuasilinearFluxes += TermdBdalphaFactor*((2.0*BHat_sub_psi[:,:] /((q*GHat + IHat)*BHat[:,:])) * (nHats[species] / dnHatdpsiHat[species]) * (q*dBHatdzeta[:,:] + dBHatdtheta[:,:]))

QuantityForQuasilinearFluxes += - TermdBdpsifactor*((2.0/BHat[:,:]) * (nHats[species] / dnHatdpsiHat[species]) * dBHatdpsiHat[:,:])

QuantityForQuasilinearFluxes += TermBShearFactor *((2.0*theta2DMod[:,:]*dqdpsiHat /((q*GHat + IHat)*BHat[:,:])) * (nHats[species] / dnHatdpsiHat[species]) * (- IHat*dBHatdzeta[:,:] + GHat*dBHatdtheta[:,:]))

if includePhi1InCalculation and includePhi1:

   QuantityForQuasilinearFluxesNoPhi1 = QuantityForQuasilinearFluxes.copy()
   
   QuantityForQuasilinearFluxes +=  TermPhi1Factor*(eta[species]*Zs[species]*alpha*Phi1Hat[:,:,iteration]/THats[species])

   QuantityForQuasilinearFluxes += TermdPhi1dalphaFactor*((BHat_sub_psi[:,:] /(q*GHat + IHat)) * (Zs[species]*alpha/THats[species]) * (nHats[species] / dnHatdpsiHat[species]) * (q*dPhi1Hatdzeta[:,:,iteration] + dPhi1Hatdtheta[:,:,iteration]))

   ##TEST##
   #temporary = TermPhi1ShearFactor*((theta2DMod[:,:]*dqdpsiHat/(q*GHat + IHat)) * (Zs[species]*alpha/THats[species]) * (nHats[species] / dnHatdpsiHat[species]) * (- IHat*dPhi1Hatdzeta[:,:,iteration] + GHat*dPhi1Hatdtheta[:,:,iteration]))
   #print ("temporary[0,0]: " + str(temporary[0,0]))
   #print ("dPhi1Hatdtheta[0,0,iteration]: " + str(dPhi1Hatdtheta[0,0,iteration]))
   #print ("temporary[1,0]: " + str(temporary[1,0]))
   #print ("dPhi1Hatdtheta[1,0,iteration]: " + str(dPhi1Hatdtheta[1,0,iteration]))
   #print ("temporary[-1,0]: " + str(temporary[-1,0]))
   #print ("dPhi1Hatdtheta[-1,0,iteration]: " + str(dPhi1Hatdtheta[-1,0,iteration]))
   #print ("temporary[0,1]: " + str(temporary[0,1]))
   #print ("dPhi1Hatdtheta[0,1,iteration]: " + str(dPhi1Hatdtheta[0,1,iteration]))
   #print ("temporary[0,-1]: " + str(temporary[0,-1]))
   #print ("dPhi1Hatdtheta[0,-1,iteration]: " + str(dPhi1Hatdtheta[0,-1,iteration]))
   #print ("temporary[0,2]: " + str(temporary[0,2]))
   #print ("temporary[0,-2]: " + str(temporary[0,-2]))
   #print ("temporary[1,1]: " + str(temporary[1,1]))
   #print ("temporary[-1,-1]: " + str(temporary[-1,-1]))
   
   #print ("temporary[0,3]: " + str(temporary[0,3]))
   #print ("temporary[0,-3]: " + str(temporary[0,-3]))
   #print ("temporary[2,1]: " + str(temporary[2,1]))
   #print ("temporary[-2,-1]: " + str(temporary[-2,-1]))
   #print ("temporary: " + str(temporary))
   #temporary[10:-10,:] = 0.0
   #temporary[:,10:-10] = 0.0
   #print ("temporary: " + str(temporary))
   #print (r'< |Phi|^2 temporary > = ' + str(flux_surface_average(phi2[:,:]*temporary[:,:])))
   #print (r'< |Phi|^2 temporary > /< |Phi|^2 > = ' + str(flux_surface_average(phi2[:,:]*temporary[:,:]) / flux_surface_average(phi2[:,:])))
   ########

   QuantityForQuasilinearFluxes += TermPhi1ShearFactor*((theta2DMod[:,:]*dqdpsiHat/(q*GHat + IHat)) * (Zs[species]*alpha/THats[species]) * (nHats[species] / dnHatdpsiHat[species]) * (- IHat*dPhi1Hatdzeta[:,:,iteration] + GHat*dPhi1Hatdtheta[:,:,iteration]))
   
   QuantityForQuasilinearFluxes = numpy.exp(- TermPhi1Exponential * Zs[species]*alpha*Phi1Hat[:,:,iteration]/THats[species]) * QuantityForQuasilinearFluxes

if includePhi1InCalculation and includePhi1:
   zLabel = r'$\exp[- Z e \Phi_1 / T_z] \{ 1 + \eta_z Z e \Phi_1 / T_z  - q^{-1} \omega_E / \omega_{\ast z} - q^{-1} \omega_B / \omega_{\ast z}\}$'
   zLabelExtra = r'$1 - q^{-1} \omega_B / \omega_{\ast z}$'
else:
   zLabel = r'$1 - q^{-1} \omega_B / \omega_{\ast z}$'

print ("zeta: " + str(zeta))
print ("zeta mod: " + str(zetaMod))
print ("theta: " + str(theta))
print ("theta mod: " + str(thetaMod))
print ("zeta2D: " + str(zeta2D))
print ("theta2D: " + str(theta2D))
print ("zeta2DMod: " + str(zeta2DMod))
print ("theta2DMod: " + str(theta2DMod))

print ("zeta shape: " + str(zeta.shape))
print ("theta shape: " + str(theta.shape))
print ("zeta2D shape: " + str(zeta2D.shape))
print ("theta2D shape: " + str(theta2D.shape))
print ("zeta2DMod shape: " + str(zeta2DMod.shape))
print ("theta2DMod shape: " + str(theta2DMod.shape))
print ("BHat shape: " + str(BHat.shape))

print ("*********************************") 
print ("READ QUANTITIES: ") 
print ("*********************************") 
print ("Zs: " + str(Zs))
print ("Z: " + str(Zs[species]))
print ("mHats: " + str(mHats))
print ("THats: " + str(THats))
print ("nHats: " + str(nHats))
print ("dTHatdpsiHat: " + str(dTHatdpsiHat))
print ("dnHatdpsiHat: " + str(dnHatdpsiHat))
print ("eta: " + str(eta))
print ("B0OverBBar: " + str(B0OverBBar))
print ("GHat: " + str(GHat))
print ("IHat: " + str(IHat))
print ("iota: " + str(iota))
print ("dIotadpsiHat: " + str(dIotadpsiHat))
print ("q: " + str(q))
print ("dqdpsiHat: " + str(dqdpsiHat))
print ("dqdpsiN: " + str(dqdpsiHat*psiAHat))
print ("psiAHat: " + str(psiAHat))
print ("nu_n: " + str(nu_n))
print ("alpha: " + str(alpha))
print ("NuPrime_aa: " + str(DensityToNuFactor*nHats))
print ("NuPrime_a = SUM_b (NuPrime_ab): " + str(DensityToNuFactor*nHats * sum(Zs**2 * nHats) / (Zs**2 * nHats)))


mindexZeta = numpy.argmin(numpy.absolute(zeta - phiCenterZeta))
mindexTheta = numpy.argmin(numpy.absolute(theta - phiCenterTheta))

print ("mindexZeta: " + str(mindexZeta))
print ("mindexTheta: " + str(mindexTheta))
print ("|phi|^2 zeta center: " + str(zeta[mindexZeta]))
print ("|phi|^2 theta center: " + str(theta[mindexTheta]))

if includePhi1InCalculation and includePhi1:
   print ("Phi1[zeta = " + str(zeta[mindexZeta]) + ", theta = " + str(theta[mindexTheta]) +  "]: " + str(TermPhi1Factor*1000*Phi1Hat[mindexZeta ,mindexTheta, iteration]))
   print ("dPhi1Hatdtheta[zeta = " + str(zeta[mindexZeta]) + ", theta = " + str(theta[mindexTheta]) +  "]: " + str(TermPhi1Factor*1000*dPhi1Hatdtheta[mindexZeta, mindexTheta, iteration]))
   print ("dPhi1Hatdtheta[zeta = " + str(zeta[mindexZeta]) + ", theta = " + str(theta[mindexTheta + 1]) +  "]: " + str(TermPhi1Factor*1000*dPhi1Hatdtheta[mindexZeta, mindexTheta + 1, iteration]))
   print ("dPhi1Hatdtheta[zeta = " + str(zeta[mindexZeta]) + ", theta = " + str(theta[mindexTheta - 1]) +  "]: " + str(TermPhi1Factor*1000*dPhi1Hatdtheta[mindexZeta, mindexTheta - 1, iteration]))
   print ("etaz * Z * alpha * Phi1[zeta = " + str(zeta[mindexZeta]) + ", theta = " + str(theta[mindexTheta]) +  "] / Tz: " + str(TermPhi1Factor*(eta[species]*Zs[species]*alpha*Phi1Hat[mindexZeta, mindexTheta, iteration]/THats[species])))
   print ("eta * Zs * alpha * Phi1[zeta = " + str(zeta[mindexZeta]) + ", theta = " + str(theta[mindexTheta]) +  "] / Ts: " + str(TermPhi1Factor*(eta*Zs*alpha*Phi1Hat[mindexZeta, mindexTheta, iteration]/THats)))
print ("*********************************")

print ("theta max: " + str(numpy.amax(theta)))
print ("zeta max: " + str(numpy.amax(zeta)))

zMinData = zFactor*numpy.amin(QuantityForQuasilinearFluxes[:,:])
zMaxData = zFactor*numpy.amax(QuantityForQuasilinearFluxes[:,:])

if includePhi1InCalculation and includePhi1:
   zMinDataNoPhi1 = zFactor*numpy.amin(QuantityForQuasilinearFluxesNoPhi1[:,:])
   zMaxDataNoPhi1 = zFactor*numpy.amax(QuantityForQuasilinearFluxesNoPhi1[:,:])

   zMinData = numpy.minimum(zMinData, zMinDataNoPhi1)
   zMaxData = numpy.maximum(zMaxData, zMaxDataNoPhi1)

if zMin == None:
   zMin = zMinData
if zMax == None:
   zMax = zMaxData

if zMax < zMin:
   print ("Error, zMax < zMin")
   sys.exit(1)

print ("zMin = " + str(zMin))
print ("zMinData = " + str(zMinData))
print ("zMax = " + str(zMax))
print ("zMaxData = " + str(zMaxData))

#sys.exit(0)

print (r'< |phi|^2 Q_{w/ Phi1} > /< |Phi|^2 > = ' + str(flux_surface_average(phi2[:,:]*QuantityForQuasilinearFluxes[:,:]) / flux_surface_average(phi2[:,:])))

if includePhi1InCalculation and includePhi1:
   print (r'< |phi|^2 Q_{w/o Phi1} > /< |Phi|^2 > = ' + str(flux_surface_average(phi2[:,:]*QuantityForQuasilinearFluxesNoPhi1[:,:]) / flux_surface_average(phi2[:,:])))

delta = (numpy.amax(QuantityForQuasilinearFluxes[:,:]) - numpy.amin(QuantityForQuasilinearFluxes[:,:])) / numLevels

if not ContourLevels:
   ContourLevels = numpy.arange(numpy.amin(QuantityForQuasilinearFluxes[:,:]), numpy.amax(QuantityForQuasilinearFluxes[:,:]) + delta/2.0, delta)
   ContourLevels = zFactor*ContourLevels

#############################################################################################################################################
##PLOT##
#############################################################################################################################################

##TRANSFORM TO PLOT SYMETRIC INTERVAL AROUND 0 IN ZETA, THETA##
zetaSortIndices = numpy.argsort(zetaMod)
thetaSortIndices = numpy.argsort(thetaMod)

zeta2DSortIndices, theta2DSortIndices = numpy.meshgrid(zetaSortIndices, thetaSortIndices)

zeta2DSortIndices = zeta2DSortIndices.transpose()
theta2DSortIndices = theta2DSortIndices.transpose()

print("zetaSortIndices: " + str(zetaSortIndices))
print("zetaMod[zetaSortIndices]: " + str(zetaMod[zetaSortIndices]))
print("thetaSortIndices: " + str(thetaSortIndices))
print("thetaMod[thetaSortIndices]: " + str(thetaMod[thetaSortIndices]))

print("QuantityForQuasilinearFluxes.shape: " + str(QuantityForQuasilinearFluxes.shape))
print("QuantityForQuasilinearFluxes[zeta2DSortIndices, theta2DSortIndices].shape: " + str(QuantityForQuasilinearFluxes[zeta2DSortIndices, theta2DSortIndices].shape))
###############################################################

FigNum = 1

ax = plt.subplot(numRows,numCols,FigNum)

#Plot1 = plt.contourf(zeta,theta,zFactor*QuantityForQuasilinearFluxes[:,:].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'), vmin=zMin, vmax=zMax)
Plot1 = plt.contourf(zetaMod[zetaSortIndices], thetaMod[thetaSortIndices], zFactor*QuantityForQuasilinearFluxes[zeta2DSortIndices, theta2DSortIndices].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'), vmin=zMin, vmax=zMax)
#Plot1 = plt.contourf(zetaMod,thetaMod,zFactor*QuantityForQuasilinearFluxes[:,:].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'))

Plot2 = plt.contour(Plot1, levels=ContourLevels, colors='k')

##TEST Plot1 = plt.contourf(zeta,theta,theta2D[:,:].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'))

##ADD CROSS AT PEAK##
#if phiMaxShow:
#   #Plot1point = plt.contour(numpy.array([phiCenterZeta]), numpy.array([phiCenterTheta]), numpy.array([1.0]), colors='k')
#   ax.text(phiCenterZeta, phiCenterTheta, r"$\max$",  color='black', zorder=10, fontsize=70, fontweight='heavy')
#   #phiCenterZeta
#   #phiCenterTheta
#####################

#plt.xlabel(r'$\zeta$' + ' [rad]')
plt.xlabel(r'$\zeta$' + " " + r'$\mathrm{[rad]}$')
#plt.ylabel(r'$\theta$'+ ' [rad]')
plt.ylabel(r'$\theta$'+ " " + r'$\mathrm{[rad]}$')
#plt.zlabel(r'$\Phi_1$'+ ' [V]')

if withTickLabels:
   plt.xticks([0,max(zeta)/4,max(zeta)/2,3*max(zeta)/4,max(zeta)])
   plt.yticks([0.0,max(theta)/4,max(theta)/2,3*max(theta)/4,max(theta)])
   plt.gca().axes.xaxis.set_ticklabels(xAxisTicks)
   plt.gca().axes.yaxis.set_ticklabels(yAxisTicks)

#plt.gca().axes.xaxis.set_label_coords(0.5,-0.09)
#plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)
plt.gca().axes.xaxis.set_label_coords(0.5,-0.09)
plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)

#ax.xaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis))
#ax.xaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))
#ax.yaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis)) 
#ax.yaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))

#plt.tight_layout()

figTitle = r''

#if show_rN:
#    #plt.title('r/a = '+str(rN))
#    figTitle = figTitle + 'r/a = '+str(rN) + ', '

if includePhi1InCalculation and includePhi1:
   #plt.title('w/ $\Phi_1$ ')
   figTitle = figTitle + '$Q$ w/ $\Phi_1$'
else:
   #plt.title('w/o $\Phi_1$ ')
   figTitle = figTitle + '$Q$ w/o $\Phi_1$'

plt.title(figTitle)

#cbar = plt.colorbar(Plot1, label=r'$\Phi_1$'+ ' [V]', ticks=ContourLevels)
#cbar = plt.colorbar(Plot1, label=r'$\Phi_1$'+ ' [V]', ticks=Plot1.levels[::2])
#cbar.add_lines(Plot2)
cbar = plt.colorbar(Plot1, format=ticker.FuncFormatter(fmt_xy_axis2), ticks=ContourLevels)
#cbar.ax.set_ylabel(zLabel, rotation=270, labelpad=35)

#with warnings.catch_warnings():
#    warnings.simplefilter("always")
#plt.clabel(Plot2, fmt='%2.1f', colors='k', fontsize=14)
plt.clabel(Plot2, fmt=ticker.FuncFormatter(fmt_xy_axis1), colors='k', fontsize=30, inline=False)

plt.subplots_adjust(wspace=0.15, top=0.95)
#plt.subplots_adjust(left=LeftMargin, right=RightMargin, top=TopMargin, bottom=BottomMargin)

######################################################################################################

if includePhi1InCalculation and includePhi1:

   FigNum += 1

   ax = plt.subplot(numRows,numCols,FigNum)

   #Plot1extra = plt.contourf(zeta,theta,zFactor*QuantityForQuasilinearFluxesNoPhi1[:,:].transpose(), numContours, cmap=plt.get_cmap('gist_rainbow'), vmin=zMin, vmax=zMax)
   Plot1extra = plt.contourf(zetaMod[zetaSortIndices], thetaMod[thetaSortIndices],zFactor*QuantityForQuasilinearFluxesNoPhi1[zeta2DSortIndices, theta2DSortIndices].transpose(), numContours, cmap=plt.get_cmap('gist_rainbow'), vmin=zMin, vmax=zMax)

   Plot2extra = plt.contour(Plot1extra, levels=ContourLevels, colors='k', vmin=zMin, vmax=zMax)

   ##ADD CROSS AT PEAK##
   #if phiMaxShow:
   #   ax.text(phiCenterZeta, phiCenterTheta, r"$\max$",  color='black', zorder=10, fontsize=70, fontweight='heavy')
   #####################

   plt.xlabel(r'$\zeta$' + " " + r'$\mathrm{[rad]}$')
   plt.ylabel(r'$\theta$'+ " " + r'$\mathrm{[rad]}$')

   if withTickLabels:
      plt.xticks([0,max(zeta)/4,max(zeta)/2,3*max(zeta)/4,max(zeta)])
      plt.yticks([0.0,max(theta)/4,max(theta)/2,3*max(theta)/4,max(theta)])
      plt.gca().axes.xaxis.set_ticklabels(xAxisTicks)
      plt.gca().axes.yaxis.set_ticklabels(yAxisTicks)

   plt.gca().axes.xaxis.set_label_coords(0.5,-0.09)
   plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)


   cbarextra = plt.colorbar(Plot1extra, format=ticker.FuncFormatter(fmt_xy_axis2), ticks=ContourLevels)
   #cbarextra.set_clim(vmin=zMin, vmax=zMax)
   #cbarextra.ax.set_ylabel(zLabelExtra, rotation=270, labelpad=35)
   plt.title(r'$Q$ w/o $\Phi_1$')

   plt.clabel(Plot2extra, fmt=ticker.FuncFormatter(fmt_xy_axis1), colors='k', fontsize=30, inline=False)

   plt.subplots_adjust(wspace=0.15, top=0.95)

######################################################################################################

FigNum += 1

ax2 = plt.subplot(numRows,numCols,FigNum)

#Plot3 = plt.contourf(zeta,theta,phi2[:,:].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'))
Plot3 = plt.contourf(zetaMod[zetaSortIndices], thetaMod[thetaSortIndices],phi2[zeta2DSortIndices, theta2DSortIndices].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'))

plt.xlabel(r'$\zeta$' + " " + r'$\mathrm{[rad]}$')
plt.ylabel(r'$\theta$'+ " " + r'$\mathrm{[rad]}$')

if withTickLabels:
   plt.xticks([0,max(zeta)/4,max(zeta)/2,3*max(zeta)/4,max(zeta)])
   plt.yticks([0.0,max(theta)/4,max(theta)/2,3*max(theta)/4,max(theta)])
   plt.gca().axes.xaxis.set_ticklabels(xAxisTicks)
   plt.gca().axes.yaxis.set_ticklabels(yAxisTicks)

plt.gca().axes.xaxis.set_label_coords(0.5,-0.09)
plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)

cbar2 = plt.colorbar(Plot3, format=ticker.FuncFormatter(fmt_xy_axis))

#cbar2.ax.set_ylabel(r'$|\phi|^2$', rotation=0, labelpad=35)
plt.title(r'$|\phi|^2$' + r' (max in ' r'$[\zeta$ = ' + str(r'${}$'.format('{:1.2f}'.format(zeta[mindexZeta]))) + r', $\theta$ = ' + str(r'${}$'.format('{:1.2f}'.format(theta[mindexTheta]))) +  r'$]$)')

#plt.clabel(Plot3, fmt=ticker.FuncFormatter(fmt_cbar), colors='k', fontsize=18, inline=False)

######################################################################################################

if includePhi1InCalculation and includePhi1:

   FigNum += 1

   ax3 = plt.subplot(numRows,numCols,FigNum)

   #Plot4 = plt.contourf(zeta,theta, Phi1Hat[:,:,iteration].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'))
   Plot4 = plt.contourf(zetaMod[zetaSortIndices], thetaMod[thetaSortIndices], Phi1Hat[zeta2DSortIndices, theta2DSortIndices, iteration].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'))

   ##ADD CROSS AT PEAK##
   #if phiMaxShow:
   #   ax3.text(phiCenterZeta, phiCenterTheta, "x",  color='black', zorder=10, fontsize=50)
   #####################

   plt.xlabel(r'$\zeta$' + " " + r'$\mathrm{[rad]}$')
   plt.ylabel(r'$\theta$'+ " " + r'$\mathrm{[rad]}$')

   if withTickLabels:
      plt.xticks([0,max(zeta)/4,max(zeta)/2,3*max(zeta)/4,max(zeta)])
      plt.yticks([0.0,max(theta)/4,max(theta)/2,3*max(theta)/4,max(theta)])
      plt.gca().axes.xaxis.set_ticklabels(xAxisTicks)
      plt.gca().axes.yaxis.set_ticklabels(yAxisTicks)

   plt.gca().axes.xaxis.set_label_coords(0.5,-0.09)
   plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)

   cbar3 = plt.colorbar(Plot4, format=ticker.FuncFormatter(fmt_xy_axis2))
   
   #cbar3.ax.set_ylabel(r'$\Phi_1~[\mathrm{kV}]$', rotation=270, labelpad=35)
   plt.title(r'$\Phi_1~[\mathrm{kV}]$')
   
   #plt.clabel(Plot3, fmt=ticker.FuncFormatter(fmt_cbar), colors='k', fontsize=18, inline=False)

######################################################################################################

if includePhi1InCalculation and includePhi1:

   FigNum += 1

   ax4 = plt.subplot(numRows,numCols,FigNum)

   #Plot5 = plt.contourf(zeta,theta, dPhi1Hatdtheta[:,:,iteration].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'))
   Plot5 = plt.contourf(zetaMod[zetaSortIndices], thetaMod[thetaSortIndices], dPhi1Hatdtheta[zeta2DSortIndices, theta2DSortIndices, iteration].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'))

   ##ADD CROSS AT PEAK##
   #if phiMaxShow:
   #   ax4.text(phiCenterZeta, phiCenterTheta, "x",  color='black', zorder=10, fontsize=50)
   #####################

   plt.xlabel(r'$\zeta$' + " " + r'$\mathrm{[rad]}$')
   plt.ylabel(r'$\theta$'+ " " + r'$\mathrm{[rad]}$')

   if withTickLabels:
      plt.xticks([0,max(zeta)/4,max(zeta)/2,3*max(zeta)/4,max(zeta)])
      plt.yticks([0.0,max(theta)/4,max(theta)/2,3*max(theta)/4,max(theta)])
      plt.gca().axes.xaxis.set_ticklabels(xAxisTicks)
      plt.gca().axes.yaxis.set_ticklabels(yAxisTicks)

   plt.gca().axes.xaxis.set_label_coords(0.5,-0.09)
   plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)

   cbar4 = plt.colorbar(Plot5, format=ticker.FuncFormatter(fmt_xy_axis2))
   
   #cbar4.ax.set_ylabel(r'$\frac{d \Phi_1}{d \theta}~[\mathrm{kV/rad}]$', rotation=270, labelpad=35)
   plt.title(r'$d \Phi_1 / d \theta~[\mathrm{kV/rad}]$')
   
   #plt.clabel(Plot3, fmt=ticker.FuncFormatter(fmt_cbar), colors='k', fontsize=18, inline=False)   

######################################################################################################

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

InfoText = r''
if show_rN:
   InfoText = InfoText + r'r/a = ' + str(rN) + '\n'
   #InfoText = InfoText + '$\;$ \n'

InfoText = InfoText + r'$Z = $ ' + str(r'${}$'.format('{:1.0f}'.format(Zs[species]))) + r',   $T_z = $ ' + str(r'${}$'.format(InfoBoxFormat.format(THats[species]))) + r' $\mathrm{keV}$' + r',   $\eta_z = $ ' + str(r'${}$'.format(InfoBoxFormat.format(eta[species]))) + r',   $q = $ ' + str(r'${}$'.format(InfoBoxFormat.format(q))) + r',   $d q / d \psi_N = $ ' + str(r'${}$'.format(InfoBoxFormat.format(dqdpsiHat*psiAHat))) + '\n'

InfoText = InfoText + r'$\mathbf{Q}$ \textbf{=} ' + zLabel + '\n'

if includePhi1InCalculation and includePhi1:

   InfoText = InfoText + r'$< |\phi|^2 \, \mathbf{Q_{\mathbf{w/} \; \Phi_1}} > /< |\phi|^2 >$ \textbf{=} ' + str(r'${}$'.format(InfoBoxFormat.format(flux_surface_average(phi2[:,:]*QuantityForQuasilinearFluxes[:,:]) / flux_surface_average(phi2[:,:])))) + '\n'

   InfoText = InfoText + r'$< |\phi|^2 \, \mathbf{Q_{\mathbf{w/o} \; \Phi_1}} > /< |\phi|^2 >$ \textbf{=} ' + str(r'${}$'.format(InfoBoxFormat.format(flux_surface_average(phi2[:,:]*QuantityForQuasilinearFluxesNoPhi1[:,:]) / flux_surface_average(phi2[:,:])))) + '\n'

   InfoText = InfoText + r'$\mathbf{Q_{\mathbf{w/} \; \Phi_1}} \,$' + r'\textbf{(}$\mathbf{\zeta}$ \textbf{=} ' + str(r'${}$'.format(InfoBoxFormat.format(zeta[mindexZeta]))) + r'\textbf{,} $\mathbf{\theta}$ \textbf{=} ' + str(r'${}$'.format(InfoBoxFormat.format(theta[mindexTheta]))) +  r'\textbf{)} \textbf{=} ' + str(r'${}$'.format(InfoBoxFormat.format(QuantityForQuasilinearFluxes[mindexZeta, mindexTheta]))) + '\n'

   InfoText = InfoText + r'$\mathbf{Q_{\mathbf{w/o} \; \Phi_1}} \,$' + r'\textbf{(}$\mathbf{\zeta}$ \textbf{=} ' + str(r'${}$'.format(InfoBoxFormat.format(zeta[mindexZeta]))) + r'\textbf{,} $\mathbf{\theta}$ \textbf{=} ' + str(r'${}$'.format(InfoBoxFormat.format(theta[mindexTheta]))) +  r'\textbf{)} \textbf{=} ' + str(r'${}$'.format(InfoBoxFormat.format(QuantityForQuasilinearFluxesNoPhi1[mindexZeta, mindexTheta]))) + '\n'

   InfoText = InfoText + r'$\Phi_1 \,$' + r'$(\zeta$ = ' + str(r'${}$'.format(InfoBoxFormat.format(zeta[mindexZeta]))) + r', $\theta$ = ' + str(r'${}$'.format(InfoBoxFormat.format(theta[mindexTheta]))) +  r'$)$ = ' + str(r'${}$'.format(InfoBoxFormat.format(TermPhi1Factor*1000*Phi1Hat[mindexZeta ,mindexTheta, iteration]))) + r' $\mathrm{V}$' + '\n'

   InfoText = InfoText + r'$d \Phi_1 / d \theta \,$' + r'$(\zeta$ = ' + str(r'${}$'.format(InfoBoxFormat.format(zeta[mindexZeta]))) + r', $\theta$ = ' + str(r'${}$'.format(InfoBoxFormat.format(theta[mindexTheta]))) +  r'$)$ = ' + str(r'${}$'.format(InfoBoxFormat.format(TermPhi1Factor*1000*dPhi1Hatdtheta[mindexZeta, mindexTheta, iteration]))) + r' $\mathrm{V/rad}$' + '\n'

   InfoText = InfoText + r'$\eta_z \, Z \, \Phi_1 \,$' + r'$(\zeta$ = ' + str(r'${}$'.format(InfoBoxFormat.format(zeta[mindexZeta]))) + r', $\theta$ = ' + str(r'${}$'.format(InfoBoxFormat.format(theta[mindexTheta]))) +  r'$) \, / \, T_z$ = ' + str(r'${}$'.format(InfoBoxFormat.format(TermPhi1Factor*eta[species]*Zs[species]*alpha*Phi1Hat[mindexZeta, mindexTheta, iteration]/THats[species]))) #+ '\n'

else:
   InfoText = InfoText + r'$< |\phi|^2 \, Q > /< |\phi|^2 >$ = ' + str(r'${}$'.format(InfoBoxFormat.format(flux_surface_average(phi2[:,:]*QuantityForQuasilinearFluxes[:,:]) / flux_surface_average(phi2[:,:])))) + '\n'


InfoText = InfoText + r'.'


if InfoBox:
   plt.text(InfoBoxXcoord, InfoBoxYcoord, InfoText, fontsize=InfoBoxLabelSize, linespacing=InfoBoxLineSpacing)

print (QuantityForQuasilinearFluxes.shape)

print ("*********************************")

if makePDF:
    print ("Saving PDF")

    if len(sys.argv)>2 : #Use the substituted name as file name
       print ("Writing plot to " + os.getcwd() + "/" + sys.argv[2] + ".pdf.")
       plt.savefig(sys.argv[2] + ".pdf", orientation = 'landscape', papertype='letter')
    else :
       head, tail = os.path.split(inspect.getfile(inspect.currentframe()))
       print ("Writing plot to " + os.getcwd() + "/" + tail + ".pdf.")
       plt.savefig(tail+'.pdf', orientation = 'landscape', papertype='letter')
else:
    plt.show()
