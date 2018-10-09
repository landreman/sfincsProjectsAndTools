#!/usr/bin/env python

import matplotlib
#import matplotlib.pyplot as plt
import h5py
import numpy
import os, sys, inspect
import warnings
import matplotlib.ticker as ticker

#show_rN = True
show_rN = False

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

FigSize = (12,10)

font = {'size':25}
matplotlib.rc('font', **font)
matplotlib.rc('lines',markeredgewidth=0,markersize=3,linewidth=2.5)
matplotlib.rc('axes',linewidth=1.5)

matplotlib.rcParams['mathtext.default'] = 'it'
matplotlib.rcParams['text.usetex'] = True

zFactor = 1
withTickLabels = False
##W7-X##
#xAxisTicks = [r'$0$', r'$\pi/10$', r'$2\pi/10$', r'$3\pi/10$', r'$4\pi/10$']
##LHD
xAxisTicks = [r'$0$', r'$\pi/20$', r'$2\pi/20$', r'$3\pi/20$', r'$4\pi/20$']

yAxisTicks = [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']

fig = plt.figure(figsize=FigSize)
fig.patch.set_facecolor('white')
numRows = 1
numCols = 1
#iteration = 0
numContours = 100
#ContourLevels = [-3.0, -1.5, 0.0, 1.5, 3.0, 4.5, 6.0]
numLevels = 5

species = 4 - 1 ##Species index starts from 0

includePhi1InCalculation = True

dIotadpsiN = 0.0333098952242 ##Manual input

TermPhi1Factor = 1.0
TermdPhi1dalphaFactor = 1.0
TermPhi1ShearFactor = 1.0
TermdBdalphaFactor = 1.0
TermdBdpsifactor = 1.0
TermBShearFactor = 1.0

zLabel = r'$1 + Ze\Phi_1 / T_z - \omega_E / \omega_{\ast z} - v_z^2 (\omega_B + \omega_\kappa) / \omega_{\ast z}$'

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

#for i in range(6):
print ("Processing file ",filename)
f = h5py.File(filename,'r')

theta = f["theta"][()]
zeta = f["zeta"][()]

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


BHat = f["BHat"][()] ## "zeta", "theta"
BHat_sub_psi = f["BHat_sub_psi"][()] ## "zeta", "theta"
dBHatdpsiHat = f["dBHatdpsiHat"][()] ## "zeta", "theta"
dBHatdtheta = f["dBHatdtheta"][()] ## "zeta", "theta"
dBHatdzeta = f["dBHatdzeta"][()] ## "zeta", "theta"
psiAHat = f["psiAHat"][()]
dIotadpsiHat = dIotadpsiN / psiAHat

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
nu_n = f["nu_n"][()]
alpha = f["alpha"][()]
iteration = f["NIterations"][()] - 1 #Results from last iteration
rN = f["rN"][()]
f.close()

DensityToNuFactor = numpy.absolute((GHat + iota*IHat)*nu_n*Zs**4 / (B0OverBBar*THats**2))

eta = nHats*dTHatdpsiHat / (THats*dnHatdpsiHat)

QuantityForQuasilinearFluxes = 1.0 + TermdBdalphaFactor*((2.0*BHat_sub_psi[:,:] /((GHat + iota*IHat)*BHat[:,:])) * (nHats[species] / dnHatdpsiHat[species]) * (dBHatdzeta[:,:] + iota*dBHatdtheta[:,:])) - TermdBdpsifactor*((2.0/BHat[:,:]) * (nHats[species] / dnHatdpsiHat[species]) * dBHatdpsiHat[:,:])

if includePhi1InCalculation and includePhi1:
   QuantityForQuasilinearFluxes = QuantityForQuasilinearFluxes + TermPhi1Factor*(eta[species]*Zs[species]*alpha*Phi1Hat[:,:,iteration]/THats[species]) + TermdPhi1dalphaFactor*((BHat_sub_psi[:,:] /(GHat + iota*IHat)) * (Zs[species]*alpha/THats[species]) * (nHats[species] / dnHatdpsiHat[species]) * (dPhi1Hatdzeta[:,:,iteration] + iota*dPhi1Hatdtheta[:,:,iteration]))

print ("*********************************") 
print ("READ QUANTITIES: ") 
print ("*********************************") 
print ("Zs: " + str(Zs))
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
print ("psiAHat: " + str(psiAHat))
print ("nu_n: " + str(nu_n))
print ("alpha: " + str(alpha))
print ("NuPrime_aa: " + str(DensityToNuFactor*nHats))
print ("NuPrime_a = SUM_b (NuPrime_ab): " + str(DensityToNuFactor*nHats * sum(Zs**2 * nHats) / (Zs**2 * nHats))) 
print ("*********************************")

print ("theta max: " + str(numpy.amax(theta)))
print ("zeta max: " + str(numpy.amax(zeta)))

zMinData = zFactor*numpy.amin(QuantityForQuasilinearFluxes[:,:])
zMaxData = zFactor*numpy.amax(QuantityForQuasilinearFluxes[:,:])
print ("zMin = " + str(zMinData))
print ("zMax = " + str(zMaxData))


delta = (numpy.amax(QuantityForQuasilinearFluxes[:,:]) - numpy.amin(QuantityForQuasilinearFluxes[:,:])) / numLevels
ContourLevels = numpy.arange(numpy.amin(QuantityForQuasilinearFluxes[:,:]), numpy.amax(QuantityForQuasilinearFluxes[:,:]) + delta/2.0, delta)
ContourLevels = zFactor*ContourLevels
    
ax = plt.subplot(numRows,numCols,1)

Plot1 = plt.contourf(zeta,theta,zFactor*QuantityForQuasilinearFluxes[:,:].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'))

Plot2 = plt.contour(Plot1,levels=ContourLevels, colors='k')

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

#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

if show_rN:
    plt.title('rN = '+str(rN))

#cbar = plt.colorbar(Plot1, label=r'$\Phi_1$'+ ' [V]', ticks=ContourLevels)
#cbar = plt.colorbar(Plot1, label=r'$\Phi_1$'+ ' [V]', ticks=Plot1.levels[::2])
#cbar.add_lines(Plot2)
cbar = plt.colorbar(Plot1, format=ticker.FuncFormatter(fmt_cbar), ticks=ContourLevels)
cbar.ax.set_ylabel(zLabel, rotation=270, labelpad=10)

#with warnings.catch_warnings():
#    warnings.simplefilter("always")
#plt.clabel(Plot2, fmt='%2.1f', colors='k', fontsize=14)
plt.clabel(Plot2, fmt=ticker.FuncFormatter(fmt_cbar), colors='k', fontsize=18, inline=False)

#plt.subplots_adjust(wspace=0.27)

print (QuantityForQuasilinearFluxes.shape)

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
