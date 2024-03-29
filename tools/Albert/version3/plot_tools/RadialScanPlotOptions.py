#!/usr/bin/env python

#####################################
####PLOT OPTIONS FOR RADIAL SCANS####
#####################################

##CAN BE OVERWRITTEN IN ORIGINAL FILE

import matplotlib
import numpy 

#matplotlib.rc_context(rc={
#'lines.linewidth': 1.5,
#'lines.dashed_pattern' : [2.8, 1.2],
#'lines.dashdot_pattern' : [4.8, 1.2, 0.8, 1.2],
#'lines.dotted_pattern' : [1.1, 1.1],
#'lines.scale_dashes': True})

#FigSize = (12,10) #LHD113208t4640, W7-X NBI case, Er scan, Flux scan, #EPS poster Er, Collisionality, Ar fluxes  #W7X_180919.055 fluxes/profiles
#FigSize = (13.5,10) #Phi1 plot
#FigSize = (13.0,10.5) #Phi1 plot
FigSize = (12,12) #Phi1 plot without colorbar
#FigSize = (20,10) #D and V plot
#FigSize = (20,20) #D and V extended plot
font = {'size':35} #Phi1 plot without colorbar
#font = {'size':38} #EPS poster
#font = {'size':50} #EPS poster
matplotlib.rc('font', **font) 
matplotlib.rc('lines',markeredgewidth=0,markersize=3,linewidth=2.5)
matplotlib.rc('axes',linewidth=1.5)
#matplotlib.rc('axes',labelpad=0.0)
#scilimits=(0,0)
matplotlib.rc('axes.formatter', limits=(0,2))
#matplotlib.rc('axes.Axes.ticklabel_format',style='sci', axis='y', scilimits=(0,0))

matplotlib.rcParams['mathtext.default'] = 'it'
matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.preamble'].append(r'\usepackage{amsmath}')

#matplotlib.style.use('grayscale')


#PlotLinespecs = ['rD-', 'bo-', 'mv-.', 'cs-.', 'gx--', 'y+--', 'kp:', 'k*:'] #LHD113208t4640
#PlotLineColors = ['r', 'b', 'm', 'c', 'g', 'y', 'k', '0.75'] #LHD113208t4640
#PlotLinespecs = ['rD-', 'mv-.', 'gx--', 'y*--'] #W7-X NBI case
#PlotLineColors = ['r', 'm', 'g', '#AAAA00'] #W7-X NBI case
#PlotLinespecs = ['mD-', 'cv-.', 'gx--'] #TEMPORARY
#PlotLineColors = ['m', 'c', 'g'] #TEMPORARY
#PlotLinespecs = ['cv-.', 'bx--'] #TEMPORARY
#PlotLineColors = ['c', 'b'] #TEMPORARY
#PlotLinespecs = ['cv', 'b--', 'r-'] #TEMPORARY
#PlotLineColors = ['c', 'b', 'r'] #TEMPORARY
#PlotLinespecs = ['rD-', 'mv-.']
#PlotLineColors = ['r', 'm']
#PlotLinespecs = ['r-', 'b--', 'g-.', 'y:']
#PlotLineColors = ['r', 'b', 'g', '#AAAA00']
#PlotLinespecs = ['rD', 'g--', 'k-'] #EPS poster Er
#PlotLineColors = ['r', 'g', '#666666'] #EPS poster Er
#PlotLinespecs = ['r-', 'b--', 'g-.'] #EPS poster Collisionality
#PlotLineColors = ['#AA0000', '#0000AA', '#00AA00'] #EPS poster Collisionality
#PlotLinespecs = ['rD', 'bv', 'g*', 'y.'] #
#PlotLineColors = ['#AA0000', '#0000AA', '#00AA00', '#AAAA00'] #eta
#PlotLinespecs = ['b--', 'g-.'] #TTF poster TJ-II Collisionality
#PlotLineColors = ['#0000AA', '#00AA00'] #TTF poster TJ-II Collisionality
#PlotLinespecs = ['rD-', 'bo-', 'gx-.', 'cs-.'] #EPS poster Ar fluxes no Phi1
#PlotLineColors = ['r', 'b', 'g', 'c'] #EPS poster Ar fluxes no Phi1
#PlotLinespecs = ['rD-', 'mv-.', 'y*--'] #EPS poster Ar fluxes with Phi1
#PlotLineColors = ['r', 'm', 'y',] #EPS poster Ar fluxes with Phi1
#PlotLinespecs = ['rD-', 'mv-.', 'b*--', 'cs:'] #TTF poster C6+ fluxes with Phi1
#PlotLineColors = ['r', 'm', 'b', 'c'] #TTF poster C6+ fluxes with Phi1
#PlotLinespecs = ['rD', 'k-'] #EPS poster Ar fluxes vs XICS
#PlotLineColors = ['r', '#666666'] #EPS poster Ar fluxes vs XICS
#PlotLinespecs = ['gx', 'cs--'] #EPS poster GENE growth rates
#PlotLineColors = ['g', 'c'] #EPS poster GENE growth rates
#PlotLinespecs = ['cs--'] # W7X_180919.055 GENE growth rates
#PlotLineColors = ['c'] # W7X_180919.055 GENE growth rates
#PlotLinespecs = ['rD', 'mv', 'k-'] #EPS poster Ar fluxes vs XICS
#PlotLineColors = ['r', 'm', '#666666'] #EPS poster Ar fluxes vs XICS
#PlotLinespecs = ['k-']
#PlotLineColors = ['#555555']
#PlotLinespecs = ['rD-', 'bo-', 'mv-.', 'cs-.', 'gx--', 'y+--', 'kp:', 'k*:', 'k2--', 'k4--']
#PlotLineColors = ['r', 'b', 'm', 'c', 'g', 'y', 'k', '0.75', '0.5', '0.25']
#PlotLinespecs = ['rD-', 'bo-.', 'gx--', 'y+:', 'kp:', 'k*:', 'k2--', 'k4--']
#PlotLineColors = ['r', 'b', 'g', 'y', 'k', '0.75', '0.5', '0.25']
#PlotLinespecs = ['rD-', 'bo-', 'mv-.', 'cs-.', 'k-'] #W7X_180919.055/046 Ar fluxes D&V vs XICS
#PlotLineColors = ['r', 'b', 'm', 'c', '#666666'] #W7X_180919.055/046 Ar fluxes D&V vs XICS
#PlotLinespecs = ['rD-', 'k-'] #W7X_180919.055/046 D&V vs XICS only one
#PlotLineColors = ['r', '#666666'] #W7X_180919.055/046 D&V vs XICS only one

##NEW W7X_180919.055##
#Profiles#
#PlotLinespecs = ['r-', 'r-', 'm--'] #W7X_180919.055 Densities
#PlotLineColors = ['#b87333', '#b87333', '#36454f'] #W7X_180919.055 Densities
#PlotLinespecs = ['r-', 'm--', 'm--'] #W7X_180919.055 Temperatures
#PlotLineColors = ['#b87333', '#36454f', '#36454f'] #W7X_180919.055 Temperatures
#PlotLinespecs = ['r-', 'm--', 'c:'] #W7X_180919.055 Collisionalities
#PlotLineColors = ['#b87333', '#36454f', '#c0c0c0'] #W7X_180919.055 Collisionalities
#PlotLinespecs = ['k-'] #W7X_180919.055 Er
#PlotLineColors = ['k'] #W7X_180919.055 Er

#Comparison#
#PlotLinespecs = ['rD-', 'bo-'] #W7X_180919.055 Neo vs Cl no Phi1
#PlotLineColors = ['r', 'b'] #W7X_180919.055 Neo vs Cl no Phi1
#PlotLinespecs = ['rD-', 'mX--', 'gp:', 'y*-.'] #W7X_180919.055 Neo no Phi1 SFINCS vs EUTERPE vs DKES
#PlotLineColors = ['r', '#ffa6b1', 'g', 'y'] #W7X_180919.055 Neo no Phi1 SFINCS vs EUTERPE vs DKES
#PlotLinespecs = ['rD-', 'mv-.', 'mP:'] #W7X_180919.055 Neo with Phi1 SFINCS vs EUTERPE
#PlotLineColors = ['r', 'm', '#614051'] #W7X_180919.055 Neo with Phi1 SFINCS vs EUTERPE
#PlotLinespecs = ['rD-', 'bo-', 'm^-.', 'cs-.', 'k-'] #W7X_180919.055 D&V vs XICS
#PlotLineColors = ['r', 'b', '#967bb6', 'c', '#666666'] #W7X_180919.055 D&V vs XICS
PlotLinespecs = ['c--', 'b-.', 'r:', 'k-'] #W7X_180919.055 D&V vs XICS NEW
PlotLineColors = ['c', 'b', 'r', '#666666'] #W7X_180919.055 D&V vs XICS NEW

##Lavender #967bb6

#PlotMarkerEdgeWidth = [3, 3, 3, 0, 3, 3, 3, 3]
#PlotMarkerEdgeWidth = [3, 3, 3, 2]
#PlotMarkerEdgeWidth = [3, 3, 5, 2] #EPS poster Ar fluxes no Phi1
#PlotMarkerEdgeWidth = [5, 3] #EPS poster GENE growth rates
#PlotMarkerEdgeWidth = [3, 3, 5, 3] #TTF poster C6+ fluxes
#PlotMarkerEdgeWidth = [3, 3, 3, 3, 3]
PlotMarkerEdgeWidth = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
#PlotMarkerSize = 10
#PlotMarkerSize = 15
PlotMarkerSize = 20 # W7X_180919.055 Neo vs Cl #W7X_180919.055 Neo no/with Phi1 Phi1 SFINCS vs EUTERPE vs DKES #W7X_180919.055 profiles
#PlotMarkerSize = 10 #W7X_180919.055/046 D & V vs XICS
PlotLineWidth=8.0 # W7X_180919.055 Neo vs Cl #W7X_180919.055 Neo no/with Phi1 SFINCS vs EUTERPE vs DKES #W7X_180919.055 profiles
#PlotLineWidth=10.0 

xAxisScale = 'linear' #'linear', 'log' or 'symlog'
#yAxisScale = 'log' #'linear', 'log' or 'symlog'
#yAxisScale = 'symlog'
yAxisScale = 'linear'

ShowGrid = True

xAxisLabel = r'$r/a$'
#xAxisLabel = r'$r [\mathrm{m}]$'
#xAxisLabel = r'$E_r [\mathrm{kV} \mathrm{m}^{-1}]$'
#yAxisLabel = r'$n_{\mathrm{C}^{6+}}$ $[10^{19} \mathrm{m}^{-3}]$'
#yAxisLabel = r'$n_{s}$ $[10^{19} \mathrm{m}^{-3}]$'
#yAxisLabel = r'$n_{s}$'
#yAxisLabel = r'$\eta_{s} = d(\ln T_{s}) / d(\ln n_{s})$'
#yAxisLabel = r'$T_{s}$ $[\mathrm{keV}]$'
#yAxisLabel = r'$\nu_{ss}^{\prime}$'
#yAxisLabel = r'$\nu_{s}^{\prime}$'
#yAxisLabel = r'$\nu_{s}^{\prime} \equiv \frac{\left(G + \iota I\right)}{v_s B_{00}} \sum_{\alpha}  \nu_{s\alpha}$'
# e   i   \mathrm{He}^{2+}   \mathrm{C}^{6+}   \mathrm{Ne}^{10+}
#yAxisLabel = r'$<\mathbf{\Gamma}_{\mathrm{C}^{6+}} \cdot \nabla r> $ $[10^{20} \mathrm{m}^{-2} \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{z} \cdot \nabla r> $ $[10^{20} \mathrm{m}^{-2} \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{\mathrm{Ar}^{16+}} \cdot \nabla r> $ $[10^{20} \mathrm{m}^{-2} \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{\mathrm{Ar}^{16+}} \cdot \nabla r> / n_{\mathrm{Ar}^{16+}} $ $[\mathrm{m} \, \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{i} \cdot \nabla r> / n_{i} $ $[\mathrm{m} \, \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{e} \cdot \nabla r> / n_{e} $ $[\mathrm{m} \, \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{z} \cdot \nabla r> / n_{z} $ $[\mathrm{m} \, \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{\mathrm{C}^{6+}} \cdot \nabla r> / n_{\mathrm{C}^{6+}} $ $[\mathrm{m} \, \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{\mathrm{C}^{6+}}^{\mathrm{Classical}} \cdot \nabla r>\!/<\mathbf{\Gamma}_{\mathrm{C}^{6+}}^{\mathrm{Neoclassical}} \cdot \nabla r> $'
#yAxisLabel = r'$E_r$ $[\mathrm{kV/m}]$'
#yAxisLabel = r'$\omega_r, \gamma$ $[c_s / a]$'
#yAxisLabel = r'$\gamma$ $[c_s / a]$'
#yAxisLabel = r'$Z_{\mathrm{eff}}$'
#yAxisLabel = r'$e \Delta \Phi_1 / 2 T_i$'
yAxisLabel = [r'$D_{\mathrm{Ar}^{16+}}$ $[\mathrm{m}^2 \, \mathrm{s}^{-1}]$', r'$V_{\mathrm{Ar}^{16+}}$ $[\mathrm{m} \, \mathrm{s}^{-1}]$']
#yAxisLabel = [r'$D_{\mathrm{Ar}^{16+}}$ $[\mathrm{m}^2 \, \mathrm{s}^{-1}]$', r'$V_{\mathrm{Ar}^{16+}}$ $[\mathrm{m} \, \mathrm{s}^{-1}]$', r'$- d \ln n_{\mathrm{Ar}^{16+}} / dr$ $[\mathrm{m}^{-1}]$', r'$- D_{\mathrm{Ar}^{16+}} \cdot d \ln n_{\mathrm{Ar}^{16+}} / dr$ $[\mathrm{m} \, \mathrm{s}^{-1}]$']
AxesLabelSize = 50
#AxesLabelSize = 40.4 ##TTF poster classical over neoclassical fluxes

TickSize = 45

AxisLimAuto = True
#AxisLimAuto = False
#xAxisLim = [0.18, 0.92] #LHD113208t4640 Er scan
#yAxisLim = [-4.5, 0.0]
#yAxisLim = [-3.2, 0.0] #LHD113208t4640 Er scan
#yAxisLim = [-15.0, 15.0] #LHD113208t4640 eta
#xAxisLim = [0.10, 1.02] #W7-X NBI case Er scan
#yAxisLim = [-25.0, 0.0] #W7-X NBI case Er scan
#xAxisLim = [0.0, 1.02] #EPS poster Er
#yAxisLim = [-30.0, 30.0] #EPS poster Er
#xAxisLim = [0.0, 0.60] #EPS poster Ar fluxes vs XICS
#yAxisLim = [-2.0, 6.0] #EPS poster Ar fluxes vs XICS
#yAxisLim = [0.0, 0.4] #TEMPORARY
#xAxisLim = [0.0, 1.0]
#yAxisLim = [0.00007, 0.15]
#yAxisLim = [0.00004, 0.40]
#yAxisLim = [0.003, 3.0]
#yAxisLim = [0.01, 5.0]
#yAxisLim = [0.0, 1.5]
#yAxisLim = [0.0, 5.5]
#yAxisLim = [1.0, 2.5]
#xAxisLim = [0.082, 0.953]  # W7X_180919.055 delta Phi1
#yAxisLim = [0.0, 0.005]  # W7X_180919.055 delta Phi1
#xAxisLim = [0.0, 1.0] # W7X_180919.055 GENE growth rates
#yAxisLim = [0.0, 0.6] # W7X_180919.055 GENE growth rates
xAxisLim = [-0.05005, 1.0314]  # W7X_180919.046 Ar fluxes vs XICS
yAxisLim = [-2.2, 4.0]  # W7X_180919.046 Ar fluxes vs XICS
#xAxisLim = [[-0.043, 0.80], [-0.043, 0.80]]  # W7X_180919.046 Ar fluxes and D & V vs XICS
#yAxisLim = [[-0.7, 1.3], [-4.0, 3.0]] # W7X_180919.046 D & V vs XICS
yAxisLim = [-0.1, 1.0]

xAxisLabelCoords = [0.5,-0.09]
#yAxisLabelCoords = [-0.11,0.5] #LHD113208t4640, W7-X NBI case, Er scan, Flux scan 
#yAxisLabelCoords = [-0.105,0.5] #EPS poster Er, TTF poster TJ-II fluxes # W7X_180919.055 GENE growth rates
#yAxisLabelCoords = [-0.09,0.5] #TTF poster TJ-II Classical over neoclassical fluxes, collisionality #W7X_180919.055 Ar fluxes vs XICS
yAxisLabelCoords = [-0.09,0.5] #CONTOUR PLOT WITHOUT COLORBAR, #EPS poster Collisionality, fluxes
#yAxisLabelCoords = [-0.15,0.5]
#yAxisLabelCoords = [-0.12,0.5] # W7X_180919.055 Neo vs Cl #W7X_180919.055 Neo no/with Phi1 SFINCS vs EUTERPE vs DKES #W7X_180919.055 profiles except collisionality
#yAxisLabelCoords = [-0.14,0.5] #W7X_180919.055 collisionality


#LeftMargin = 0.15 #LHD113208t4640, W7-X NBI case, Er scan, Flux scan #EPS poster Er, fluxes, TTF poster TJ-II Classical over neoclassical fluxes
#LeftMargin = 0.17 #EPS poster Collisionality, TTF poster TJ-II collisionality # W7X_180919.055 Neo vs Cl #W7X_180919.055 Neo no/with Phi1 SFINCS vs EUTERPE vs DKES #W7X_180919.055 profiles except collisionality
#LeftMargin = 0.135 #CONTOUR PLOT WITHOUT COLORBAR
LeftMargin = 0.145 #CONTOUR PLOT WITHOUT COLORBAR WITH TITLE
#LeftMargin = 0.205 #W7X_180919.055 collisionality
#LeftMargin = 0.10 #W7X_180919.055/046 D & V vs XICS

#RightMargin = 0.95 #LHD113208t4640, W7-X NBI case, Er scan, Flux scan #EPS poster Er, fluxes, TTF poster TJ-II Classical over neoclassical fluxes #W7X_180919.055 Neo no/with Phi1 SFINCS vs EUTERPE vs DKES #W7X_180919.055 profiles
#RightMargin = 0.97 #EPS poster Collisionality, TTF poster TJ-II collisionality
#RightMargin = 0.935 #CONTOUR PLOT WITHOUT COLORBAR
RightMargin = 0.925 #CONTOUR PLOT WITHOUT COLORBAR WITH TITLE
#RightMargin = 0.98 #W7X_180919.055/046 D & V vs XICS

#TopMargin = 0.95 #LHD113208t4640, W7-X NBI case, Er scan, Flux scan #EPS poster Er, fluxes, TTF poster TJ-II
#TopMargin = 0.97 #CONTOUR PLOT WITHOUT COLORBAR
TopMargin = 0.95 #CONTOUR PLOT WITHOUT COLORBAR WITH TITLE
#TopMargin = 0.935 # W7X_180919.055 Neo vs Cl #W7X_180919.055 Neo no/with Phi1 SFINCS vs EUTERPE vs DKES #W7X_180919.055 profiles

#BottomMargin = 0.15 #LHD113208t4640, W7-X NBI case, Er scan, Flux scan #EPS poster Er, fluxes, TTF poster TJ-II
#BottomMargin = 0.145 # W7X_180919.055 Neo vs Cl #W7X_180919.055 Neo no/with Phi1 SFINCS vs EUTERPE vs DKES #W7X_180919.055 profiles
BottomMargin = 0.17 #CONTOUR PLOT WITHOUT COLORBAR

ShowLegend = False

#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Pitch-angle scattering w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$', r'SFINCS Pitch-angle scattering w/ $\Phi_1$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Pitch-angle scattering w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$', r'SFINCS Pitch-angle scattering w/ $\Phi_1$', r'DKES + momentum correction', r'DKES (no momentum correction)']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$', r'DKES + momentum correction', r'DKES (no momentum correction)']
#PlotLegendLabels = [r'$n_e$', r'$n_i$', r'$n_{\mathrm{He}^{2+}}$', r'$n_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$n_e$~$[10^{19} \mathrm{m}^{-3}]$', r'$n_i$~$[10^{19} \mathrm{m}^{-3}]$', r'$n_{\mathrm{Ar}^{16+}}$~$[\mathrm{AU}]$']
#PlotLegendLabels = [r'$n_e$', r'$n_i$', r'$n_{\mathrm{He}^{2+}}$', r'$10 \times n_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$T_e$', r'$T_i$', r'$T_{\mathrm{He}^{2+}}$', r'$T_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$\eta_{e}$', r'$\eta_{i}$', r'$\eta_{\mathrm{He}^{2+}}$', r'$\eta_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$T_e$', r'$T_i$', r'$T_{\mathrm{Ne}^{10+}}$']
#PlotLegendLabels = [r'$T_e$', r'$T_i$', r'$T_{\mathrm{Ar}^{16+}}$']
#PlotLegendLabels = [r'$T_e$', r'$T_i = T_{\mathrm{He}^{2+}} = T_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$\nu_{ee}^{\prime}$', r'$\nu_{ii}^{\prime}$', r'$\nu_{HeHe}^{\prime}$', r'$\nu_{CC}^{\prime}$']
#PlotLegendLabels = [r'$\nu_{e}^{\prime}$', r'$\nu_{i}^{\prime}$', r'$\nu_{He}^{\prime}$', r'$\nu_{C}^{\prime}$']
#PlotLegendLabels = [r'$\nu_{ee}^{\prime}$', r'$\nu_{ii}^{\prime}$', r'$\nu_{NeNe}^{\prime}$']
#PlotLegendLabels = [r'$\nu_{e}^{\prime}$', r'$\nu_{i}^{\prime}$', r'$\nu_{Ne}^{\prime}$']
#PlotLegendLabels = [r'$\nu_{e}^{\prime}$', r'$\nu_{i}^{\prime}$', r'$\nu_{Ar}^{\prime}$'] #EPS poster Collisionality
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/ $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$ + magnetic drifts']
#PlotLegendLabels = [r'LHD113208t4640 Fokker-Planck w/o $\Phi_1$', r'Inward-shifted Fokker-Planck w/o $\Phi_1$', r'LHD113208t4640 Fokker-Planck w/ $\Phi_1$', r'Inward-shifted Fokker-Planck w/ $\Phi_1$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'DKES', r'XICS']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'XICS']
#PlotLegendLabels = [r'$- \omega_r$', r'$\gamma$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$',  r'SFINCS Fokker-Planck w/ $\Phi_1$', r'XICS']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Pitch-angle scattering w/o $\Phi_1$']

#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Pitch-angle scattering w/o $\Phi_1$', r'EUTERPE mom. conservation w/o $\Phi_1$', r'EUTERPE Pitch-angle scattering w/o $\Phi_1$'] #EPS poster Ar fluxes no Phi1
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$'] #EPS poster Ar fluxes with Phi1
#PlotLegendLabels = [r'SFINCS Pitch-angle scattering w/o $\Phi_1$', r'SFINCS Pitch-angle scattering w/ $\Phi_1$', r'EUTERPE Pitch-angle scattering w/o $\Phi_1$', r'EUTERPE Pitch-angle scattering w/ $\Phi_1$'] ##TTF poster C6+ fluxes
#PlotLegendLabels = [r'$\nu_{i}^{\prime}$', r'$\nu_{C}^{\prime}$'] #TTF poster TJ-II Collisionality
#PlotLegendLabels = [r'SFINCS w/o $\Phi_1$', r'SFINCS w/ $\Phi_1$', r'SFINCS w/ $\Phi_1$ EUTERPE QN'] #EPS poster Ar fluxes with Phi with EUTERPE QN
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Classical w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$', r'SFINCS Classical w/ $\Phi_1$', r'DKES Neoclassical', r'DKES Classical']
#PlotLegendLabels = [r'SFINCS Neoclassical w/o $\Phi_1$', r'SFINCS Classical w/o $\Phi_1$', r'SFINCS Neoclassical + Classical w/o $\Phi_1$', r'SFINCS Neoclassical w/ $\Phi_1$', r'SFINCS Classical w/ $\Phi_1$', r'SFINCS Neoclassical + Classical w/ $\Phi_1$', r'DKES Neoclassical']
#PlotLegendLabels = [r'SFINCS Neoclassical w/o $\Phi_1$', r'SFINCS Neoclassical + Classical w/o $\Phi_1$', r'SFINCS Neoclassical PAS w/o $\Phi_1$', r'SFINCS Neoclassical + Classical PAS w/o $\Phi_1$', r'SFINCS Neoclassical w/ $\Phi_1$', r'SFINCS Neoclassical + Classical w/ $\Phi_1$', r'DKES Neoclassical', r'DKES Neoclassical + Classical']
#PlotLegendLabels = [r'SFINCS Neoclassical w/o $\Phi_1$', r'SFINCS Neoclassical + Classical w/o $\Phi_1$', r'SFINCS Neoclassical PAS w/o $\Phi_1$', r'SFINCS Neoclassical + Classical PAS w/o $\Phi_1$', r'DKES Neoclassical', r'Classical + DKES Neoclassical']
#PlotLegendLabels = [r'SFINCS Neoclassical w/o $\Phi_1$', r'SFINCS Neoclassical + Classical w/o $\Phi_1$', r'DKES Neoclassical with mom. corr.', r'Classical + DKES Neoclassical with mom. corr.']
#PlotLegendLabels = [r'SFINCS Neoclassical w/o $\Phi_1$', r'SFINCS Neoclassical + Classical w/o $\Phi_1$', r'SFINCS Neoclassical w/ $\Phi_1$', r'SFINCS Neoclassical + Classical w/ $\Phi_1$', r'DKES Neoclassical with mom. corr.', r'Classical + DKES Neoclassical with mom. corr.']
#PlotLegendLabels = [r'Neoclassical, real $n_z$', r'Neoclassical + Classical, real $n_z$', r'Neoclassical, const. $Z_{\mathrm{eff}}$', r'Neoclassical + Classical, const. $Z_{\mathrm{eff}}$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck, kinetic $e^{-}$', r'SFINCS Pitch-angle scattering, kinetic $e^{-}$', r'KNOSOS Pitch-angle scattering, kinetic $e^{-}$', r'KNOSOS Pitch-angle scattering, adiabatic $e^{-}$', r'EUTERPE Pitch-angle scattering + mom. corr., adiabatic $e^{-}$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck, kinetic $e^{-}$', r'SFINCS Pitch-angle scattering, kinetic $e^{-}$', r'KNOSOS Pitch-angle scattering, kinetic $e^{-}$', r'KNOSOS Pitch-angle scattering, kinetic $e^{-}$, $\mathbf{B}$-drift', r'EUTERPE Pitch-angle scattering + mom. corr., adiabatic $e^{-}$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck, kinetic $e^{-}$', r'KNOSOS Pitch-angle scattering, kinetic $e^{-}$', r'EUTERPE Pitch-angle scattering + mom. corr., adiabatic $e^{-}$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$', r'EUTERPE Pitch-angle scattering + mom. corr. w/ $\Phi_1$', r'DKES Pitch-angle scattering w/o $\Phi_1$']
#PlotLegendLabels = [r'$10 \times$ SFINCS w/o $\Phi_1$ Neoclassical only', r'$10 \times$ SFINCS w/o $\Phi_1$',  r'$10 \times$ SFINCS w/ $\Phi_1$', r'SFINCS w/ $\Phi_1$', r'XICS']
#PlotLegendLabels = [r'$10 \times$ SFINCS w/o $\Phi_1$ Neocl. only', r'$10 \times$ SFINCS w/o $\Phi_1$',  r'$10 \times$ SFINCS w/ $\Phi_1$', r'SFINCS w/ $\Phi_1$', r'XICS']
#PlotLegendLabels = [r'Standard version', r'Fourier version']
#PlotLegendLabels = [r'SFINCS', r'XICS']
#PlotLegendLabels = [r'SFINCS', r'EUTERPE', r'DKES', r'DKES + momentum corr.']
#PlotLegendLabels = [r'SFINCS', r'DKES', r'DKES + momentum corr.']
#PlotLegendLabels = [r'$\mathrm{Si}^{12+}$', r'$\mathrm{W}^{44+}$']
#PlotLegendLabels = [r'$\mathrm{Si}^{12+}$ Neoclassical', r'$\mathrm{Si}^{12+}$ Neoclassical + Classical', r'$\mathrm{W}^{44+}$ Neoclassical', r'$\mathrm{W}^{44+}$ Neoclassical + Classical']

##NEW W7X_180919.055##
#Profiles#
#PlotLegendLabels = [r'$n_e$~$[10^{19} \mathrm{m}^{-3}]$', r'$n_i$~$[10^{19} \mathrm{m}^{-3}]$', r'$n_{\mathrm{Ar}^{16+}}$~$[\mathrm{AU}]$']
#PlotLegendLabels = [r'$T_e$', r'$T_i$', r'$T_{\mathrm{Ar}^{16+}}$']
#PlotLegendLabels = [r'$\nu_{e}^{\prime}$', r'$\nu_{i}^{\prime}$', r'$\nu_{Ar}^{\prime}$']

#Comparison#
#PlotLegendLabels = [r'SFINCS w/o $\Phi_1$ Neoclassical only', r'SFINCS w/o $\Phi_1$ Classical + Neoclassical'] #W7X_180919.055 Neo vs Cl no Phi1
#PlotLegendLabels = [r'SFINCS w/o $\Phi_1$', r'EUTERPE w/o $\Phi_1$', r'DKES no momentum corr.', r'DKES + momentum corr.'] #W7X_180919.055 Neo no Phi1 SFINCS vs EUTERPE vs DKES
#PlotLegendLabels = [r'SFINCS w/o $\Phi_1$', r'SFINCS w/ $\Phi_1$', r'EUTERPE w/ $\Phi_1$'] #W7X_180919.055 Neo with Phi1 SFINCS vs EUTERPE vs DKES
#PlotLegendLabels = [r'$10 \times$ SFINCS w/o $\Phi_1$ Neocl. only', r'$10 \times$ SFINCS w/o $\Phi_1$',  r'$10 \times$ SFINCS w/ $\Phi_1$', r'SFINCS w/ $\Phi_1$', r'XICS']
PlotLegendLabels = [r'SFINCS', r'$\mathbf{10} \mathbf{\times}$ SFINCS', r'$\mathbf{10} \mathbf{\times}$ SFINCS Neocl. only', r'XICS']

#PlotLegendLabels = []

#LegendFontSize = 15
#LegendProperties = {'weight':'bold'}
#LegendProperties ={'size':'20', 'weight':'heavy'}
#LegendProperties ={'size':'15', 'weight':'heavy'}
#LegendProperties ={'size':'30', 'weight':'heavy'} #LHD113208t4640, W7-X NBI case, Er scan, #EPS poster fluxes no Phi1 XICS #W7X_180919.055 Neo vs Cl no Phi1 # Neo no/with Phi1 SFINCS vs EUTERPE vs DKES #W7X_180919.055 profiles
LegendProperties ={'size':'25', 'weight':'heavy'} #LHD113208t4640 Helium scan, #EPS poster fluxes no Phi1
#LegendProperties ={'size':'27', 'weight':'heavy'} #EPS poster Er
#LegendProperties ={'size':'35', 'weight':'heavy'} #EPS poster Collisionality, fluxes with Phi1, TTF poster TJ-II collisionality
#LegendProperties ={'size':'30.5', 'weight':'heavy'} #TTF poster TJ-II classical over neoclassical
#LegendProperties ={'size':'22', 'weight':'heavy'} #W7X_180919.055 normalized Ar fluxes vs XICS
#LegendProperties ={'size':'15', 'weight':'heavy'}
LegendPosition = 3
LegendNumberColumns = 1

#LegendBBoxToAnchor = [0.01, 0.01, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.765, 1., .102]
#LegendBBoxToAnchor = [0.435, 0.765, 1., .102]
#LegendBBoxToAnchor = [0.30, 0.01, 1., .102]
#LegendBBoxToAnchor = [0.13, 0.01, 1., .102]
#LegendBBoxToAnchor = [0.01, 0.655, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.592, 1., .102] #LHD113208t4640 Er scan
#LegendBBoxToAnchor = [0.005, 0.35, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.335, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.375, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.47, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.005, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.58, 1., .102]
#LegendBBoxToAnchor = [0.777, 0.005, 1., .102]
#LegendBBoxToAnchor = [0.32, 0.005, 1., .102]
#LegendBBoxToAnchor = [0.165, 0.005, 1., .102] #W7-X NBI case Electron scan
#LegendBBoxToAnchor = [0.005, 0.005, 1., .102] #EPS poster Er, TTF poster TJ-II fluxes
#LegendBBoxToAnchor = [0.01, 0.43, 1., .102] #EPS poster Collisionality, TTF poster TJ-II collisionality
#LegendBBoxToAnchor = [0.28, 0.005, 1., .102] #EPS poster fluxes no Phi1
#LegendBBoxToAnchor = [0.20, 0.005, 1., .102] #EPS poster fluxes with Phi1
#LegendBBoxToAnchor = [0.005, 0.828, 1., .102] #EPS poster fluxes no Phi1 XICS
#LegendBBoxToAnchor = [0.005, 0.805, 1., .102] #EPS poster GENE growth rates
#LegendBBoxToAnchor = [0.005, 0.85, 1., .102]
#LegendBBoxToAnchor = [0.168, 0.005, 1., .102] ##TTF poster TJ-II classical over neoclassical fluxes PAS
#LegendBBoxToAnchor = [0.30, 0.005, 1., .102] ##TTF poster TJ-II classical over neoclassical fluxes FP
#LegendBBoxToAnchor = [0.55, 0.005, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.55, 1., .102]
#LegendBBoxToAnchor = [0.23, 0.88, 1., .102]
#LegendBBoxToAnchor = [0.2, 0.8, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.68, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.83, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.805, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.76, 1., .102]
#LegendBBoxToAnchor = [0.315, 0.55, 1., .102]
#LegendBBoxToAnchor = [0.51, 0.68, 1., .102]
#LegendBBoxToAnchor = [0.49, 0.005, 1., .102]
#LegendBBoxToAnchor = [0.32, 0.67, 1., .102]
#LegendBBoxToAnchor = [0.01, 0.69, 1., .102] # W7X_180919.055/046 collisionality
#LegendBBoxToAnchor = [0.56, 0.005, 1., .102]  # W7X_180919.055 densities
#LegendBBoxToAnchor = [0.005, 0.005, 1., .102]  # W7X_180919.046 densities # W7X_180919.055/046 temperatures # W7X_180919.055 delta Phi1 # W7X_180919.055 Neo vs Cl
#LegendBBoxToAnchor = [0.19, 0.75, 1., .102]  # W7X_180919.055 comparison normalized neoclassical fluxes
#LegendBBoxToAnchor = [0.005, 0.70, 1., .102]  #W7X_180919.055/046 normalized Ar fluxes vs XICS, W7X_180919.055/046 D & V
#LegendBBoxToAnchor = [0.005, 0.85, 1., .102]  #W7X_180919.055/046 normalized Ar fluxes vs XICS, W7X_180919.055/046 D & V extended
#LegendBBoxToAnchor = [0.005, 0.67, 1., .102] #W7X_180919.055 Neo no Phi1 SFINCS vs EUTERPE vs DKES
#LegendBBoxToAnchor = [0.365, 0.67, 1., .102] #W7X_180919.055 Neo no Phi1 SFINCS vs EUTERPE vs DKES main ions
#LegendBBoxToAnchor = [0.005, 0.005, 1., .102] #W7X_180919.055 Neo with Phi1 SFINCS vs EUTERPE #W7X_180919.055 temperatures
#LegendBBoxToAnchor = [0.525, 0.005, 1., .102] #W7X_180919.055 Neo with Phi1 SFINCS vs EUTERPE vs DKES main ions
#LegendBBoxToAnchor = [0.615, 0.005, 1., .102] #W7X_180919.055 densities
#LegendBBoxToAnchor = [0.765, 0.005, 1., .102] #W7X_180919.055 collisionality
LegendBBoxToAnchor = [0.355, 0.005, 1., .102]

ShowSubPlotLabel = False
SubPlotLabel = '(e)'
SubPlotLabelSize = 40
#SubPlotLabelXcoord = 0.085 #W7-X NBI case Electron scan, Ion scan, Neon scan
#SubPlotLabelXcoord = 0.40
#SubPlotLabelXcoord = 0.01
#SubPlotLabelXcoord = 0.003 #CONTOUR PLOT WITHOUT COLORBAR LHD discharge 113208 at t = 4.64 s
SubPlotLabelXcoord = 0.006 #CONTOUR PLOT WITHOUT COLORBAR W7-X_NBI_case_Q34Q78_Z10_Zeff2p0 #W7X_180919.055 CONTOUR PLOT WITHOUT COLORBAR
#SubPlotLabelXcoord = 0.28 #LHD113208t4640 Carbon scan
#SubPlotLabelXcoord = 0.198 #LHD113208t4640 Electron scan, Ion scan, Helium scan
#SubPlotLabelXcoord = 0.91 #W7X_180919.055 Neo no Phi1 SFINCS vs EUTERPE vs DKES
#SubPlotLabelXcoord = 0.007 #W7X_180919.055 Neo no Phi1 SFINCS vs EUTERPE vs DKES main ion # W7X_180919.055 Neo vs Cl #W7X_180919.055 profiles

#SubPlotLabelYcoord = -0.00071 #W7-X NBI case Neon scan
#SubPlotLabelYcoord = 0.1302 #W7-X NBI case Ion scan
#SubPlotLabelYcoord = 0.08
#SubPlotLabelYcoord = 0.0945 #W7-X NBI case Electron scan
#SubPlotLabelYcoord = 2.090
#SubPlotLabelYcoord = 3.6
#SubPlotLabelYcoord = 1.422
#SubPlotLabelYcoord = 1.222 #LHD113208t4640 Electron scan
#SubPlotLabelYcoord = 1.14 #LHD113208t4640 Ion scan
#SubPlotLabelYcoord = 0.0746 #LHD113208t4640 Helium scan
#SubPlotLabelYcoord = 0.001 #LHD113208t4640 Carbon scan
#SubPlotLabelYcoord = 5.2
#SubPlotLabelYcoord = 5.84 #CONTOUR PLOT WITHOUT COLORBAR LHD discharge 113208 at t = 4.64 s
#SubPlotLabelYcoord = 5.715 #CONTOUR PLOT WITHOUT COLORBAR W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
#SubPlotLabelYcoord = 5.55 #CONTOUR PLOT WITHOUT COLORBAR W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
#SubPlotLabelYcoord = 2.422
#SubPlotLabelYcoord = 0.1
#SubPlotLabelYcoord = 0.248
#SubPlotLabelYcoord = 0.4 #W7X_180919.055 Neo no Phi1 SFINCS vs EUTERPE vs DKES
#SubPlotLabelYcoord = 0.225 #W7X_180919.055 Neo no Phi1 SFINCS vs EUTERPE vs DKES main ion
#SubPlotLabelYcoord = 0.415 # W7X_180919.055 Neo vs Cl
#SubPlotLabelYcoord = 0.226 # W7X_180919.055 Neo vs Cl main ion
#SubPlotLabelYcoord = 5.75 #W7X_180919.055 densities
#SubPlotLabelYcoord = 3.55 #W7X_180919.055 temperatures
#SubPlotLabelYcoord = 300 #W7X_180919.055 collisionality
#SubPlotLabelYcoord = 3.4 #W7X_180919.055 Er
#SubPlotLabelYcoord = 5.78 #W7X_180919.055 CONTOUR PLOT WITHOUT COLORBAR SFINCS
#SubPlotLabelYcoord = 5.97 #W7X_180919.055 CONTOUR PLOT WITHOUT COLORBAR EUTERPE
SubPlotLabelYcoord = 5.873 #W7X_180919.055 CONTOUR PLOT WITHOUT COLORBAR KNOSOS

FilledErrors = False

ErrorBars = [False, False, False, True]
ErrorBarAlpha = 0.3

ShowLineAtXzero = True
ShowLineAtYzero = False

NoScientificAxes = False

ChangeXaxisTicks = False
#NewXaxisTicks = np.arange(0.0, 1.1, 0.2)

ChangeYaxisTicks = False
#NewYaxisTicks = np.arange(-10.0, 4.0, 2.0)

#For log W7X_180919.055 collisionality#
#NewYaxisTicks = np.append(np.arange(0.001, 0.01, 0.001), np.arange(0.01, 0.1, 0.01))
#NewYaxisTicks = np.append( NewYaxisTicks, np.arange(0.1, 1.0, 0.1))
#NewYaxisTicks = np.append( NewYaxisTicks, np.arange(1.0, 10.0, 1.0))
#NewYaxisTicks = np.append( NewYaxisTicks, np.arange(10.0, 100.0, 10.0))
#NewYaxisTicks = np.append( NewYaxisTicks, np.arange(100.0, 1000.0, 100.0))
##
