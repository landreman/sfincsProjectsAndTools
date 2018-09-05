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

#FigSize = (12,10) #LHD113208t4640, W7-X NBI case, Er scan, Flux scan, #EPS poster Er, Collisionality, Ar fluxes
#FigSize = (13.5,10) #Phi1 plot
#FigSize = (13.0,10.5) #Phi1 plot
FigSize = (12,12) #Phi1 plot without colorbar
font = {'size':35} #Phi1 plot without colorbar
#font = {'size':38} #EPS poster
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
PlotLinespecs = ['rD', 'bv', 'g*', 'y.'] #
PlotLineColors = ['#AA0000', '#0000AA', '#00AA00', '#AAAA00'] #eta
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
#PlotLinespecs = ['rD', 'mv', 'k-'] #EPS poster Ar fluxes vs XICS
#PlotLineColors = ['r', 'm', '#666666'] #EPS poster Ar fluxes vs XICS
#PlotLinespecs = ['k-']
#PlotLineColors = ['#555555']
#PlotMarkerEdgeWidth = [3, 3, 3, 0, 3, 3, 3, 3]
#PlotMarkerEdgeWidth = [3, 3, 3, 2]
#PlotMarkerEdgeWidth = [3, 3, 5, 2] #EPS poster Ar fluxes no Phi1
#PlotMarkerEdgeWidth = [5, 3] #EPS poster GENE growth rates
PlotMarkerEdgeWidth = [3, 3, 5, 3] #TTF poster C6+ fluxes
#PlotMarkerSize = 10
PlotMarkerSize = 15
PlotLineWidth=7.0

xAxisScale = 'linear' #'linear', 'log' or 'symlog'
#yAxisScale = 'log' #'linear', 'log' or 'symlog'
#yAxisScale = 'symlog'
yAxisScale = 'linear'

ShowGrid = True

xAxisLabel = r'$r/a$'
#yAxisLabel = r'$n_{\mathrm{C}^{6+}}$ $[10^{19} \mathrm{m}^{-3}]$'
#yAxisLabel = r'$n$ $[10^{19} \mathrm{m}^{-3}]$'
yAxisLabel = r'$\eta_{s} = d(\ln T_{s}) / d(\ln n_{s})$'
#yAxisLabel = r'$T$ $[\mathrm{keV}]$'
#yAxisLabel = r'$\nu_{ss}^{\prime}$'
#yAxisLabel = r'$\nu_{s}^{\prime}$'
#yAxisLabel = r'$\nu_{s}^{\prime} \equiv \frac{\left(G + \iota I\right)}{v_s B_{00}} \sum_{\alpha}  \nu_{s\alpha}$'
# e   i   \mathrm{He}^{2+}   \mathrm{C}^{6+}   \mathrm{Ne}^{10+}
#yAxisLabel = r'$<\mathbf{\Gamma}_{\mathrm{C}^{6+}} \cdot \nabla r> $ $[10^{20} \mathrm{m}^{-2} \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{\mathrm{Ar}^{16+}} \cdot \nabla r> $ $[10^{20} \mathrm{m}^{-2} \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{\mathrm{Ar}^{16+}} \cdot \nabla r> / n_{\mathrm{Ar}^{16+}} $ $[\mathrm{m} \, \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{\mathrm{C}^{6+}} \cdot \nabla r> / n_{\mathrm{C}^{6+}} $ $[\mathrm{m} \, \mathrm{s}^{-1}]$'
#yAxisLabel = r'$<\mathbf{\Gamma}_{\mathrm{C}^{6+}}^{\mathrm{Classical}} \cdot \nabla r>\!/<\mathbf{\Gamma}_{\mathrm{C}^{6+}}^{\mathrm{Neoclassical}} \cdot \nabla r> $'
#yAxisLabel = r'$E_r$ $[\mathrm{kV/m}]$'
#yAxisLabel = r'$\omega_r, \gamma$ $[c_s / a]$'
#yAxisLabel = r'$Z_{\mathrm{eff}}$'
#AxesLabelSize = 50
AxesLabelSize = 43 ##TTF poster classical over neoclassical fluxes

TickSize = 45

#AxisLimAuto = True
AxisLimAuto = False
xAxisLim = [0.18, 0.92] #LHD113208t4640 Er scan
#yAxisLim = [-4.5, 0.0]
#yAxisLim = [-3.2, 0.0] #LHD113208t4640 Er scan
yAxisLim = [-15.0, 15.0] #LHD113208t4640 eta
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

xAxisLabelCoords = [0.5,-0.09]
#yAxisLabelCoords = [-0.11,0.5] #LHD113208t4640, W7-X NBI case, Er scan, Flux scan
#yAxisLabelCoords = [-0.105,0.5] #EPS poster Er, TTF poster TJ-II fluxes
#yAxisLabelCoords = [-0.09,0.5] #TTF poster TJ-II Classical over neoclassical fluxes, collisionality
yAxisLabelCoords = [-0.09,0.5] #CONTOUR PLOT WITHOUT COLORBAR, #EPS poster Collisionality, fluxes
#yAxisLabelCoords = [-0.15,0.5]

#LeftMargin = 0.15 #LHD113208t4640, W7-X NBI case, Er scan, Flux scan #EPS poster Er, fluxes, TTF poster TJ-II Classical over neoclassical fluxes
#LeftMargin = 0.17 #EPS poster Collisionality, TTF poster TJ-II collisionality
#LeftMargin = 0.135 #CONTOUR PLOT WITHOUT COLORBAR
LeftMargin = 0.145 #CONTOUR PLOT WITHOUT COLORBAR WITH TITLE
#LeftMargin = 0.25

#RightMargin = 0.95 #LHD113208t4640, W7-X NBI case, Er scan, Flux scan #EPS poster Er, fluxes, TTF poster TJ-II Classical over neoclassical fluxes
#RightMargin = 0.97 #EPS poster Collisionality, TTF poster TJ-II collisionality
#RightMargin = 0.935 #CONTOUR PLOT WITHOUT COLORBAR
RightMargin = 0.925 #CONTOUR PLOT WITHOUT COLORBAR WITH TITLE

#TopMargin = 0.95 #LHD113208t4640, W7-X NBI case, Er scan, Flux scan #EPS poster Er, fluxes, TTF poster TJ-II
#TopMargin = 0.97 #CONTOUR PLOT WITHOUT COLORBAR
TopMargin = 0.95 #CONTOUR PLOT WITHOUT COLORBAR WITH TITLE

#BottomMargin = 0.15 #LHD113208t4640, W7-X NBI case, Er scan, Flux scan #EPS poster Er, fluxes, TTF poster TJ-II
BottomMargin = 0.17 #CONTOUR PLOT WITHOUT COLORBAR

ShowLegend = True

#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Pitch-angle scattering w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$', r'SFINCS Pitch-angle scattering w/ $\Phi_1$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Pitch-angle scattering w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$', r'SFINCS Pitch-angle scattering w/ $\Phi_1$', r'DKES + momentum correction', r'DKES (no momentum correction)']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$', r'DKES + momentum correction', r'DKES (no momentum correction)']
#PlotLegendLabels = [r'$n_e$', r'$n_i$', r'$n_{\mathrm{He}^{2+}}$', r'$n_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$n_e$', r'$n_i$', r'$n_{\mathrm{He}^{2+}}$', r'$10 \times n_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$T_e$', r'$T_i$', r'$T_{\mathrm{He}^{2+}}$', r'$T_{\mathrm{C}^{6+}}$']
PlotLegendLabels = [r'$\eta_{e}$', r'$\eta_{i}$', r'$\eta_{\mathrm{He}^{2+}}$', r'$\eta_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$T_e$', r'$T_i$', r'$T_{\mathrm{Ne}^{10+}}$']
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

#PlotLegendLabels = []

#LegendFontSize = 15
#LegendProperties = {'weight':'bold'}
#LegendProperties ={'size':'20', 'weight':'heavy'}
LegendProperties ={'size':'30', 'weight':'heavy'} #LHD113208t4640, W7-X NBI case, Er scan, #EPS poster fluxes no Phi1 XICS
#LegendProperties ={'size':'25', 'weight':'heavy'} #LHD113208t4640 Helium scan, #EPS poster fluxes no Phi1
#LegendProperties ={'size':'27', 'weight':'heavy'} #EPS poster Er
#LegendProperties ={'size':'35', 'weight':'heavy'} #EPS poster Collisionality, fluxes with Phi1, TTF poster TJ-II collisionality
#LegendProperties ={'size':'30.5', 'weight':'heavy'} #TTF poster TJ-II classical over neoclassical
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
LegendBBoxToAnchor = [0.005, 0.005, 1., .102] #EPS poster Er, TTF poster TJ-II fluxes
#LegendBBoxToAnchor = [0.01, 0.43, 1., .102] #EPS poster Collisionality, TTF poster TJ-II collisionality
#LegendBBoxToAnchor = [0.28, 0.005, 1., .102] #EPS poster fluxes no Phi1
#LegendBBoxToAnchor = [0.20, 0.005, 1., .102] #EPS poster fluxes with Phi1
#LegendBBoxToAnchor = [0.005, 0.828, 1., .102] #EPS poster fluxes no Phi1 XICS
#LegendBBoxToAnchor = [0.005, 0.805, 1., .102] #EPS poster GENE growth rates
#LegendBBoxToAnchor = [0.005, 0.85, 1., .102]
#LegendBBoxToAnchor = [0.168, 0.005, 1., .102] ##TTF poster TJ-II classical over neoclassical fluxes

ShowSubPlotLabel = False
SubPlotLabel = '(c)'
SubPlotLabelSize = 40
SubPlotLabelXcoord = 0.085 #W7-X NBI case Electron scan, Ion scan, Neon scan
#SubPlotLabelXcoord = 0.40
#SubPlotLabelXcoord = 0.01
#SubPlotLabelXcoord = 0.003 #CONTOUR PLOT WITHOUT COLORBAR LHD discharge 113208 at t = 4.64 s
#SubPlotLabelXcoord = 0.006 #CONTOUR PLOT WITHOUT COLORBAR W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
#SubPlotLabelXcoord = 0.28 #LHD113208t4640 Carbon scan
#SubPlotLabelXcoord = 0.198 #LHD113208t4640 Electron scan, Ion scan, Helium scan
SubPlotLabelYcoord = -0.00071 #W7-X NBI case Neon scan
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

FilledErrors = False

ErrorBars = [False, False, True, True]
ErrorBarAlpha = 0.3

ShowLineAtXzero = True
ShowLineAtYzero = False

NoScientificAxes = True
