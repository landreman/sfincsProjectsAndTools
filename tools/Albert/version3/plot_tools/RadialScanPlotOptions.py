#!/usr/bin/env python

#####################################
####PLOT OPTIONS FOR RADIAL SCANS####
#####################################

##CAN BE OVERWRITTEN IN ORIGINAL FILE

import matplotlib
import numpy 

FigSize = (12,10)
#FigSize = (13.5,10) 
font = {'size':35}
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

PlotLinespecs = ['rD-', 'bo-', 'mv-.', 'cs-.', 'yx--', 'g+--', 'kp:', 'k*:'] 
PlotLineColors = ['r', 'b', 'm', 'c', 'y', 'g', 'k', '0.75']
#PlotLinespecs = ['rD-', 'mv-.']
#PlotLineColors = ['r', 'm']
PlotMarkerEdgeWidth = [3, 3, 3, 0, 3, 3, 3, 3] 
PlotMarkerSize = 10 

xAxisScale = 'linear' #'linear', 'log' or 'symlog'
yAxisScale = 'linear' #'linear', 'log' or 'symlog' 

ShowGrid = True

xAxisLabel = r'$r/a$'
#yAxisLabel = r'$n_{\mathrm{C}^{6+}}$ $[10^{19} \mathrm{m}^{-3}]$'
#yAxisLabel = r'$n$ $[10^{19} \mathrm{m}^{-3}]$'
yAxisLabel = r'$T$ $[\mathrm{keV}]$'
#yAxisLabel = r'$\nu^{\prime}$'
# e   i   \mathrm{He}^{2+}   \mathrm{C}^{6+}   \mathrm{Ne}^{10+}
#yAxisLabel = r'$<\mathbf{\Gamma}_{e} \cdot \nabla r> $ $[10^{20} \mathrm{m}^{-2} \mathrm{s}^{-1}]$' 
#yAxisLabel = r'$E_r$ $[\mathrm{kV/m}]$'
AxesLabelSize = 40

AxisLimAuto = True 
xAxisLim = [0.18, 0.92]
yAxisLim = [-3.2, 0.0]
#xAxisLim = [0.10, 1.02]
#yAxisLim = [-25.0, 0.0]

xAxisLabelCoords = [0.5,-0.09]
yAxisLabelCoords = [-0.11,0.5]

LeftMargin = 0.15
RightMargin = 0.95
TopMargin = 0.95
BottomMargin = 0.15

ShowLegend = True

#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Pitch-angle scattering w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$', r'SFINCS Pitch-angle scattering w/ $\Phi_1$'] 
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$']
#PlotLegendLabels = [r'$n_e$', r'$n_i$', r'$n_{\mathrm{He}^{2+}}$', r'$n_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$T_e$', r'$T_i$', r'$T_{\mathrm{He}^{2+}}$', r'$T_{\mathrm{C}^{6+}}$']
PlotLegendLabels = [r'$T_e$', r'$T_i$', r'$T_{\mathrm{Ne}^{10+}}$']
#PlotLegendLabels = [r'$\nu_{ee}^{\prime}$', r'$\nu_{ii}^{\prime}$', r'$\nu_{HeHe}^{\prime}$', r'$\nu_{CC}^{\prime}$']
#PlotLegendLabels = [r'$\nu_{ee}^{\prime}$', r'$\nu_{ii}^{\prime}$', r'$\nu_{NeNe}^{\prime}$']
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$']

#PlotLegendLabels = []

#LegendFontSize = 15
#LegendProperties = {'weight':'bold'}
LegendProperties ={'size':'30', 'weight':'heavy'}
LegendPosition = 3
LegendNumberColumns = 1
#LegendBBoxToAnchor = [0.01, 0.01, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.765, 1., .102]
#LegendBBoxToAnchor = [0.435, 0.765, 1., .102]
#LegendBBoxToAnchor = [0.30, 0.01, 1., .102]
#LegendBBoxToAnchor = [0.13, 0.01, 1., .102]
#LegendBBoxToAnchor = [0.01, 0.655, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.665, 1., .102]
LegendBBoxToAnchor = [0.005, 0.35, 1., .102]

ShowSubPlotLabel = True
SubPlotLabel = '(b)'
SubPlotLabelXcoord = 0.085
#SubPlotLabelXcoord = 0.40
#SubPlotLabelYcoord = -0.00043
#SubPlotLabelYcoord = 0.1283
SubPlotLabelYcoord = 2.090
