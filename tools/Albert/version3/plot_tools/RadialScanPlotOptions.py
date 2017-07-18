#!/usr/bin/env python

#####################################
####PLOT OPTIONS FOR RADIAL SCANS####
#####################################

##CAN BE OVERWRITTEN IN ORIGINAL FILE

import matplotlib

FigSize = (12,10)
#FigSize = (13.5,10) 
font = {'size':25}
matplotlib.rc('font', **font) 
matplotlib.rc('lines',markeredgewidth=0,markersize=3,linewidth=2.5)
matplotlib.rc('axes',linewidth=1.5)
matplotlib.rcParams['mathtext.default'] = 'it'
matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.preamble'].append(r'\usepackage{amsmath}')
#matplotlib.rcParams['text.latex.preamble'].append(r'\usepackage{sfmath}')

PlotLinespecs = ['rD-', 'bo-', 'mv-.', 'cs-.', 'yx--', 'g+--', 'kp:', 'k*:'] 
PlotLineColors = ['r', 'b', 'm', 'c', 'y', 'g', 'k', '0.75'] 
PlotMarkerEdgeWidth = [3, 3, 3, 0, 3, 3, 3, 3] 
PlotMarkerSize = 10 

xAxisScale = 'linear' #'linear', 'log' or 'symlog'
yAxisScale = 'linear' #'linear', 'log' or 'symlog' 

ShowGrid = True

xAxisLabel = r'$r/a$' 
yAxisLabel = r'$<\mathbf{\Gamma}_{e} \cdot \nabla r> $ $[10^{20} \mathrm{m}^{-2} \mathrm{s}^{-1}]$' 
#yAxisLabel = r'$E_r$ $[\mathrm{kV/m}]$'
AxesLabelSize = 30

AxisLimAuto = True 
xAxisLim = [0.18, 0.92]
yAxisLim = [-3.2, 0.0]

PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Pitch-angle scattering w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$', r'SFINCS Pitch-angle scattering w/ $\Phi_1$'] 
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$']

#PlotLegendLabels = []

#LegendFontSize = 15
#LegendProperties = {'weight':'bold'}
LegendProperties ={'size':'20', 'weight':'heavy'}
LegendPosition = 3
LegendNumberColumns = 1
#LegendBBoxToAnchor = [0.01, 0.01, 1., .102]
#LegendBBoxToAnchor = [0.435, 0.765, 1., .102]
LegendBBoxToAnchor = [0.435, 0.01, 1., .102]

