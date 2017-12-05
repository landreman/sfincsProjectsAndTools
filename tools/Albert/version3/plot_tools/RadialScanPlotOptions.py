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

#FigSize = (12,10)
#FigSize = (13.5,10)
FigSize = (12,12)
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

#matplotlib.style.use('grayscale')


#PlotLinespecs = ['rD-', 'bo-', 'mv-.', 'cs-.', 'yx--', 'g+--', 'kp:', 'k*:']
#PlotLineColors = ['r', 'b', 'm', 'c', 'y', 'g', 'k', '0.75']
#PlotLinespecs = ['rD-', 'mv-.']
#PlotLineColors = ['r', 'm']
PlotLinespecs = ['r-', 'b--', 'g-.', 'y:']
PlotLineColors = ['r', 'b', 'g', '#AAAA00']
#PlotLinespecs = ['k-']
#PlotLineColors = ['#555555']
PlotMarkerEdgeWidth = [3, 3, 3, 0, 3, 3, 3, 3] 
PlotMarkerSize = 10
PlotLineWidth=5.0

xAxisScale = 'linear' #'linear', 'log' or 'symlog'
yAxisScale = 'log' #'linear', 'log' or 'symlog' 

ShowGrid = True

xAxisLabel = r'$r/a$'
#yAxisLabel = r'$n_{\mathrm{C}^{6+}}$ $[10^{19} \mathrm{m}^{-3}]$'
#yAxisLabel = r'$n$ $[10^{19} \mathrm{m}^{-3}]$'
#yAxisLabel = r'$T$ $[\mathrm{keV}]$'
yAxisLabel = r'$\nu_{ss}^{\prime}$'
# e   i   \mathrm{He}^{2+}   \mathrm{C}^{6+}   \mathrm{Ne}^{10+}
#yAxisLabel = r'$<\mathbf{\Gamma}_{e} \cdot \nabla r> $ $[10^{20} \mathrm{m}^{-2} \mathrm{s}^{-1}]$' 
#yAxisLabel = r'$E_r$ $[\mathrm{kV/m}]$'
#yAxisLabel = r'$Z_{\mathrm{eff}}$'
AxesLabelSize = 40

TickSize = 35

AxisLimAuto = False 
#xAxisLim = [0.18, 0.92]
#yAxisLim = [-3.2, 0.0]
#xAxisLim = [0.10, 1.02]
#yAxisLim = [-25.0, 0.0]
xAxisLim = [0.0, 1.0]
yAxisLim = [0.00007, 0.15]
#yAxisLim = [0.0, 1.5]
#yAxisLim = [0.0, 5.5]
#yAxisLim = [1.0, 2.5]

xAxisLabelCoords = [0.5,-0.09]
#yAxisLabelCoords = [-0.11,0.5]
yAxisLabelCoords = [-0.09,0.5] #CONTOUR PLOT WITHOUT COLORBAR

#LeftMargin = 0.15
LeftMargin = 0.135 #CONTOUR PLOT WITHOUT COLORBAR
#RightMargin = 0.95
RightMargin = 0.935 #CONTOUR PLOT WITHOUT COLORBAR
#TopMargin = 0.95
TopMargin = 0.97 #CONTOUR PLOT WITHOUT COLORBAR
#BottomMargin = 0.15
BottomMargin = 0.17 #CONTOUR PLOT WITHOUT COLORBAR

ShowLegend = False

#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Pitch-angle scattering w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$', r'SFINCS Pitch-angle scattering w/ $\Phi_1$'] 
#PlotLegendLabels = [r'SFINCS Fokker-Planck w/o $\Phi_1$', r'SFINCS Fokker-Planck w/ $\Phi_1$']
#PlotLegendLabels = [r'$n_e$', r'$n_i$', r'$n_{\mathrm{He}^{2+}}$', r'$n_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$n_e$', r'$n_i$', r'$n_{\mathrm{He}^{2+}}$', r'$10 \times n_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$T_e$', r'$T_i$', r'$T_{\mathrm{He}^{2+}}$', r'$T_{\mathrm{C}^{6+}}$']
#PlotLegendLabels = [r'$T_e$', r'$T_i$', r'$T_{\mathrm{Ne}^{10+}}$']
#PlotLegendLabels = [r'$T_e$', r'$T_i = T_{\mathrm{He}^{2+}} = T_{\mathrm{C}^{6+}}$']
PlotLegendLabels = [r'$\nu_{ee}^{\prime}$', r'$\nu_{ii}^{\prime}$', r'$\nu_{HeHe}^{\prime}$', r'$\nu_{CC}^{\prime}$']
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
#LegendBBoxToAnchor = [0.005, 0.35, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.47, 1., .102]
#LegendBBoxToAnchor = [0.005, 0.005, 1., .102]
LegendBBoxToAnchor = [0.005, 0.58, 1., .102]

ShowSubPlotLabel = False
SubPlotLabel = '(f)'
SubPlotLabelSize = 40
#SubPlotLabelXcoord = 0.085
#SubPlotLabelXcoord = 0.40
#SubPlotLabelXcoord = 0.01
#SubPlotLabelXcoord = 0.003 #CONTOUR PLOT WITHOUT COLORBAR LHD discharge 113208 at t = 4.64 s
SubPlotLabelXcoord = 0.006 #CONTOUR PLOT WITHOUT COLORBAR W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
#SubPlotLabelXcoord = 0.28
#SubPlotLabelYcoord = -0.00043
#SubPlotLabelYcoord = 0.1283
#SubPlotLabelYcoord = 0.08
#SubPlotLabelYcoord = 2.090
#SubPlotLabelYcoord = 1.422
#SubPlotLabelYcoord = 5.2
#SubPlotLabelYcoord = 5.84 #CONTOUR PLOT WITHOUT COLORBAR LHD discharge 113208 at t = 4.64 s
SubPlotLabelYcoord = 5.715 #CONTOUR PLOT WITHOUT COLORBAR W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
#SubPlotLabelYcoord = 5.55 #CONTOUR PLOT WITHOUT COLORBAR W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
#SubPlotLabelYcoord = 2.422
#SubPlotLabelYcoord = 0.1
