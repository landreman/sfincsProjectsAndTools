function PlotMultiSfincsScan_20140925()

% Name of .h5 HDF5 file from SFINCS:
h5filename='sfincsOutput.h5';

excludeRunsThatDidntConverge = true;
%excludeRunsThatDidntConverge = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quantities to show
%QuantitiesToPlot = {'ArrayFirstSpeciesParticleFluxCoefficients'; 'particleFlux'};
%QuantitiesToPlotArray = cellstr({'ArrayFirstSpeciesParticleFluxCoefficients'; 'particleFlux'})
QuantitiesToPlotArray = cellstr({'ArrayFirstSpeciesParticleFluxCoefficients'})
%QuantitiesOnXaxis = {'nHats'; 'd(PhiHat)d(psi_N)'};
%QuantitiesOnXaxisArray = cellstr({'nHats'; 'd(PhiHat)d(psi_N)'})
QuantitiesOnXaxisArray = cellstr({'nHats'})
%XindicesToPlot = {1 1};
XindicesToPlot = {2}; %This is an array of integers which controls which array element to use in the plots for parameters that are vectors in SFINCS.
%XindicesToPlot must be of the same length as QuantitiesOnXaxisArray and only
%contain positive integers. '1' should be used if the parameter is a scalar
%in SFINCS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot options
set(0, 'defaultTextInterpreter', 'tex'); 

PlotLabel = {'(impurityFlux) = (sfincsFactor) \{ L_{11}^{zz} [(1/n_{z}) (dn_{z}/d\psi) + (Ze/T_{z}) (d\Phi/d\psi) - (3/2T_{z}) (dT_{z}/d\psi)] +', ...
    '+ L_{11}^{zi} [(1/n_{i}) (dn_{i}/d\psi) + (e/T_{i}) (d\Phi/d\psi) - (3/2T_{i}) (dT_{i}/d\psi)] + L_{12}^{zz} [(1/T_{z})(dT_{z}/d\psi)] + L_{12}^{zi} [(1/T_{i})(dT_{i}/d\psi)]\}'};

FigureWindowSize = 15; 
figureOffset = 0;

PlotColors = [1,0,0;  0.8,0.6,0;  0,0.7,0;  0,0.8,0.9;  0,0,1;  1,0,1;  0.6,0.6,0.6;  0,0,0];

PlotLinespecs = {'d-r', 'o-g', 'x-b', '.-m', '*-c', '.-r', '.-r', '.-b', '.-m'};

PlotXLogs = {'log', 'log', 'log'}; %'linear' or 'log'
PlotYLogs = {'log', 'log', 'log'}; %'linear' or 'log'

yAxesLabels = {'L_{11}^{zz}', 'L_{11}^{zi}', 'L_{12}^{zz} + L_{12}^{zi}', '\Gamma_{i}', '\Gamma_{z}'};
xAxesLabels = {'n_{z} (= C \cdot n_{i})', 'd \Phi / d \psi'};


leftMargin = 0.1;
rightMargin = 0.02;
topMargin = 0.065;
bottomMargin = 0.1;
interPlotHorizontalSpacing = 0.08;
interPlotVerticalSpacing = 0.1;


xLabelSize = 15;
yLabelSize = 15;
PlotLineWidth = 2.5;
PlotMarkerSize = 8;
TickFontSize = 13;

%Not used at the moment%
LegendSize = 11;
LabelSize = 16;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if length(QuantitiesOnXaxisArray) ~= length(XindicesToPlot)
    error('Arrays QuantitiesOnXaxis and XindicesToPlot must have the same length');
end

assert(length(QuantitiesToPlotArray) > 0, 'Quantities to plot empty');
assert(length(QuantitiesOnXaxisArray) > 0, 'Parameter quantities empty');

numRuns = 0;

files = dir();
dumpedFieldsYet = false;
didItConverges = [];
outputs = [];
Nspecies = -1;

validFilenames = {};

for iFile = 1:size(files,1)
    if ~ files(iFile).isdir
        continue
    end

    % Skip the . and .. directories
    if strcmp(files(iFile).name,'.') | strcmp(files(iFile).name,'..')
        continue
    end

    % Try to open an HDF5 file
    filename = [files(iFile).name, '/', h5filename];
    try
        info = h5info(filename);
        fprintf('Successfully opened h5 file %s\n',filename)
    catch
        fprintf('Did not succeed in opening h5 file %s\n',filename)
        continue
    end

    if ~ dumpedFieldsYet
        fprintf('Fields saved in the HDF5 files:\n')
        for i=1:numel(info.Groups(1).Datasets)
            fprintf('  %s\n',info.Groups(1).Datasets(i).Name)
        end
        dumpedFieldsYet = true;
    end

    programMode = h5read(filename,'/programMode');
    if programMode ~= 1
        fprintf('Ignoring this run since programMode is not 1.\n')
        continue
    end

    location  = getLocationString(1);
    integerToRepresentTrue = h5read(filename,[location,'integerToRepresentTrue']);
    didItConverge = h5read(filename,[location,'didItConverge']);
    if excludeRunsThatDidntConverge && (didItConverge ~= integerToRepresentTrue)
        fprintf('Ignoring this run since it did not converge.\n')
        continue
    end


    % If we made it this far, then let's count the run.
    numRuns = numRuns + 1;
    %validFilenames(numRuns) = filename
    validFilenames{end+1} = filename;
    didItConverges(numRuns) = (didItConverge == integerToRepresentTrue);
    if didItConverge ~= integerToRepresentTrue
        beep
        fprintf('Warning: run %d did not converge.\n',i)
    end


    Nspecies_new = h5read(filename,[location,'Nspecies']);
    if Nspecies < 0
        Nspecies = Nspecies_new;
    else
        if Nspecies ~= Nspecies_new
            error('Number of species is not consistent among runs')
        end
    end
end

if numRuns < 1
    error('No runs found!')
end

PlotYarray = cell(length(QuantitiesToPlotArray), length(QuantitiesOnXaxisArray), numRuns);
PlotXarray = cell(length(QuantitiesToPlotArray), length(QuantitiesOnXaxisArray), numRuns);

NumberOfFigs = 0;

for j=1:length(QuantitiesToPlotArray)
    for k=1:length(QuantitiesOnXaxisArray)
        for runNum=1:numRuns
            filename = validFilenames{runNum};
            %data = h5read(filename,[location,char(QuantitiesToPlotArray(j))]);
            %PlotYarray(runNum) = (h5read(filename,[location,char(QuantitiesToPlotArray(j))]))(1);
            PlotYarray{j,k,runNum} = h5read(filename,[location,char(QuantitiesToPlotArray(j))]);
            %char(QuantitiesToPlotArray(j))
            %test = h5read(filename,[location,'ArrayFirstSpeciesParticleFluxCoefficients'])
            
            %PlotXarray(runNum) = h5read(filename,[location,QuantitiesOnXaxisArray(k)]);
            PlotXarray{j,k,runNum} = h5read(filename,[location,char(QuantitiesOnXaxisArray(k))]);
        end
        NumberOfFigs = NumberOfFigs + length(PlotYarray{j,k,1});
    end
end

%NumberOfPlotRows = ceil(NumberOfFigs / sqrt(NumberOfFigs));
%NumberOfPlotCols = ceil(NumberOfFigs / NumberOfPlotRows);

NumberOfPlotCols = ceil(NumberOfFigs / sqrt(NumberOfFigs));
NumberOfPlotRows = ceil(NumberOfFigs / NumberOfPlotCols);


%%PLOT 
FigPlot = figure(1)
clf

set(gcf,'Color','w')
set(gcf,'Units','inches','Position',[1,1,1.5*FigureWindowSize*NumberOfPlotCols/(NumberOfPlotCols + NumberOfPlotRows),1.0*FigureWindowSize*NumberOfPlotRows/(NumberOfPlotCols + NumberOfPlotRows)])

plotHeight = (1-topMargin-bottomMargin-(NumberOfPlotRows - 1)*interPlotVerticalSpacing)/NumberOfPlotRows;
%plotBottom_topRow = bottomMargin + plotHeight + interPlotVerticalSpacing;

plotWidth = (1-leftMargin-rightMargin-(NumberOfPlotCols - 1)*interPlotHorizontalSpacing)/NumberOfPlotCols;
%plotLeft_middleColumn = leftMargin + plotWidth + interPlotHorizontalSpacing;
%plotLeft_rightColumn = leftMargin + 1*(plotWidth + interPlotHorizontalSpacing);

plotNumber = 0;

Xarray = [];
Yarray = [];

for runNum=1:numRuns
    Xarray(runNum) = PlotXarray{1,1,runNum}(2);
    Yarray(runNum) = PlotYarray{1,1,runNum}(2);
end


PosYaxesLabel = 0;
for j=1:length(QuantitiesToPlotArray)
    for k=1:length(QuantitiesOnXaxisArray)
        for l=1:length(PlotYarray{j,k,1})
            Xarray = [];
            Yarray = [];
            for runNum=1:numRuns
                Xarray(runNum) = PlotXarray{j,k,runNum}(XindicesToPlot{k});
                Yarray(runNum) = PlotYarray{j,k,runNum}(l);
            end
            
            [XarraySorted, SortIndex] = sort(Xarray);
            YarraySorted = Yarray(SortIndex);
            
            plotNumber = plotNumber + 1;
            plotCol = rem(plotNumber - 1,NumberOfPlotCols) + 1;
            plotRow = ceil(plotNumber / NumberOfPlotCols);
            
            PlotColor = [0 0 0];
            PlotLine = 'o-k';
            PlotXLog = 'linear';
            PlotYLog = 'linear';
            
            xAxesLabel = '';
            yAxesLabel = '';
            
            if plotNumber <= length(PlotColors) 
                PlotColor = PlotColors(plotNumber,:);
            end
            if plotNumber <= length(PlotLinespecs) 
                PlotLine = PlotLinespecs{plotNumber};
            end
            if plotNumber <= length(PlotXLogs) 
                PlotXLog = PlotXLogs{plotNumber};
            end
            if plotNumber <= length(PlotYLogs) 
                PlotYLog = PlotYLogs{plotNumber};
            end
            
            if k <= length(xAxesLabels) 
                xAxesLabel = xAxesLabels{k};
            end
            if PosYaxesLabel + l <= length(yAxesLabels) 
                yAxesLabel = yAxesLabels{PosYaxesLabel + l};
            end
            
            subplot('Position',[leftMargin + (plotCol - 1)*(plotWidth + interPlotHorizontalSpacing), bottomMargin + (NumberOfPlotRows - plotRow)*(plotHeight + interPlotVerticalSpacing), plotWidth, plotHeight])
            plot(XarraySorted, YarraySorted, PlotLine, 'linewidth', PlotLineWidth, 'MarkerSize', PlotMarkerSize, 'Color', PlotColor)
            
            set(gca,'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off', 'FontSize', TickFontSize, 'XScale', PlotXLog, 'YScale', PlotYLog)
            
            xlabel(xAxesLabel,'FontSize',xLabelSize)
            ylabel(yAxesLabel,'FontSize',yLabelSize)
            
        end
    end
    PosYaxesLabel = PosYaxesLabel + length(PlotYarray{j,1,1});
end

suptitle(PlotLabel);

return;
  
% --------------------------------------------------------   

    function l = getLocationString(runNum)
        l = sprintf('/run%3d/',runNum);
    end

end