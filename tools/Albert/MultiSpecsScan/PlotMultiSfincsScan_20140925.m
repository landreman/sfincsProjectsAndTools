function PlotMultiSfincsScan_20140925()

% Name of .h5 HDF5 file from SFINCS:
h5filename='sfincsOutput.h5';

excludeRunsThatDidntConverge = true;
%excludeRunsThatDidntConverge = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quantities to show
%QuantitiesToPlot = {'ArrayFirstSpeciesParticleFluxCoefficients'; 'particleFlux'};
QuantitiesToPlotArray = cellstr({'ArrayFirstSpeciesParticleFluxCoefficients'; 'particleFlux'})
%QuantitiesOnXaxis = {'nHats'; 'd(PhiHat)d(psi_N)'};
QuantitiesOnXaxisArray = cellstr({'nHats'; 'd(PhiHat)d(psi_N)'})
XindicesToPlot = {1 1}; %This is an array of integers which controls which array element to use in the plots for parameters that are vectors in SFINCS.
%XindicesToPlot must be of the same length as QuantitiesOnXaxisArray and only
%contain positive integers. '1' should be used if the parameter is a scalar
%in SFINCS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot options
FigureWindowSize = 15; 
figureOffset = 0;

PlotColors = [1,0,0;  0.8,0.6,0;  0,0.7,0;  0,0.8,0.9;  0,0,1;  1,0,1;  0.6,0.6,0.6;  0,0,0];

PlotLinespecs = {'.-r', 'o-g', 'x-b', '.-m', '*-c', '.-r', '.-r', '.-b', '.-m'};

yAxesLabels = {'L_{11}^{zz}', 'L_{11}^{zi}', 'L_{12}^{zz} + L_{12}^{zi}', '\Gamma_{i}', '\Gamma_{z}'};
xAxesLabels = {'n_{i}', 'd \Phi / d \psi'};


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
%RHSMode1s = 1;
%RHSMode1s = 0;
%RHSMode2s = 0;
files = dir();
dumpedFieldsYet = false;
didItConverges = [];
%Nthetas = [];
%Nzetas = [];
%Nxis = [];
%NLs = [];
%Nxs = [];
%NxPotentialsPerVths = [];
%xMaxs = [];
%log10tols = [];
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

    %{
    RHSMode = h5read(filename,[location,'RHSMode']);
    switch RHSMode
        case 1
            RHSMode1s = RHSMode1s + 1;
        case 2
            RHSMode2s = RHSMode2s + 1;            
        otherwise
            error('Unrecognized RHSMode')
    end
    if RHSMode1s > 0 && RHSMode2s > 0
        error('Runs must either all have RHSMode=1 or RHSMode=2')
    end
     %}

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

%return
%PlotYarray
%PlotYarray{1,1,1}
%PlotYarray{1,1,2}
%PlotYarray{1,1,3}
%PlotYarray{2,1,1}
%PlotYarray{2,1,2}
%PlotYarray{2,1,3}
%PlotYarray{1,2,1}
%PlotYarray{1,2,2}
%PlotYarray{1,2,3}
%PlotYarray{2,2,1}
%PlotYarray{2,2,2}
%PlotYarray{2,2,3}
%PlotXarray
%PlotXarray{1,1,1}
%PlotXarray{1,1,2}
%PlotXarray{1,1,3}
%return
    
    %Nthetas(runNum) = h5read(filename,[location,'Ntheta']);
    %Nzetas(runNum) = h5read(filename,[location,'Nzeta']);
    %Nxis(runNum) = h5read(filename,[location,'Nxi']);
    %NLs(runNum) = h5read(filename,[location,'NL']);
    %Nxs(runNum) = h5read(filename,[location,'Nx']);
    %NxPotentialsPerVths(runNum) = h5read(filename,[location,'NxPotentialsPerVth']);
    %xMaxs(runNum) = h5read(filename,[location,'xMax']);
    %log10tols(runNum) = -log10(h5read(filename,[location,'solverTolerance']));
    %fprintf('%d %d %d %d %d %g %g %g\n',Nthetas(runNum), ...
    %        Nzetas(runNum),Nxis(runNum),NLs(runNum),Nxs(runNum), ...
    %        NxPotentialsPerVths(runNum), xMaxs(runNum), log10tols(runNum))

    %outputs(runNum,((1:Nspecies)-1)*3+1) = h5read(filename,[location,'particleFlux']);
    %outputs(runNum,((1:Nspecies)-1)*3+2) = h5read(filename,[location,'heatFlux']);
    %outputs(runNum,((1:Nspecies)-1)*3+3) = h5read(filename,[location,'FSABFlow']);
    %if Nspecies == 1
    %    outputs(runNum,4) = h5read(filename,[location,'didItConverge']);
    %    outputs(runNum,5) = h5read(filename,[location,'elapsed time (s)']);
    %end
    
    %{
    if RHSMode1s > 0
        outputs(runNum,1) = h5read(filename,[location,'particleFlux']);
        outputs(runNum,2) = h5read(filename,[location,'heatFlux']);
        outputs(runNum,3) = h5read(filename,[location,'FSABFlow']);
        outputs(runNum,4) = h5read(filename,[location,'didItConverge']);
        outputs(runNum,5) = h5read(filename,[location,'elapsed time (s)']);
    else
        transportMatrix = h5read(filename,[location,'transportMatrix']);
        outputs(runNum,1) = transportMatrix(1,1);
        outputs(runNum,2) = transportMatrix(1,2);
        outputs(runNum,3) = transportMatrix(1,3);
        outputs(runNum,4) = transportMatrix(2,1);
        outputs(runNum,5) = transportMatrix(2,2);
        outputs(runNum,6) = transportMatrix(2,3);
        outputs(runNum,7) = transportMatrix(3,1);
        outputs(runNum,8) = transportMatrix(3,2);
        outputs(runNum,9) = transportMatrix(3,3);
        
    end
    %}
    
%end

PlotYarray
PlotXarray

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



%NumberOfPlotRows
%NumberOfPlotCols

plotNumber = 0;

PlotYarray{1,1,:}
PlotYarray{1,1,2}(3)
PlotArrr = PlotYarray{1,1,:}
PlotArrr(1,:)
PlotXarray{1,1,:}

%return

Xarray = [];
Yarray = [];

for runNum=1:numRuns
    Xarray(runNum) = PlotXarray{1,1,runNum}(2);
    Yarray(runNum) = PlotYarray{1,1,runNum}(2);
end

Xarray
Yarray

%return

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
            
            plotNumber = plotNumber + 1;
            plotCol = rem(plotNumber - 1,NumberOfPlotCols) + 1;
            plotRow = ceil(plotNumber / NumberOfPlotCols);
            
            PlotColor = [0 0 0];
            PlotLine = 'o-k';
            
            xAxesLabel = '';
            yAxesLabel = '';
            
            if plotNumber <= length(PlotColors) 
                PlotColor = PlotColors(plotNumber,:);
            end
            if plotNumber <= length(PlotLinespecs) 
                PlotLine = PlotLinespecs{plotNumber};
            end
            
            if k <= length(xAxesLabels) 
                xAxesLabel = xAxesLabels{k};
            end
            if PosYaxesLabel + l <= length(yAxesLabels) 
                yAxesLabel = yAxesLabels{PosYaxesLabel + l};
            end
            
            subplot('Position',[leftMargin + (plotCol - 1)*(plotWidth + interPlotHorizontalSpacing), bottomMargin + (NumberOfPlotRows - plotRow)*(plotHeight + interPlotVerticalSpacing), plotWidth, plotHeight])
            plot(Xarray, Yarray, PlotLine, 'linewidth', PlotLineWidth, 'MarkerSize', PlotMarkerSize, 'Color', PlotColor)
            
            set(gca,'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off', 'FontSize', TickFontSize)
            
            xlabel(xAxesLabel,'FontSize',xLabelSize)
            ylabel(yAxesLabel,'FontSize',yLabelSize)
            
            %length(PlotColors)
            %length(PlotLinespecs)
            %PlotColors = [1,0,0;  0.8,0.6,0;  0,0.7,0;  0,0.8,0.9;  0,0,1;  1,0,1;  0.6,0.6,0.6;  0,0,0];

            %PlotLinespecs = {'.-r', '.-g', '.-b', '.-m', '.-c', '.-r', '.-r', '.-b', '.-m'};

            %yAxesLabels = {};
            %xAxesLabels = {};
        end
    end
    PosYaxesLabel = PosYaxesLabel + length(PlotYarray{j,1,1});
end


return;


%maxs=ones(numQuantities,1)*(-1e10);
%mins=ones(numQuantities,1)*(1e10);
maxs=ones(numQuantities,1)*(-Inf);
mins=ones(numQuantities,1)*Inf;
for iParameter = 1:numParameters
    maxs = max([maxs, quantities{iParameter}'],[],2);
    mins = min([mins, quantities{iParameter}'],[],2);
end



figure(1+figureOffset)
clf
set(gcf,'Color','w')

numCols = numParameters;

for iQuantity = 1:numQuantities
    if maxs(iQuantity) <= mins(iQuantity)
        maxs(iQuantity) = mins(iQuantity)+1;
    end
    for iParameter = 1:numParameters
        subplot(numRows, numCols, iParameter  + (plotRows(iQuantity) - 1)*numParameters)
        plot(1./abscissae{iParameter}, quantities{iParameter}(:,iQuantity)', PlotLinespecs{iQuantity})
        hold on
        plot(1./[convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
        ylim([mins(iQuantity), maxs(iQuantity)])
        xlabel(['1/',parametersToVary{iParameter}])
        ylabel(yAxesLabels{plotRows(iQuantity)})
    end
end

temp=dbstack;
nameOfThisProgram=sprintf('%s',temp(1).file);
stringForTop = ['Convergence scan from fortran multi-species version of SFINCS, plotted using ',nameOfThisProgram];


annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
    'Interpreter','none','VerticalAlignment','bottom',...
    'FontSize',12,'LineStyle','none','String',stringForTop);

stringForBottom = ['Run in: ',pwd];

annotation('textbox',[0 0 1 .04],'HorizontalAlignment','center',...
           'Interpreter','none','VerticalAlignment','top',...
           'FontSize',12,'LineStyle','none','String', ...
           stringForBottom);
       
% --------------------------------------------------------   

    function l = getLocationString(runNum)
        l = sprintf('/run%3d/',runNum);
    end

end