clear

savePDFs = true;
%savePDFs = false;

figureFilename = 'SFINCS_DKES_comparison';

addpath('/global/homes/l/landrema/MATLAB')

titleFontSize = 18;

colors = [1,0,0;
    0.7,0.5,0;
    0,0.7,0;
    0,0,1];

% Warning: this script extracts r/a from the directory names and doesn't actually
% confirm that r/a for each run matches the directory name.

h5filename='sfincsOutput.h5';

%profilesFilename='C:\Users\landreman\Box Sync\MATLAB\dataForSfincs\profiles.txt';
profilesFilename='/global/scratch2/sd/landrema/W7X_files_from_Yuriy/profiles.txt';

Nr = 0;
directoryEntries1 = dir();
r_over_as = [];
B0s = [];
ionFluxes = {};
electronFluxes = {};
ionHeatFluxes = {};
electronHeatFluxes = {};
ionFlows = {};
electronFlows = {};
dPhiHatdpsiNs = {};
location = '/run  1/';
a = 0.52625; % W7-X minor radius
psiAHat = 0;
B00_0 = 2.380; % (0,0) Boozer component of |B| on the magnetic axis.
B00_1 = 2.487; % (0,0) Boozer component of |B| at the last closed magnetic surface.


% ---------------------------------------------
% Read in profile info from Yuriy:

data=importdata(profilesFilename);

r_over_as_Yuriy = data.data(:,1);
ne_Yuriy = data.data(:,3);
ni_Yuriy = data.data(:,4);
Te_Yuriy = data.data(:,5);
Ti_Yuriy = data.data(:,6);
Er_Yuriy = data.data(:,7);
particleFlux_Yuriy = data.data(:,8);
Qe_Yuriy = data.data(:,10);
Qi_Yuriy = data.data(:,11);
je_Yuriy = data.data(:,12);
ji_Yuriy = data.data(:,13);

for idir1 = 1:numel(directoryEntries1)
    if ~ directoryEntries1(idir1).isdir
        continue
    end

    % Skip the . and .. directories
    radiusDirName = directoryEntries1(idir1).name;
    if strcmp(radiusDirName,'.') || strcmp(radiusDirName,'..')
        continue
    end

    % If we made it this far, then count this directory.

    cd(radiusDirName)
    didAnyRunsWorkAtThisRadius = false;

    for iPhysics = 1:4
        switch iPhysics
          case 1
            physicsDirName='FP_full';
          case 2
            physicsDirName='FP_DKES';
          case 3
            physicsDirName='PAS_full';
          case 4
            physicsDirName='PAS_DKES';
          otherwise
            error('Should not get here')
        end
        try
            cd(physicsDirName)
        catch
            fprintf('Warning: directory %s/%s does not exist.\n',radiusDirName,physicsDirName) 
            continue
        end

        directoryEntries2 = dir();
        iEr = 0;
        for idir2 = 1:numel(directoryEntries2)
            if ~ directoryEntries2(idir2).isdir
                continue
            end
            
            % Skip the . and .. directories
            ErDirName = directoryEntries2(idir2).name;
            if strcmp(ErDirName,'.') || strcmp(ErDirName,'..')
                continue
            end
            
            filename = [ErDirName, '/', h5filename];
            try
                info = h5info(filename);
                fprintf('Succeeded in reading %s/%s/%s\n',radiusDirName,physicsDirName,ErDirName);
            catch
                fprintf('Unable to read %s/%s/%s\n',radiusDirName,physicsDirName,ErDirName);
                continue
            end

            iEr = iEr + 1;

            if ~didAnyRunsWorkAtThisRadius
                % If we get here, this must be the first successful run at this radius.
                didAnyRunsWorkAtThisRadius = true;
                Nr = Nr + 1;
                r_over_a = sscanf(radiusDirName(10:end),'%g');
                r_over_as(Nr) = r_over_a;
                dPhiHatdpsiNs{Nr} = cell(4,1);
                ionFluxes{Nr} = cell(4,1);
                electronFluxes{Nr} = cell(4,1);
                ionHeatFluxes{Nr} = cell(4,1);
                electronHeatFluxes{Nr} = cell(4,1);
                ionFlows{Nr} = cell(4,1);
                electronFlows{Nr} = cell(4,1);
                for iPhysics2 = 1:4
                    dPhiHatdpsiNs{Nr}{iPhysics2} = [];
                    ionFluxes{Nr}{iPhysics2} = [];
                    electronFluxes{Nr}{iPhysics2} = [];
                    ionHeatFluxes{Nr}{iPhysics2} = [];
                    electronHeatFluxes{Nr}{iPhysics2} = [];
                    ionFlows{Nr}{iPhysics2} = [];
                    electronFlows{Nr}{iPhysics2} = [];
                end
                B0s(Nr) = h5read(filename,[location,'B0OverBBar']);
            end

            if psiAHat == 0
                psiAHat = h5read(filename,[location,'psiAHat']);
            end
            dPhiHatdpsiNs{Nr}{iPhysics}(iEr) = h5read(filename,[location,'d(PhiHat)d(psi_N)']);

            temp = h5read(filename,[location,'particleFlux']);
            ionFluxes{Nr}{iPhysics}(iEr) = temp(1);
            electronFluxes{Nr}{iPhysics}(iEr) = temp(2);
            
            temp = h5read(filename,[location,'heatFlux']);
            ionHeatFluxes{Nr}{iPhysics}(iEr) = temp(1);
            electronHeatFluxes{Nr}{iPhysics}(iEr) = temp(2);

            temp = h5read(filename,[location,'FSABFlow']);
            ionFlows{Nr}{iPhysics}(iEr) = temp(1);
            electronFlows{Nr}{iPhysics}(iEr) = temp(2);
        end

        cd('..')
    end

    cd('..')
end

%r_over_as

linespecsIons={'o-','+-','x-','.-'};
linespecsElectrons={'o:','+:','x:','.:'};
linespecs={'o','+','x','.'};

numRoots = zeros(4,Nr);
ambipolar_dPhiHatdpsiNs = zeros(4,Nr);
ambipolar_Ers = zeros(4,Nr);
particleFluxes_sfincs = zeros(4,Nr);
particleFluxes_sfincsDimensional = zeros(4,Nr);
ionHeatFluxes_sfincs = zeros(4,Nr);
ionHeatFluxes_sfincsDimensional = zeros(4,Nr);
electronHeatFluxes_sfincs = zeros(4,Nr);
electronHeatFluxes_sfincsDimensional = zeros(4,Nr);
ji_sfincs = zeros(4,Nr);
ji_sfincsDimensional = zeros(4,Nr);
je_sfincs = zeros(4,Nr);
je_sfincsDimensional = zeros(4,Nr);

vBar = 437695;  % = sqrt(2*(1.6e-19)*1000/(1.67e-27))
mBar = 1.67262178e-27; % proton mass
e = 1.60217657e-19;

for ir=1:Nr
    
    ErFactor = -2 / a * r_over_as(ir);

    %{
    % If r used to define the radial fluxes is a * sqrt(psi / psi_LCFS):
    particleFluxFactor = a * vBar / (2 * psiAHat * r_over_as(ir));
    heatFluxFactor = (1e-6) * (1e20) * mBar * a * vBar * vBar * vBar / (2 * psiAHat * r_over_as(ir));
    %}
    
    % If r used to define the radial fluxes is sqrt(psi / [pi * B_00(0)]):
    particleFluxFactor = vBar / (sqrt(2*psiAHat*B00_0) * r_over_as(ir));
    heatFluxFactor = (1e-6) * (1e20) * mBar * vBar * vBar * vBar / (sqrt(2*psiAHat*B00_0) * r_over_as(ir));
    
    %{
    % If r used to define the radial fluxes is sqrt(psi / [pi * B_00(1)]):
    particleFluxFactor = vBar / (sqrt(2*psiAHat*B00_1) * r_over_as(ir));
    heatFluxFactor = (1e-6) * (1e20) * mBar * vBar * vBar * vBar / (sqrt(2*psiAHat*B00_1) * r_over_as(ir));
    %}

    jFactor = -e*(1e20)*vBar/B0s(ir);
    %jFactor = -e*(1e20)*vBar/B0s(ir)*0.88117;
    %jFactor = -e*(1e20)*vBar/2.38;
    % In jFactor, I'm not sure why the minus sign must be included, but it seems necessary to get reasonable agreement with Turkin.
    
    for iPhysics = 1:4
    %for iPhysics = 2
        if numel(dPhiHatdpsiNs{ir}{iPhysics}) < 1
            continue
        end

        % Sort results in order of increasing E_r:
        [~,permutation] = sort(dPhiHatdpsiNs{ir}{iPhysics});
        dPhiHatdpsiNs{ir}{iPhysics} = dPhiHatdpsiNs{ir}{iPhysics}(permutation);
        ionFluxes{ir}{iPhysics} = ionFluxes{ir}{iPhysics}(permutation);
        electronFluxes{ir}{iPhysics} = electronFluxes{ir}{iPhysics}(permutation);
        ionHeatFluxes{ir}{iPhysics} = ionHeatFluxes{ir}{iPhysics}(permutation);
        electronHeatFluxes{ir}{iPhysics} = electronHeatFluxes{ir}{iPhysics}(permutation);
        ionFlows{ir}{iPhysics} = ionFlows{ir}{iPhysics}(permutation);
        electronFlows{ir}{iPhysics} = electronFlows{ir}{iPhysics}(permutation);


        dPhiHatdpsiNs_fine = linspace(min(dPhiHatdpsiNs{ir}{iPhysics}), max(dPhiHatdpsiNs{ir}{iPhysics}), 100);
        radialCurrent = interp1(dPhiHatdpsiNs{ir}{iPhysics}, ionFluxes{ir}{iPhysics} - electronFluxes{ir}{iPhysics}, dPhiHatdpsiNs_fine);
        signRadialCurrent = sign(radialCurrent);
        changes = abs(signRadialCurrent(2:end) - signRadialCurrent(1:(end-1)));
        numRoots(ir,iPhysics) = numel(changes);
        if sum(changes)>1
            %fprintf('Multiple roots for Er!!\n')
        end
        fprintf('sum(changes) = %d\n',sum(changes))
        
        % Solve for E_r:
        x0 = [min(dPhiHatdpsiNs{ir}{iPhysics}), max(dPhiHatdpsiNs{ir}{iPhysics})];
        radialCurrent_x = dPhiHatdpsiNs{ir}{iPhysics};
        radialCurrent_y = ionFluxes{ir}{iPhysics} - electronFluxes{ir}{iPhysics};
        try
            ambipolar_dPhiHatdpsiN = fzero(@(dPhiHatdpsi) interp1(radialCurrent_x, radialCurrent_y, dPhiHatdpsi), x0);
        catch
            fprintf('Warning: unable to solve for ambipolar Er.\n')
            ambipolar_dPhiHatdpsiN = mean(radialCurrent_x);
        end
        ambipolar_dPhiHatdpsiNs(iPhysics,ir) = ambipolar_dPhiHatdpsiN;
        ambipolar_Ers(iPhysics,ir) = ambipolar_dPhiHatdpsiNs(iPhysics,ir) * ErFactor;
        
        % Interpolate to get other quantities at the ambipolar E_r:
        particleFluxes_sfincs(iPhysics,ir) = interp1(dPhiHatdpsiNs{ir}{iPhysics}, ionFluxes{ir}{iPhysics}, ambipolar_dPhiHatdpsiN);
        particleFluxes_sfincsDimensional(iPhysics,ir) = particleFluxes_sfincs(iPhysics,ir) * particleFluxFactor;

        ionHeatFluxes_sfincs(iPhysics,ir) = interp1(dPhiHatdpsiNs{ir}{iPhysics}, ionHeatFluxes{ir}{iPhysics}, ambipolar_dPhiHatdpsiN);
        ionHeatFluxes_sfincsDimensional(iPhysics,ir) = ionHeatFluxes_sfincs(iPhysics,ir) * heatFluxFactor;

        electronHeatFluxes_sfincs(iPhysics,ir) = interp1(dPhiHatdpsiNs{ir}{iPhysics}, electronHeatFluxes{ir}{iPhysics}, ambipolar_dPhiHatdpsiN);
        electronHeatFluxes_sfincsDimensional(iPhysics,ir) = electronHeatFluxes_sfincs(iPhysics,ir) * heatFluxFactor;

        ji_sfincs(iPhysics,ir) = interp1(dPhiHatdpsiNs{ir}{iPhysics}, ionFlows{ir}{iPhysics}, ambipolar_dPhiHatdpsiN);
        ji_sfincsDimensional(iPhysics,ir) = ji_sfincs(iPhysics,ir) * jFactor;

        je_sfincs(iPhysics,ir) = interp1(dPhiHatdpsiNs{ir}{iPhysics}, electronFlows{ir}{iPhysics}, ambipolar_dPhiHatdpsiN);
        je_sfincsDimensional(iPhysics,ir) = je_sfincs(iPhysics,ir) * (-jFactor);
    end

    
    % ************************************************
    % Plot results in sfincs units:
    % ************************************************
    
    figure(ir)
    clf
    set(gcf,'Color','w')
    numRows = 2;
    numCols = 4;
    
    subplot(numRows,numCols,1)
    for iPhysics=1:4
        plot(dPhiHatdpsiNs{ir}{iPhysics}, ionFluxes{ir}{iPhysics}, linespecsIons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('dPhiHatdpsiN')
    ylabel('ion particle flux (sfincs units)')
    
    subplot(numRows,numCols,5)
    for iPhysics=1:4
        plot(dPhiHatdpsiNs{ir}{iPhysics}, electronFluxes{ir}{iPhysics}, linespecsElectrons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('dPhiHatdpsiN')
    ylabel('electron particle flux (sfincs units)')
    
    subplot(numRows,numCols,4)
    for iPhysics=1:4
        if numel(dPhiHatdpsiNs{ir}{iPhysics}) > 1
            [~,permutation] = sort(dPhiHatdpsiNs{ir}{iPhysics});
            plot(dPhiHatdpsiNs{ir}{iPhysics}(permutation), ionFluxes{ir}{iPhysics}(permutation), linespecs{iPhysics}, 'Color',colors(iPhysics,:))
            hold on
            
            dPhiHatdpsiN_fine = linspace(min(dPhiHatdpsiNs{ir}{iPhysics}), max(dPhiHatdpsiNs{ir}{iPhysics}), 100);
            plot(dPhiHatdpsiN_fine, interp1(dPhiHatdpsiNs{ir}{iPhysics}, ionFluxes{ir}{iPhysics}, dPhiHatdpsiN_fine, 'cubic'), '-', 'Color',colors(iPhysics,:))
            
            plot(dPhiHatdpsiNs{ir}{iPhysics}(permutation), electronFluxes{ir}{iPhysics}(permutation), linespecs{iPhysics}, 'Color',colors(iPhysics,:))
            plot(dPhiHatdpsiN_fine, interp1(dPhiHatdpsiNs{ir}{iPhysics}, electronFluxes{ir}{iPhysics}, dPhiHatdpsiN_fine, 'cubic'), ':', 'Color',colors(iPhysics,:))

            plot(ambipolar_dPhiHatdpsiNs(iPhysics,ir), particleFluxes_sfincs(iPhysics,ir), 's', 'Color',colors(iPhysics,:))
        end
    end
    xlabel('dPhiHatdpsiN')
    ylabel('particle flux (sfincs units)')
    title({'Solution for ambipolarity','(Ions solid, electrons dotted)'})
    
    subplot(numRows,numCols,2)
    for iPhysics=1:4
        plot(dPhiHatdpsiNs{ir}{iPhysics}, ionHeatFluxes{ir}{iPhysics}, linespecsIons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('dPhiHatdpsiN')
    ylabel('ion heat flux (sfincs units)')
    
    subplot(numRows,numCols,6)
    for iPhysics=1:4
        plot(dPhiHatdpsiNs{ir}{iPhysics}, electronHeatFluxes{ir}{iPhysics}, linespecsElectrons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('dPhiHatdpsiN')
    ylabel('electron heat flux (sfincs units)')
    
    subplot(numRows,numCols,3)
    for iPhysics=1:4
        plot(dPhiHatdpsiNs{ir}{iPhysics}, ionFlows{ir}{iPhysics}, linespecsIons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('dPhiHatdpsiN')
    ylabel('ion flow (sfincs units)')
    
    subplot(numRows,numCols,7)
    for iPhysics=1:4
        plot(dPhiHatdpsiNs{ir}{iPhysics}, electronFlows{ir}{iPhysics}, linespecsElectrons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('dPhiHatdpsiN')
    ylabel('electron flow (sfincs units)')

    subplot(numRows,numCols,8)
    for iPhysics=1:4
        plot(dPhiHatdpsiNs{ir}{iPhysics}, ionFlows{ir}{iPhysics} - electronFlows{ir}{iPhysics}, linespecsIons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('dPhiHatdpsiN')
    ylabel('ion-electron flow (sfincs units)')
    
    stringForTop = ['Scan of E_r for r/a = ',num2str(r_over_as(ir))];
    annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
               'Interpreter','none','VerticalAlignment','bottom',...
               'FontSize',titleFontSize,'LineStyle','none','String',stringForTop,...
               'FontName','Luxi Sans');

    
    % **********************************************
    % Plot the same results again in Yuriy's units:
    % **********************************************
    
    figure(ir+100)
    clf
    set(gcf,'Color','w')
    numRows = 2;
    numCols = 4;
    
    subplot(numRows,numCols,1)
    for iPhysics=1:4
        plot(ErFactor * dPhiHatdpsiNs{ir}{iPhysics}, particleFluxFactor * ionFluxes{ir}{iPhysics}, linespecsIons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('E_r [kV/m]')
    ylabel('Ion particle flux [10^{20} / (m^2 s)]')
    
    subplot(numRows,numCols,5)
    for iPhysics=1:4
        plot(ErFactor * dPhiHatdpsiNs{ir}{iPhysics}, particleFluxFactor * electronFluxes{ir}{iPhysics}, linespecsElectrons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('E_r [kV/m]')
    ylabel('Electron particle flux [10^{20} / (m^2 s)]')
    
    subplot(numRows,numCols,4)
    for iPhysics=1:4
        if numel(dPhiHatdpsiNs{ir}{iPhysics}) > 1
            [~,permutation] = sort(dPhiHatdpsiNs{ir}{iPhysics});
            plot(ErFactor * dPhiHatdpsiNs{ir}{iPhysics}(permutation), particleFluxFactor * ionFluxes{ir}{iPhysics}(permutation), linespecs{iPhysics}, 'Color',colors(iPhysics,:))
            hold on
            
            dPhiHatdpsiN_fine = linspace(min(dPhiHatdpsiNs{ir}{iPhysics}), max(dPhiHatdpsiNs{ir}{iPhysics}), 100);
            plot(ErFactor * dPhiHatdpsiN_fine, particleFluxFactor * interp1(dPhiHatdpsiNs{ir}{iPhysics}, ionFluxes{ir}{iPhysics}, dPhiHatdpsiN_fine, 'cubic'), '-', 'Color',colors(iPhysics,:))
            
            plot(ErFactor * dPhiHatdpsiNs{ir}{iPhysics}(permutation), particleFluxFactor * electronFluxes{ir}{iPhysics}(permutation), linespecs{iPhysics}, 'Color',colors(iPhysics,:))
            plot(ErFactor * dPhiHatdpsiN_fine, particleFluxFactor * interp1(dPhiHatdpsiNs{ir}{iPhysics}, electronFluxes{ir}{iPhysics}, dPhiHatdpsiN_fine, 'cubic'), ':', 'Color',colors(iPhysics,:))

            plot(ambipolar_Ers(iPhysics,ir), particleFluxes_sfincsDimensional(iPhysics,ir), 's', 'Color',colors(iPhysics,:))
        end
    end
    xlabel('E_r [kV/m]')
    ylabel('Particle flux [10^{20} / (m^2 s)]')
    title({'Solution for ambipolarity','(Ions solid, electrons dotted)'})
    %    title('Solution for ambipolarity')
    
    subplot(numRows,numCols,2)
    for iPhysics=1:4
        plot(ErFactor * dPhiHatdpsiNs{ir}{iPhysics}, heatFluxFactor * ionHeatFluxes{ir}{iPhysics}, linespecsIons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('E_r [kV/m]')
    ylabel('Ion heat flux [MW / m^2]')
    
    subplot(numRows,numCols,6)
    for iPhysics=1:4
        plot(ErFactor * dPhiHatdpsiNs{ir}{iPhysics}, heatFluxFactor * electronHeatFluxes{ir}{iPhysics}, linespecsElectrons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('E_r [kV/m]')
    ylabel('Electron heat flux [MW / m^2]')
    
    subplot(numRows,numCols,3)
    for iPhysics=1:4
        plot(ErFactor * dPhiHatdpsiNs{ir}{iPhysics}, jFactor * ionFlows{ir}{iPhysics}, linespecsIons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('E_r [kV/m]')
    ylabel('Ion contribution to j_{bs} [A / m^2]')
    
    subplot(numRows,numCols,7)
    for iPhysics=1:4
        plot(ErFactor * dPhiHatdpsiNs{ir}{iPhysics}, jFactor * electronFlows{ir}{iPhysics}, linespecsElectrons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('E_r [kV/m]')
    ylabel('Electron contribution to j_{bs} [A / m^2]')

    subplot(numRows,numCols,8)
    for iPhysics=1:4
        plot(ErFactor * dPhiHatdpsiNs{ir}{iPhysics}, jFactor * (ionFlows{ir}{iPhysics} - electronFlows{ir}{iPhysics}), linespecsIons{iPhysics}, 'Color',colors(iPhysics,:))
        hold on
    end
    xlabel('E_r [kV/m]')
    ylabel('Total bootstrap current [A / m^2]')
    
    annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
               'Interpreter','none','VerticalAlignment','bottom',...
               'FontSize',titleFontSize,'LineStyle','none','String',stringForTop,...
               'FontName','Luxi Sans');

    if savePDFs
        export_fig(sprintf('%s_%03d',figureFilename,ir), '-pdf')
    end
end

% **********************************************************
% Now make a figure that summarizes the results
% **********************************************************

figure(100)
clf
set(gcf,'Color','w')
numRows = 3;
numCols = 4;

subplot(numRows,numCols,1)
plot(r_over_as_Yuriy, ne_Yuriy, '-')
xlabel('r/a')
ylabel('n [10^{20} / m^3]')
title('Input: Density profile')

subplot(numRows,numCols,2)
plot(r_over_as_Yuriy, Te_Yuriy, '-')
hold on
plot(r_over_as_Yuriy, Ti_Yuriy, '-r')
xlabel('r/a')
ylabel('T [keV]')
title('Input: Temperature profiles')
legend('T_e','T_i')

subplot(numRows,numCols,3)
plot(r_over_as_Yuriy, Er_Yuriy, '-k')
hold on
for iPhysics = 1:4
    plot(r_over_as, ambipolar_Ers(iPhysics,:), linespecs{iPhysics}, 'Color',colors(iPhysics,:))
end
xlabel('r/a')
ylabel('E_r [kV/m]')
title('Output: E_r')


subplot(numRows,numCols,4)
plot(r_over_as_Yuriy, particleFlux_Yuriy, '-k')
hold on
for iPhysics = 1:4
    plot(r_over_as, particleFluxes_sfincsDimensional(iPhysics,:), linespecs{iPhysics}, 'Color',colors(iPhysics,:))
end
xlabel('r/a')
ylabel('Particle flux [10^{20} / (m^2 s)]')
title('Output: Neoclassical particle flux')

subplot(numRows,numCols,5)
plot(r_over_as_Yuriy, Qi_Yuriy, '-k')
hold on
for iPhysics = 1:4
    plot(r_over_as, ionHeatFluxes_sfincsDimensional(iPhysics,:), linespecs{iPhysics}, 'Color',colors(iPhysics,:))
end
xlabel('r/a')
ylabel('Heat flux [MW / m^2]')
title('Output: Neoclassical ion heat flux')

subplot(numRows,numCols,6)
plot(r_over_as_Yuriy, Qe_Yuriy, '-k')
hold on
for iPhysics = 1:4
    plot(r_over_as, electronHeatFluxes_sfincsDimensional(iPhysics,:), linespecs{iPhysics}, 'Color',colors(iPhysics,:))
end
xlabel('r/a')
ylabel('Heat flux [MW / m^2]')
title('Output: Neoclassical electron heat flux')

subplot(numRows,numCols,7)
plot(r_over_as_Yuriy, ji_Yuriy, '-k')
hold on
for iPhysics = 1:4
    plot(r_over_as, ji_sfincsDimensional(iPhysics,:), linespecs{iPhysics}, 'Color',colors(iPhysics,:))
end
xlabel('r/a')
ylabel('<j_i dot B> / B_0    [A / m^2]')
title({'Output: Ion component','of bootstrap current'})

subplot(numRows,numCols,8)
plot(r_over_as_Yuriy, je_Yuriy, '-k')
hold on
for iPhysics = 1:4
    plot(r_over_as, je_sfincsDimensional(iPhysics,:), linespecs{iPhysics}, 'Color',colors(iPhysics,:))
end
xlabel('r/a')
ylabel('<j_e dot B> / B_0    [A / m^2]')
title({'Output: Electron component','of bootstrap current'})

subplot(numRows,numCols,9)
plot(r_over_as_Yuriy, ji_Yuriy + je_Yuriy, '-k')
hold on
for iPhysics = 1:4
    plot(r_over_as, ji_sfincsDimensional(iPhysics,:) + je_sfincsDimensional(iPhysics,:), linespecs{iPhysics}, 'Color',colors(iPhysics,:))
end
xlabel('r/a')
ylabel('<j dot B> / B_0   [A / m^2]')
title('Output: Total bootstrap current')

legendHandle = legend('DKES calculations by Turkin',...
                      'SFINCS: Fokker-Planck, full trajectories',...
                      'SFINCS: Fokker-Planck, monoenergetic trajectories',...
                      'SFINCS: Pitch-angle scattering, full trajectories',...
                      'SFINCS: Pitch-angle scattering, monoenergetic trajectories');
set(legendHandle,'Interpreter','none',...
                 'Position',[0.348470052083333 0.175821456848773 0.3076171875 0.0964964370546318]);

stringForTop = 'Comparison of SFINCS and momentum-corrected DKES for W7-X';

annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
           'Interpreter','none','VerticalAlignment','bottom',...
           'FontSize',titleFontSize,'LineStyle','none','String',stringForTop,...
           'FontName','Luxi Sans');

temp=dbstack;
stringForBottom = sprintf('Plotted using %s, run in %s',temp(1).file, pwd);

annotation('textbox',[0 0 1 .04],'HorizontalAlignment','center',...
           'Interpreter','none','VerticalAlignment','top',...
           'FontSize',12,'LineStyle','none','String', ...
           stringForBottom);

if savePDFs
    export_fig([figureFilename,'_000'], '-pdf')
end

B0s

save('processedSfincsScanResults')