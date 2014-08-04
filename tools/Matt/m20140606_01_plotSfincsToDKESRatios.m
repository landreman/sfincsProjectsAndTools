filename='processedSfincsScanResults';
load(filename)

xRange=[0,1];
yRange = [0.95,1.2];

figure(200)
clf
set(gcf,'Color','w')
numRows = 2;
numCols = 2;

subplot(numRows,numCols,1)
particleFlux_Yuriy_interpolated = interp1(r_over_as_Yuriy, particleFlux_Yuriy, r_over_as);
for iPhysics = 1:4
    plot(r_over_as, particleFluxes_sfincsDimensional(iPhysics,:) ./ particleFlux_Yuriy_interpolated,...
         linespecs{iPhysics}, 'Color',colors(iPhysics,:))
    hold on
end
plot([0,1],[1,1],':k')
xlabel('r/a')
title('Particle flux: SFINCS / Yuriy')
xlim(xRange)
ylim(yRange)

subplot(numRows,numCols,2)
Qi_Yuriy_interpolated = interp1(r_over_as_Yuriy, Qi_Yuriy, r_over_as);
for iPhysics = 1:4
    plot(r_over_as, ionHeatFluxes_sfincsDimensional(iPhysics,:) ./ Qi_Yuriy_interpolated,...
         linespecs{iPhysics}, 'Color',colors(iPhysics,:))
    hold on
end
plot([0,1],[1,1],':k')
xlabel('r/a')
title('Ion heat flux: SFINCS / Yuriy')
xlim(xRange)
ylim(yRange)

subplot(numRows,numCols,3)
Qe_Yuriy_interpolated = interp1(r_over_as_Yuriy, Qe_Yuriy, r_over_as);
for iPhysics = 1:4
    plot(r_over_as, electronHeatFluxes_sfincsDimensional(iPhysics,:) ./ Qe_Yuriy_interpolated,...
         linespecs{iPhysics}, 'Color',colors(iPhysics,:))
    hold on
end
xlabel('r/a')
title('Electron heat flux: SFINCS / Yuriy')
xlim(xRange)
ylim(yRange)

legendHandle = legend('SFINCS: Fokker-Planck, full trajectories',...
                      'SFINCS: Fokker-Planck, monoenergetic trajectories',...
                      'SFINCS: Pitch-angle scattering, full trajectories',...
                      'SFINCS: Pitch-angle scattering, monoenergetic trajectories');
set(legendHandle,'Interpreter','none',...
                 'Position',[0.6 0.175821456848773 0.3076171875 0.0964964370546318]);
plot([0,1],[1,1],':k')
