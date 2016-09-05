function PlotPhi1_PEST_Coordinates()


addpath ~/sfincsProjectsAndTools/tools/Hakan/BoozerFilesAndGeom
addpath ~/sfincsProjectsAndTools/tools/Hakan/version3Scans
% addpath /hydra/u/almo/Multispecies_benchmark/VMEC_data/W7X_ECRH_benchmark
addpath /draco/u/almo/Phi1/LHD/lhd2_A_III/Input

%sfx=h5load('/hydra/u/almo/Multispecies_benchmark/Phi1/Convergence_Scan/Phi1_multispecies_benchmark_baseline3/baseCase/sfincsOutput.h5');
sfx=h5load('/draco/u/almo/Phi1/LHD/lhd2_A_III/runs/lhd2_A_III_EUTERPE_s0p1718750_5/Nx13/sfincsOutput.h5');
%woutfile='wout_w7x.1000_1000_1000_1000_+0000_+0000.01.00jh_l.nc';
woutfile='wout_lhd2.nc'; 

[Pest,Vmec]=makePESTfromVmec(woutfile,double(sfx.psiN),double(sfx.Ntheta),double(sfx.Nzeta),-1);

Vmec.Phi1=squeeze(double(sfx.Phi1Hat(:,:,end)))';

Pest.Phi1=interp2_cyclic(Vmec.vmecu,Vmec.vmecw,Vmec.Phi1,Pest.vmecu,Pest.vmecw,double(sfx.NPeriods));

figure(1)
surf(Pest.pzeta,Pest.ptheta,Pest.Phi1)
shading flat; view(0,90);axis([0,2*pi/5,0,2*pi]);colorbar
xlabel('\zeta_{PEST}');ylabel('\theta_{PEST}')
title('\Phi_1')

figure(2)
surf(Vmec.vmecw,Vmec.vmecu,Vmec.Phi1)
shading flat; view(0,90);axis([0,2*pi/5,0,2*pi]);colorbar
xlabel('v');ylabel('u')
title('\Phi_1') 
