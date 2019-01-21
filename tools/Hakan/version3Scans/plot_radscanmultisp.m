function [runs,miss]=plot_radscanmultisp(directory)

if directory(end-3:end)=='.dat'
  dirlist={};
  if not(exist(directory,'file'))
    if directory(1)=='/'
      directory=[getenv('SFINCS_HOME'),'/fortran/version3',directory];
    else    
      directory=[getenv('SFINCS_HOME'),'/fortran/version3/',directory];
    end
  end
  fid = fopen(directory);
  tline = fgetl(fid);
  while ischar(tline)
    if tline(1)~='%' && tline(1)~='!'
      dirlist={dirlist{:},tline};
    end
    tline = fgetl(fid);
  end
  fclose(fid);
  [runs,miss]=getresults(dirlist);
else
  [runs,miss]=getresults(directory);
end

if runs.NumElements==0
  error('Nothing found!')
end

e=1.6022e-19;
mp=1.6726e-27;
mbar=mp;
nbar=1e20;
Tbar=e*1e3;
Rbar=1;
Bbar=1;
Tbar=1.6022e-19*1e3;
Phibar=1e3;
pbar=nbar*Tbar;
psiAHat=runs.psiAHat(1);
vbar=sqrt(1e3*e*2/mp);
iota=runs.iota';
G=runs.GHat'*Bbar*Rbar;
I=runs.IHat'*Bbar*Rbar;
iota=runs.iota';
B00=runs.B0OverBBar'*Bbar;

psiAHat=runs.psiAHat(1);
Nspec=size(runs.Zs,2);

ion=find(runs.Zs(1,:)~=-1);
TikeV=runs.THats(:,ion);
vTi=sqrt(TikeV*e*1e3*2/mp./runs.mHats(:,ion));
ni20=runs.nHats(:,ion);
Z=runs.Zs(:,ion);
mi=runs.mHats(:,ion)*mp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fz=14;
addpath ~/sfincs/sfincs/equilibria/
%boozerfile='w7x-lim-op1_1.bc';
boozerfile='w7x-sc1.bc';
Geom=readBoozerfile(boozerfile,0.1,3,3,'StelSym','signcorr2');
a=Geom.minorradiusW7AS;
%dVdsoverNper=interp1(Geom.rnorm,Geom.dVdsoverNper,runs.rN);
%dVdr=dVdsoverNper*Geom.Nperiods*2.*runs.rN/Geom.minorradiusW7AS;
dpsiNdr=2*runs.rN/a;

dVdpsiN=abs(runs.VPrimeHat*Rbar/Bbar.*runs.psiAHat);
%dVdpsiN=runs.VPrimeHat*Rbar/Bbar.*runs.psiAHat;

%runs.VPrimeHat*Rbar/Bbar.*runs.psiAHat

if runs.includePhi1
   partFluxpers= runs.particleFlux_vd_psiN*vbar*nbar/Rbar.*(ones(Nspec,1)*dVdpsiN)';
   heatFluxMW  = runs.heatFlux_vd_psiN*vbar^3*nbar*mbar/Rbar.*(ones(Nspec,1)*dVdpsiN)'/1e6;
else
   partFluxpers= runs.particleFlux_vm_psiN*vbar*nbar/Rbar.*(ones(Nspec,1)*dVdpsiN)';
   partFluxpersm2=runs.particleFlux_vm_psiN*vbar*nbar/Rbar.*(ones(Nspec,1)*(1./dpsiNdr))';
   heatFluxMW  = runs.heatFlux_vm_psiN*vbar^3*nbar*mbar/Rbar.*(ones(Nspec,1)*dVdpsiN)'/1e6;
end
if Nspec==3
  partFluxpers_spec1and2=partFluxpers(:,1:2)
  partFluxpers_spec3=partFluxpers(:,3)
  partFluxpersm2_spec1and2=partFluxpersm2(:,1:2)
  partFluxpersm2_spec3=partFluxpersm2(:,3)
else  
  partFluxpers_specnoeand1=[NaN*ones(runs.NumElements,1),partFluxpers(:,1)]
  partFluxpers_spec2=partFluxpers(:,2)  
end
if 0%isfield(runs,'classicalParticleFlux_psiN')
  classpartFluxpers=runs.classicalParticleFlux_psiN*vbar*nbar/Rbar.*(ones(Nspec,1)*dVdpsiN)';
  classpartFluxpers_spec1and2=classpartFluxpers(:,1:2);
  classpartFluxpers_spec3=classpartFluxpers(:,3);
  %cl2ncl_ratio_spec1and2=classpartFluxpers_spec1and2./partFluxpers_spec1and2
  %cl2ncl_ratio_spec3=classpartFluxpers_spec3./partFluxpers_spec3
end
%partFluxpers_f0= runs.particleFlux_vm0_psiN*vbar*nbar/Rbar.*(ones(Nspec,1)*dVdpsiN)';
%partFluxpers_df=partFluxpers-partFluxpers_f0;
%heatFluxMW=runs.heatFlux_vm_psiN*vbar^3*nbar*mbar/Rbar.*(ones(Nspec,1)*dVdpsiN)'/1e6;
%heatFluxMW_f0=runs.heatFlux_vm0_psiN*vbar^3*nbar*mbar/Rbar.*(ones(Nspec,1)*dVdpsiN)'/1e6;
%heatFluxMW_df=heatFluxMW-heatFluxMW_f0;
flowpersm2=runs.FSABFlow*vbar*nbar*Bbar./B00;

runs.rN(find(runs.rN<0.1))=NaN;
sumZpartFluxpers=runs.Zs(1,:)*partFluxpers';


if 1
  fig(1)
  if 1
    plot(runs.rN,runs.Zs.*partFluxpers,'-',...
         runs.rN,sumZpartFluxpers,'k--')%,runs.rN,partFluxpers_df,'--')
    ylabel('Z \Gamma [1/s]')
  else
    plot(runs.rN,partFluxpers(:,1),...
         runs.rN,partFluxpers(:,2)*100)
    legend('\Gamma_1','\Gamma_2*100')
    ylabel('\Gamma [1/s]')
  end
  xlabel('r/a')


  %axis([0,1,-2e20,3e21])


  fig(2)
  plot(runs.rN,heatFluxMW,'-')%,runs.rN,heatFluxMW_df,'--')
  xlabel('r/a')
  ylabel('Q [MW]')
  legend('Q_1','Q_2')

  fig(3)
  plot(runs.rN,runs.FSABFlow)
  xlabel('r/a')
  legend('spec 1','spec 2')
  ylabel('flow [1/(sm^2)]')
end