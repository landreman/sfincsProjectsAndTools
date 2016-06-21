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
Nspec=size(runs.NTV,2);
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
%addpath ~/sfincs/sfincs/equilibria/
%boozerfile='w7x-lim-op1_1.bc';
%Geom=readBoozerfile(boozerfile);
%a=Geom.minorradiusW7AS;
%dVdsoverNper=interp1(Geom.rnorm,Geom.dVdsoverNper,runs.rN);
%dVdr=dVdsoverNper*Geom.Nperiods*2.*runs.rN/Geom.minorradiusW7AS;
dVdpsiN=runs.VPrimeHat*Rbar/Bbar.*runs.psiAHat;

partFluxpers= runs.particleFlux_vm_psiN*vbar*nbar/Rbar.*([1;1]*dVdpsiN)';
heatFluxMW=runs.heatFlux_vm_psiN*vbar^3*nbar*mbar/Rbar.*([1;1]*dVdpsiN)'/1e6;

fig(1)
plot(runs.rN,partFluxpers)
xlabel('r/a')
ylabel('\Gamma [1/s]')
legend('e','i')

fig(2)
plot(runs.rN,heatFluxMW)
xlabel('r/a')
ylabel('Q [MW]')
legend('e','i')

fig(3)
plot(runs.rN,runs.FSABFlow)