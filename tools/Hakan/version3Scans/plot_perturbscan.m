function [runs,miss]=plot_perturbscan(directory,perturbfac)

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
eps0=8.8542e-12;
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
nu_n=runs.nu_n';



nu_bar=nu_n*vbar/Rbar;
nu_ii=nu_bar.*ni20./sqrt(runs.mHats(:,ion).*TikeV.^3);

nuPrimei=(G+iota.*I).*nu_ii./vTi./B00;


Flux_psi=runs.particleFlux_vm_psiN(:,ion)*nbar*vbar/Rbar*psiAHat;

A1=(runs.dnHatdpsiN(:,ion)./runs.nHats(:,ion)...
    +(runs.dPhiHatdpsiN'*Phibar).*runs.Zs(:,ion).*e./(runs.THats(:,ion)*Tbar)...
    -3/2*runs.dTHatdpsiN(:,ion)./runs.THats(:,ion))/psiAHat;
A2=runs.dTHatdpsiN(:,ion)./runs.THats(:,ion)/psiAHat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nameliststr=runs.run.input_namelist;
nameliststr=nameliststr(find(nameliststr~=char(10))); %remove eol.
eind=findstr(nameliststr,'equilibriumFile');
startind=eind+15;
while nameliststr(startind)~='"';
  startind=startind+1;
end
startind=startind+1;
stopind=startind;
while nameliststr(stopind)~='"';
  stopind=stopind+1;
end
stopind=stopind-1;
fullequilibriumpath=nameliststr(startind:stopind);
strind=findstr(fullequilibriumpath,'equilibria/');
justtheequilibriumfilename=fullequilibriumpath(strind+11:end);
Geom=readBoozerfile(justtheequilibriumfilename);
%eqnamnet=nameliststr(startind:stopind)

rnorm=runs.rN(1);
rind=findnearest(Geom.rnorm,rnorm);
Ntheta=151;
Nzeta=55;%Ntheta;
[Ham,Booz]=makeHamada(Geom,rind,Ntheta,Nzeta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate D11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if all(A2==0)
  D11=-Flux_psi/Booz.FSAgpsipsi./(ni20*1e20)./A1;
  rLi=vTi./(e*Z.*B00./mi);
  Diip=pi./iota.*vTi/16/Geom.majorradius.*rLi.^2;
  D11overDiip=D11./Diip;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fac=[1/8, 1/4, 1/2, 1, 2, 4]';%Hardcoded here
%fac=[1/16, 1/8, 1/4, 1/2, 1, 2, 4]';%Hardcoded here
fac=perturbfac;

pow12=log(Flux_psi(2)/Flux_psi(1))/log(2)
pow13=log(Flux_psi(3)/Flux_psi(1))/log(4)

Flux_r=Flux_psi/sqrt(Booz.FSAgpsipsi)

fac_Flux_re18=[fac,Flux_r/1e18]

Flux_psi(1)/sqrt(Booz.FSAgpsipsi).*(fac/fac(1)).^2./...
    (Flux_psi/sqrt(Booz.FSAgpsipsi))
if 1
fig(1)
loglog(fac,Flux_psi/sqrt(Booz.FSAgpsipsi),'+--g',...
       fac,Flux_psi(1)/sqrt(Booz.FSAgpsipsi).*(fac/fac(1)).^2,'k-')%,...
%       fac,Flux_psi(1)/sqrt(Booz.FSAgpsipsi).*(fac/fac(1)).^1.8,'b:',...
%       fac,Flux_psi(1)/sqrt(Booz.FSAgpsipsi).*(fac/fac(1)).^2.2,'r:')%,...
%       fac,Flux_psi(1)/sqrt(Booz.FSAgpsipsi)*(fac/fac(1)).^pow,'m--')
set(gca,'FontSize',14)
xlabel('perturbation amplitude')
ylabel('radial particle flux [1/(sm^2)]')
%xlim([0.1,5])
%axis([0.1 4,0.8e16,1e19])
%axis([0.05 4,0.5e16,1e19])
%set(gca,'YTick',[10^15,10^16,10^17,10^18,10^19,10^20])
%set(gca,'XTick',[1/16, 0.125, 0.25, 0.5, 1, 2, 4])
legend('SFINCS','(perturbation amplitude)^2',2)
%legend('calculated','\propto F^{1.8}','\propto F^2','\propto F^{2.2}',2)
%legend('calculated','\propto F^2',...
%       ['\propto F^{',num2str(round(pow*100)/100)],'}',2,'interpreter','latex')
%title(['ripple case,  \rho_{pol}=0.5,  d\Phi/ds = ',num2str(runs.dPhiHatdpsiN(1)),' kV'])
%title('half radius')

fig(2)
flux_rel_asymp=Flux_psi/sqrt(Booz.FSAgpsipsi)./...
       (Flux_psi(1)/sqrt(Booz.FSAgpsipsi).*(fac/fac(1)).^2)
semilogx(fac,Flux_psi/sqrt(Booz.FSAgpsipsi)./...
       (Flux_psi(1)/sqrt(Booz.FSAgpsipsi).*(fac/fac(1)).^2),'k-')%,...
%       fac,Flux_psi(1)/sqrt(Booz.FSAgpsipsi).*(fac/fac(1)).^1.8,'b:',...
%       fac,Flux_psi(1)/sqrt(Booz.FSAgpsipsi).*(fac/fac(1)).^2.2,'r:')%,...
%       fac,Flux_psi(1)/sqrt(Booz.FSAgpsipsi)*(fac/fac(1)).^pow,'m--')
set(gca,'FontSize',14)
xlabel('perturbation amplitude')
ylabel('particle flux relative to asymptote')


end

%figure(1);print -depsc rip_perturbscan_Er0.eps
%figure(1);print -depsc rip_perturbscan_realEr_fortalk.eps
%figure(1);print -depsc rmp90_perturbscan.eps
%figure(2);print -depsc rmp90_perturbscan_relative.eps

if all(A2==0)
  fig(2)
  loglog(fac,D11overDiip,'+--b', fac,D11overDiip(1)*(fac/fac(1)).^2,'r--')
  xlabel('perturbation amplification factor F')
  ylabel('D_{11} / D_p')
  xlim([0.1,5])
  legend('calculated','\propto F^2',2)
  title('ripple case, \rho_{pol}=0.5')
end