function [runs,miss]=plot_nhatscan(directory)

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
  [runs,miss]=getresults(directory,'nHats');
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


[tauP,tauE]=calcAlltau(runs,151,55);

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
Geom=readBoozerfile(nameliststr(startind:stopind));


rnorm=runs.rN(1);
rind=findnearest(Geom.rnorm,rnorm);
Ntheta=151;
Nzeta=55;%Ntheta;
[Ham,Booz]=makeHamada(Geom,rind,Ntheta,Nzeta);
%FSAgpsipsi=Booz.FSAgpsipsi



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate D11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D11=-Flux_psi/Booz.FSAgpsipsi./(ni20*1e20)./A1;
rLi=vTi./(e*Z.*B00./mi);
Diip=abs(pi./iota).*vTi/16/Geom.majorradius.*rLi.^2;

%nuPrime_appr=nu_ii./vTi./abs(iota).*Geom.majorradius

ind=findnearest(ni20,0.3046);


log_Lambda_i=17.5257;
lc_i =  ( 3/( 4*sqrt(pi)) )*(TikeV*e*1e3).^2./(ni20*1e20.*(Z.^4*e^4)*log_Lambda_i)*(4*pi*eps0)^2;
nuStarAM = 2./abs(iota).*Geom.R00(rind) ./ lc_i;
%nuStarAM=0.23084e-2*nuPrimei/nuPrimei(ind)

Diip_AM=2.2466*(TikeV/TikeV(ind)).^1.5;
FSAgpsipsi_AM=0.2032;
D11_AM=-Flux_psi/FSAgpsipsi_AM./(ni20*1e20)./A1;




fig(1)
loglog(nuStarAM,D11./Diip,'b--+',nuStarAM(ind),D11(ind)./Diip(ind),'bo')
%loglog(nuPrimei,D11./Diip,'b--+',nuPrimei(ind),D11(ind)./Diip(ind),'bo')
ylim([1e-4,1])
xlabel('\nu''')

fig(2)
pos=get(gcf,'Position');
pos(4)=150;
set(gcf,'Position',pos);
loglog(nuStarAM,D11_AM./Diip_AM,'r--+',nuStarAM(ind),D11_AM(ind)./Diip_AM(ind),'ro')
set(gca,'FontSize',14)
%ylim([1e-4,1])
%axis([1e-5 4e-1 1e-4 1])
axis([1e-4 4e-1 1e-4 1e-3])
set(gca,'XTick',[1e-5,1e-4,1e-3,1e-2,1e-1,1e0])
xlabel('\nu^*')
ylabel('D_{11} / D_p')
%figure(2);print -depsc D11overDp_rhopol=0.5_rip.eps
%figure(2);print -depsc D11overDp_rhopol=0.5_rmp90.eps


fig(3)
loglog(nuStarAM,D11,'b--+',nuStarAM(ind),D11(ind),'bo',...
      nuStarAM,D11(ind)*(nuStarAM/nuStarAM(ind)).^0.5,'k:'),...
      %nuStarAM,D11(ind)*(nuStarAM/nuStarAM(ind)).^1,'g:' )
xlim([1e-5,4e-1])
%axis([1e-5 4e-1 1e-4 1])
set(gca,'XTick',[1e-5,1e-4,1e-3,1e-2,1e-1])
xlabel('\nu^*')
ylabel('D_{11} [m^2/s]')
%legend('D_{11}','\nu^{1/2}')
%figure(3);print -depsc D11_rhotor=0.9.eps

fig(4)
loglog(nuStarAM,abs(tauP),'r-',nuStarAM,abs(tauE),'b-',nuStarAM,abs(tauE+tauP),'k--')
xlim([1e-5,4e-1])
set(gca,'XTick',[1e-5,1e-4,1e-3,1e-2,1e-1])
xlabel('\nu^*')
ylabel('torque [Nm/m^3]')
legend('\tau_P','\tau_E','\tau_E+\tau_P')



if 0
  for rind=1:25
    Ntheta=151;
    Nzeta=55;%Ntheta;
    [Ham,Booz]=makeHamada(Geom,rind,Ntheta,Nzeta);
    FSAgpsipsi(rind)=Booz.FSAgpsipsi;
  end
  
  fig(2)
  plot(Geom.rnorm,FSAgpsipsi)
end

