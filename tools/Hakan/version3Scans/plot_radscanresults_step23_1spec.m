function [runs2,runs3]=plot_radscanresults_step23_1spec(step2dir,step3dir)

if step2dir(end-3:end)=='.dat'
  dirlist={};
  if not(exist(step2dir,'file'))
    if step2dir(1)=='/'
      step2dir=[getenv('SFINCS_HOME'),'/fortran/version3',step2dir];
    else    
      step2dir=[getenv('SFINCS_HOME'),'/fortran/version3/',step2dir];
    end
  end
  fid = fopen(step2dir);
  tline = fgetl(fid);
  while ischar(tline)
    if tline(1)~='%' && tline(1)~='!'
      dirlist={dirlist{:},tline};
    end
    tline = fgetl(fid);
  end
  fclose(fid);
  [runs2,miss2]=getresults(dirlist);
else
  [runs2,miss2]=getresults(step2dir);
end

if step3dir(end-3:end)=='.dat'
  dirlist={};
  if not(exist(step3dir,'file'))
    if step3dir(1)=='/'
      step3dir=[getenv('SFINCS_HOME'),'/fortran/version3',step3dir];
    else    
      step3dir=[getenv('SFINCS_HOME'),'/fortran/version3/',step3dir];
    end
  end
  fid = fopen(step3dir);
  tline = fgetl(fid);
  while ischar(tline)
    if length(tline)>0
      if tline(1)~='%' && tline(1)~='!'
        dirlist={dirlist{:},tline};
      end
    end
    tline = fgetl(fid);
  end
  fclose(fid);
  [runs3,miss3]=getresults(dirlist);
else
  [runs3,miss3]=getresults(step3dir);
end

[tauP2,tauE2]=calcAlltau(runs2);
[tauP3,tauE3]=calcAlltau(runs3);
%[tauP2,tauE2]=calcAlltau(runs2,151,35);
%[tauP3,tauE3]=calcAlltau(runs3,151,35);

tauE=tauE2'+tauE3'
%tauE2
%tauE3
tauP=tauP2'+tauP3'

e=1.6022e-19;
mp=1.6726e-27;
nbar=1e20;
Tbar=e*1e3;
Rbar=1;
Bbar=1;
pbar=nbar*Tbar;
ion=1;
mi=runs2.mHats(:,ion)*mp;
if any(diff(runs2.mHats(:,ion)))
  error('Different masses')
else
  mi=mi(1);
end
miHat=runs2.mHats(:,ion);
ni=runs2.nHats(:,ion)*nbar;
ni20=runs2.nHats(:,ion);
TikeV=runs2.THats(:,ion);
Ti=TikeV*1e3*e;
vTi=sqrt(Ti*2./mi);
psiAHat=runs2.psiAHat(1);
alpha=1;
Z=runs2.Zs(:,ion);
Deltahalf=runs2.Delta' / 2;
vbar=sqrt(1e3*e*2/mp);
iota=runs2.iota';
G=runs2.GHat';
I=runs2.IHat';
B00=runs2.B0OverBBar';


%runs2.FSABFlow(:,ion)
%runs3.FSABFlow(:,ion)

FSABflowi=runs2.FSABFlow(:,ion)+runs3.FSABFlow(:,ion);


tau1Direct=-runs2.NTV*pbar; %tau is defined as -NTV.
tauinDirect=-runs3.NTV*pbar;

tau1FromFlux=-runs2.NTVfromFlux*pbar;
tauinFromFlux=-runs3.NTVfromFlux*pbar;

tau1FromFlux_over_tau1Direct=(tau1FromFlux./tau1Direct)';
tauinFromFlux_over_tauinDirect=(tauinFromFlux./tauinDirect)';

if runs2.NumElements~=runs3.NumElements
  error('The number of finished runs in step 2 and 3 are different!')
end
if any(runs2.rN~=runs3.rN)
  error('The radii in step 2 and 3 do not match')
end

%We can extract omega_meas from the input data of step 2.
%the density gradient was calculated using
%dni20dPsiN_step2= -Z  *ni20./TikeV.*...
%      (alpha*dPotentialkVdPsiN + iota./Deltahalf.*omega_torrot/vbar);
%dni20dPsiN_step2= -Z  *ni20./TikeV.*...
%      (alpha*dPotentialkVdPsiN + iota./Deltahalf.*omega_torrot/vbar*psiAHat);

omega_torrot=vbar.*Deltahalf./iota/psiAHat .* ( ...
    -runs2.dnHatdpsiN(:,ion)./runs2.Zs(:,ion)...
    ./runs2.nHats(:,ion).*runs2.THats(:,ion)...
    -alpha*runs2.dPhiHatdpsiN');


dni20dPsiN_step2=runs2.dnHatdpsiN(:,ion);
dni20dPsiN_step3=runs3.dnHatdpsiN(:,ion);
dni20dPsiN=dni20dPsiN_step2+dni20dPsiN_step3;
dTikeVdPsiN=runs3.dTHatdpsiN(:,ion);
dPotentialkVdPsiN=runs2.dPhiHatdpsiN'; %=runs3.dPhiHatdpsiN', The LHS value


% Calculate thermodyn. forces
A1istep2=(dni20dPsiN_step2./ni20 ...
          +Z./TikeV.*alpha.*dPotentialkVdPsiN)/psiAHat;
A2istep2=zeros(size(A1istep2));

A1istep3=(dni20dPsiN_step3./ni20 ...
          +Z./TikeV.*alpha.*dPotentialkVdPsiN ...
          -3/2*dTikeVdPsiN./TikeV)/psiAHat;
A2istep3=(dTikeVdPsiN./TikeV)/psiAHat;

A1i=A1istep2+A1istep3;
A2i=A2istep2+A2istep3;

%pretend the electrons have no influence:
factor=iota.^2.*B00.*(G+iota.*I)./G.^2./TikeV.^2.*Z.^2 ...
    .*sqrt(TikeV./miHat)./Deltahalf.^2./ni20;
L11=runs2.particleFlux_vm_psiN(:,ion)*psiAHat.*factor./A1istep2;
L12=(runs3.particleFlux_vm_psiN(:,ion)*psiAHat.*factor-L11.*A1istep3)./A2istep3;

L12overL11_2p85inSqrtnuRegime=L12'./L11'
L31=runs2.FSABFlow(:,ion)./runs2.nHats(:,ion)*vbar*Bbar^2*Rbar.*...
    iota.*Z.*e./Ti./G./A1istep2;
L31_should_be_minus_one=L31'

kappaiFSAB2_step2=vbar*Bbar*runs2.FSABFlow(:,ion)./runs2.nHats(:,ion)+...
      runs2.THats(:,ion)*Tbar.*G./runs2.Zs(:,ion)/e./iota.*...
      (A1istep2+5/2*A2istep2);
kappaiFSAB2_step3=vbar*Bbar*runs3.FSABFlow(:,ion)./runs3.nHats(:,ion)+...
      runs3.THats(:,ion)*Tbar.*G./runs3.Zs(:,ion)/e./iota.*...
      (A1istep3+5/2*A2istep3);

if 1
  FSAomegai_step2=kappaiFSAB2_step2./(G+iota.*I) ...
      -runs2.THats(:,ion)*Tbar./iota./runs2.Zs(:,ion)/e.*...
       (A1istep2+5/2*A2istep2);
  FSAomegai_step3=kappaiFSAB2_step3./(G+iota.*I) ...
      -runs3.THats(:,ion)*Tbar./iota./runs3.Zs(:,ion)/e.*...
      (A1istep3+5/2*A2istep3);
else
  FSAomegai_step2=1./(G+iota.*I).*...
      (vbar*Bbar*runs2.FSABFlow(:,ion)./runs2.nHats(:,ion)...
       -runs2.THats(:,ion)*Tbar.*I./runs2.Zs(:,ion)/e.*...
       (A1istep2+5/2*A2istep2));
  FSAomegai_step3=1./(G+iota.*I).*...
      (vbar*Bbar*runs3.FSABFlow(:,ion)./runs3.nHats(:,ion)...
       -runs3.THats(:,ion)*Tbar.*I./runs3.Zs(:,ion)/e.*...
       (A1istep3+5/2*A2istep3));
end
FSAomegai=FSAomegai_step2+FSAomegai_step3;





%Calculate dPotentialkVdPsiN from the added result FSABflowi of step 2  and 3
%This only works for ErFromRot, where omega_torrotLHS=omega_torrot
dTikeVni20dPsiN=ni20.*dTikeVdPsiN+TikeV.*dni20dPsiN;
if 1
  dPotentialkVdPsiN_step23= ...
    -psiAHat/alpha./Deltahalf/vbar.*iota.*FSAomegai...
      -iota./(G+iota.*I)./ni20/alpha.*...
      (I./Z.*dTikeVni20dPsiN - psiAHat./Deltahalf.*FSABflowi) ...
    +G./(G+iota.*I).*dPotentialkVdPsiN;
  dPotentialkVdPsiN_step23_omegatorrot=-psiAHat/alpha./Deltahalf/vbar.*iota.*FSAomegai;
else
    dPotentialkVdPsiN_step23= ...
    -psiAHat/alpha./Deltahalf/vbar.*iota.*omega_torrot...
      -iota./(G+iota.*I)./ni20/alpha.*...
      (I./Z.*dTikeVni20dPsiN - psiAHat./Deltahalf.*FSABflowi) ...
    +G./(G+iota.*I).*dPotentialkVdPsiN;
      dPotentialkVdPsiN_step23_omegatorrot=-psiAHat/alpha./Deltahalf/vbar.*iota.*omega_torrot;
end
%Calculate some terms for plotting
dPotentialkVdPsiN_step23_omegatorrot=-psiAHat/alpha./Deltahalf/vbar.*iota.*omega_torrot;
dPotentialkVdPsiN_step23_offs1=-iota./(G+iota.*I)./ni20/alpha.*...
    (I./Z.*dTikeVni20dPsiN);
dPotentialkVdPsiN_step23_offs2=-iota./(G+iota.*I)./ni20/alpha.*...
    ( -psiAHat./Deltahalf.*FSABflowi);
dPotentialkVdPsiN_step23_Erin=G./(G+iota.*I).*dPotentialkVdPsiN;



%%%%%%%%%%%%%%%%%%%%%%%%% Calculate FSAg_phiphi %%%%%%%%%%%%%%%%%%
nameliststr=runs2.run.input_namelist;
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

if 0 %approximately
  FSAg_phiphi=(runs2.GHat./runs2.B0OverBBar)'.^2; 
else %calculated
  
  disp('Calculating Hamada coordinates')
  FSAg_phiphi=zeros(size(runs2.rN))';
  FSAgpsipsi=zeros(size(runs2.rN))';
  for ind=1:length(runs2.rN)
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b%3i/%3i',ind,length(runs2.rN))
    rnorm=runs2.rN(ind);
    rind=findnearest(Geom.rnorm,rnorm);
    [Ham,Booz]=makeHamada(Geom,rind,85,35,'forceSize');
    FSAg_phiphi(ind)=Booz.FSAg_phiphi;
    FSAgpsipsi(ind)=Booz.FSAgpsipsi;
  end
  fprintf(1,'\n')
end

nutDirect = -tau1Direct(:,ion)./mi./ni./FSAg_phiphi./omega_torrot;
omegainDirect = -tauinDirect(:,ion)./tau1Direct(:,ion).*omega_torrot;

nutFromFlux = -tau1FromFlux(:,ion)./mi./ni./FSAg_phiphi./omega_torrot;
omegainFromFlux = -tauinFromFlux(:,ion)./tau1FromFlux(:,ion).*omega_torrot;

nutFromL=-Ti.*G.^2.*L11./mi./FSAg_phiphi./vTi./B00./(G+iota.*I);
omegainFromL=(TikeV*1e3*e)./(Z*e.*iota).*(A1istep3+A2i.*L12./L11);



tauiDirect=(tau1Direct(:,ion)+tauinDirect(:,ion))';

tauiFromFlux=(tau1FromFlux(:,ion)+tauinFromFlux(:,ion))';
tau_totDirect=tauiDirect;
tau_totFromFlux=tauiFromFlux;

tau_fromPresAnis=tauP2+tauP3

integrDirect=tau_totDirect...
    .*4*pi^2./runs2.FSABHat2.*(runs2.GHat+runs2.iota.*runs2.IHat)*...
    psiAHat*Rbar^3;
integrFromFlux=tau_totFromFlux...
    .*4*pi^2./runs2.FSABHat2.*(runs2.GHat+runs2.iota.*runs2.IHat)*...
    psiAHat*Rbar^3;

s=runs2.rN.^2;
tau_NmDirect=trapz(s,integrDirect)
tau_NmFromFlux=trapz(s,integrFromFlux)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Approximate analytical calculations in the ripple plateau
% for the AUG equilibrium rip_n16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
make_rippleplateau=1;
if make_rippleplateau
  
  FSAB2=zeros(length(runs2.rN),1);
  FSAg_phiphi=zeros(length(runs2.rN),1);
  integr=zeros(length(runs2.rN),1);
  FSAdelta2overB=zeros(length(runs2.rN),1);
  
  for ind=1:length(runs2.rN)
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b%3i/%3i',ind,length(runs2.rN))
    rnorm=runs2.rN(ind);
    rind=findnearest(Geom.rnorm,rnorm);
    Ntheta=85;
    Nzeta=35;%Ntheta;
    [Ham,Booz]=makeHamada(Geom,rind,Ntheta,Nzeta,'forceSize');
    FSAg_phiphi(ind)=Booz.FSAg_phiphi;
    FSAB2(ind)=Booz.FSAB2;
    Btormean=mean(Booz.B,2)*ones(1,Nzeta);
    deltatw=(Booz.B-Btormean)./Btormean;
    integr(ind)=4*pi^2/(Ntheta*Nzeta)*sum(sum(deltatw.^2./Btormean.^3))*2;
    FSAdelta2overB(ind)=2*sum(sum(deltatw.^2./Btormean.^3))/sum(sum(1./Booz.B.^2));
  end
  fprintf(1,'\n')
  
  
  plat.L11=Geom.Nperiods*FSAB2.*integr/4/pi^2.*B00*sqrt(pi)...
           .*(G+iota.*I).^2./G.^2;
  plat.L12=plat.L11*3;
  %FSAg_phiphi_appr=R00.^2;
  %ripple plateau predictions
  plat.nut=Ti.*G.^2.*plat.L11/(mi)./FSAg_phiphi./vTi./B00./(G+iota.*I);
  plat.omegain=omega_torrot +Ti./(Z*e.*iota).*(A1i+A2i.*plat.L12./plat.L11);
  %omegain=Ti./(Z*e.*iota).*(A1istep3+A2i.*L12./L11); %equal to the one above
  plat.tauhat=plat.nut*(mi).*FSAg_phiphi./(Z*e.*iota).*(A1i+A2i.*plat.L12./plat.L11);
  plat.tau=(TikeV*1e3*e).*(ni20*1e20).*plat.tauhat;
  
  
  plat.tau2=Geom.Nperiods*mi^2.*vTi.^3.*(ni20*1e20)./(Z*e)./iota*...
            sqrt(pi)/4.*(G+iota.*I).*FSAdelta2overB.*(A1i+A2i*3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Andreas type values for D11  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Flux_psi=runs2.particleFlux_vm_psiN(:,ion)*nbar*vbar/Rbar*psiAHat;
D11=-Flux_psi./FSAgpsipsi./(ni20*1e20)./A1istep2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rNplot=runs2.rN;
rNplot(end)=NaN;
fz=14;

fig(1)
if make_rippleplateau
  plot(rNplot,omega_torrot,'k-',...
       rNplot,omegainDirect,'r--+',...
       rNplot,omegainFromFlux,'b--d',...
       rNplot,plat.omegain,'g--o')
else
    plot(rNplot,omega_torrot,'k-',...
       rNplot,omegainDirect,'r--+',...
       rNplot,omegainFromFlux,'b--d')
end
set(gca,'FontSize',fz)
xlabel('\rho_{tor}')
ylabel('\omega   [ rad / s ]')
legend('\omega_{meas}','\omega_{in} from anisotropy','\omega_{in} from flux','analytic approx.')%'\omega_{in} from L',4)
title('rotation frequency')

fig(2)
if make_rippleplateau
plot(rNplot,1./nutDirect,'r--+',...
     rNplot,1./nutFromFlux,'b--d',...
     rNplot,1./plat.nut,'g--o')
else
  plot(rNplot,1./nutDirect,'r--+',...
     rNplot,1./nutFromFlux,'b--d')
end  
set(gca,'FontSize',fz)
xlabel('\rho_{tor}')
ylabel('1 / \nu_t  [s]')
legend('from anisotropy','from flux','analytic approx.')
%legend('from anisotropy','from flux','from L')
title('braking time')

fig(3)
if make_rippleplateau
    plot(rNplot,tauiDirect,'r--+',...
       rNplot,tauiDirect+tauE,'r--o',...
       rNplot,tauiFromFlux,'b--d',...
       rNplot,plat.tau2,'g--o')%,...
else
    plot(rNplot,tauiDirect,'r--+',...
         rNplot,tauiDirect+tauE,'r--o',...
         rNplot,tauiFromFlux,'b--d')
end
set(gca,'FontSize',fz)
title('torque')
xlabel('\rho_{tor}')
ylabel('Nm / m^3')
legend('from anisotropy','from anisotropy and Er','from flux','analytic approx.',2)
ax=axis;
ax(3)=0;
axis(ax);
%plat.tau2./tauiFromFlux'


fig(4)
plot(rNplot,-dPotentialkVdPsiN_step23,...
     rNplot,-dPotentialkVdPsiN_step23_omegatorrot,...
     rNplot,-dPotentialkVdPsiN_step23_offs1,...
     rNplot,-dPotentialkVdPsiN_step23_offs2,'r--',...
     rNplot,-dPotentialkVdPsiN_step23_Erin,'-',...
     rNplot,-dPotentialkVdPsiN,'--')
set(gca,'FontSize',fz)
title('radial electric field')
xlabel('\rho_{tor}')
ylabel('- d\Phi / d\psi_N')
legend('step 2&3: total','2&3: \omega','2&3: dp/dr','2&3: flow','2&3: E_r in','from step 1')

fig(7)
plot(rNplot,kappaiFSAB2_step2,runs3.rN,kappaiFSAB2_step3,...
     rNplot,kappaiFSAB2_step2+kappaiFSAB2_step3)
set(gca,'FontSize',fz)
title('\kappa_i <B^2>')
xlabel('\rho_{tor}')
legend('step 2','step 3','2 + 3')
%axis([0 1 0 16e4])

fig(8)
plot(rNplot,FSAomegai_step2,'m',...
     rNplot,FSAomegai_step3,'c',...
     rNplot,FSAomegai_step2+FSAomegai_step3,'r--',...
     rNplot,omega_torrot,'k')
set(gca,'FontSize',fz)
title('toroidal flow <\omega_i>')
xlabel('\rho_{tor}')
legend('step 2','step 3','2 + 3','input')



%figure(1);print -depsc 30835_rip90_omega_in.eps
%figure(2);print -depsc 30835_rip90_nut.eps
%figure(3);print -depsc 30835_rip90_torque.eps

%figure(1);print -depsc 30835_rmp90_omega_in.eps
%figure(2);print -depsc 30835_rmp90_nut.eps
%figure(3);print -depsc 30835_rmp90_torque.eps

if 0
  fig(9)
  plot(rNplot,omega_torrot,...
       rNplot,-A1istep2.*(TikeV*1e3)./Z./iota)
  set(gca,'FontSize',fz)
end

fig(20)
semilogy(rNplot,D11)
xlabel('\rho_{tor}')
ylabel('D_{11} [m^2/s]')

if make_rippleplateau
  fig(10)
  plot(rNplot,L11,rNplot,-plat.L11)
  set(gca,'FontSize',fz)
  xlabel('\rho_{tor}')
  title('L_{11}')
end