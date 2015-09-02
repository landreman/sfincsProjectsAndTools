function plot_radscanresults_step23(step2dir,step3dir)

if step2dir(end-3:end)=='.dat'
  dirlist={};
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

e=1.6022e-19;
mp=1.6726e-27;
nbar=1e20;
Tbar=e*1e3;
Rbar=1;
Bbar=1;
pbar=nbar*Tbar;
ion=2; %index for ions in output data
elec=1;%index for electrons in output data
mi=runs2.mHats(:,ion)*mp;
miHat=runs2.mHats(:,ion);
ni=runs2.nHats(:,ion)*nbar;
ne20=runs2.nHats(:,elec);
ni20=runs2.nHats(:,ion);
TekeV=runs2.THats(:,elec);
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

runs2.FSABFlow(:,ion)
runs3.FSABFlow(:,ion)

FSABflowi=runs2.FSABFlow(:,ion)+runs3.FSABFlow(:,ion)


tau1Direct=-runs2.NTV*pbar; %tau is defined as -NTV.
tauinDirect=-runs3.NTV*pbar;

tau1FromFlux=-runs2.NTVfromFlux*pbar;
tauinFromFlux=-runs3.NTVfromFlux*pbar;


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
    -runs2.dnHatdpsiN(:,ion)./runs2.Zs(:,ion)./runs2.nHats(:,ion).*runs2.THats(:,ion)...
    -alpha*runs2.dPhiHatdpsiN');


dne20dPsiN_step2=runs2.dnHatdpsiN(:,elec);
dni20dPsiN_step2=runs2.dnHatdpsiN(:,ion);
dne20dPsiN_step3=runs3.dnHatdpsiN(:,elec);
dni20dPsiN_step3=runs3.dnHatdpsiN(:,ion);
dni20dPsiN=dni20dPsiN_step2+dni20dPsiN_step3;
dTekeVdPsiN=runs3.dTHatdpsiN(:,elec);
dTikeVdPsiN=runs3.dTHatdpsiN(:,ion);
dPotentialkVdPsiN=runs2.dPhiHatdpsiN';


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

FSAomegai_step2=1./(G+iota.*I).*...
    (vbar*Bbar*runs2.FSABFlow(:,ion)./runs2.nHats(:,ion)...
     -runs2.THats(:,ion)*Tbar.*I./runs2.Zs(:,ion)/e.*...
     (A1istep2+5/2*A2istep2));
FSAomegai_step3=1./(G+iota.*I).*...
    (vbar*Bbar*runs3.FSABFlow(:,ion)./runs3.nHats(:,ion)...
     -runs3.THats(:,ion)*Tbar.*I./runs3.Zs(:,ion)/e.*...
     (A1istep2+5/2*A2istep2));
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

dPotentialkVdPsiN_step23_offs1=-iota./(G+iota.*I)./ni20/alpha.*...
    (I./Z.*dTikeVni20dPsiN);
dPotentialkVdPsiN_step23_offs2=-iota./(G+iota.*I)./ni20/alpha.*...
    ( -psiAHat./Deltahalf.*FSABflowi);
dPotentialkVdPsiN_step23_Erin=G./(G+iota.*I).*dPotentialkVdPsiN;



%%%%%%%%%%%%%%%%%%%%%%%%% Calculate FSAg_phiphi %%%%%%%%%%%%%%%%%%
if 1 %approximately
  FSAg_phiphi=(runs2.GHat./runs2.B0OverBBar)'.^2; 
else %calculated
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
  
  disp('Calculating Hamada coordinates')
  FSAg_phiphi=zeros(size(runs2.rN))';
  for ind=1:length(runs2.rN)
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b%3i/%3i',ind,length(runs2.rN))
    rnorm=runs2.rN(ind);
    rind=findnearest(Geom.rnorm,rnorm);
    [Ham,Booz]=makeHamada(Geom,rind,128,128);
    FSAg_phiphi(ind)=Booz.FSAg_phiphi;
  end
  fprintf(1,'\n')
end
%calcu=FSAg_phiphi
%appr=(runs2.GHat./runs2.B0OverBBar)'.^2

nutDirect = -tau1Direct(:,ion)./mi./ni./FSAg_phiphi./omega_torrot;
omegainDirect = -tauinDirect(:,ion)./tau1Direct(:,ion).*omega_torrot;

nutFromFlux = -tau1FromFlux(:,ion)./mi./ni./FSAg_phiphi./omega_torrot;
omegainFromFlux = -tauinFromFlux(:,ion)./tau1FromFlux(:,ion).*omega_torrot;

nutFromL=-Ti.*G.^2.*L11./mi./FSAg_phiphi./vTi./B00./(G+iota.*I);
omegainFromL=(TikeV*1e3*e)./(Z*e.*iota).*(A1istep3+A2i.*L12./L11);


taueDirect=(tau1Direct(:,elec)+tauinDirect(:,elec))';
tauiDirect=(tau1Direct(:,ion)+tauinDirect(:,ion))';
taueFromFlux=(tau1FromFlux(:,elec)+tauinFromFlux(:,elec))';
tauiFromFlux=(tau1FromFlux(:,ion)+tauinFromFlux(:,ion))';
tau_totDirect=taueDirect+tauiDirect;
tau_totFromFlux=taueFromFlux+tauiFromFlux;

integrDirect=tau_totDirect...
    .*4*pi^2./runs2.FSABHat2.*(runs2.GHat+runs2.iota.*runs2.IHat)*...
    psiAHat*Rbar^3;
integrFromFlux=tau_totFromFlux...
    .*4*pi^2./runs2.FSABHat2.*(runs2.GHat+runs2.iota.*runs2.IHat)*...
    psiAHat*Rbar^3;

s=runs2.rN.^2;
tau_NmDirect=trapz(s,integrDirect)
tau_NmFromFlux=trapz(s,integrFromFlux)


fig(1)
plot(runs2.rN,omega_torrot,'k-',runs2.rN,omegainDirect,'g--',runs2.rN,omegainFromFlux,'b--',runs2.rN,omegainFromL,'r-.')
xlabel('r / a')
ylabel('\omega  [s^{-1}]')
legend('\omega_{meas}','\omega_{in} direct','\omega_{in} from flux','\omega_{in} from L')
title('rotation frequency')

fig(2)
plot(runs2.rN,1./nutDirect,'g--',runs2.rN,1./nutFromFlux,'b-',runs2.rN,1./nutFromL,'r-.')
xlabel('r / a')
ylabel('1 / \nu_t  [s]')
legend('direct','from flux','from L')
title('braking time')

fig(3)
plot(runs2.rN,taueDirect,  runs2.rN,tauiDirect,...
     runs2.rN,taueFromFlux,runs2.rN,tauiFromFlux)
title('torque')
xlabel('r / a')
ylabel('Nm / m^3')
legend('e direct','i direct','e from flux','i from flux')


fig(4)
plot(runs2.rN,dPotentialkVdPsiN_step23,...
     runs2.rN,dPotentialkVdPsiN_step23_omegatorrot,...
     runs2.rN,dPotentialkVdPsiN_step23_offs1,...
     runs2.rN,dPotentialkVdPsiN_step23_offs2,'r--',...
     runs2.rN,dPotentialkVdPsiN_step23_Erin,'-',...
     runs2.rN,dPotentialkVdPsiN,'--')
title('radial electric field')
xlabel('r / a')
ylabel('d\Phi / d\psi_N')
legend('step 2&3: total','2&3: \omega','2&3: dp/dr','2&3: flow','2&3: E_r in','from step 1')

fig(7)
plot(runs2.rN,kappaiFSAB2_step2,runs3.rN,kappaiFSAB2_step3,...
     runs2.rN,kappaiFSAB2_step2+kappaiFSAB2_step3)
title('\kappa_i <B^2>')
xlabel('r / a')
legend('step 2','step 3')
%axis([0 1 0 16e4])

fig(8)
plot(runs2.rN,FSAomegai_step2,'m',...
     runs3.rN,FSAomegai_step3,'c',...
     runs2.rN,FSAomegai_step2+FSAomegai_step3,'r--',...
     runs2.rN,omega_torrot,'k')
title('<\omega_i>')
xlabel('r / a')
legend('step 2','step 3','2 + 3','input')


%figure(1);print -depsc 30835_rip90_omega_in.eps
%figure(2);print -depsc 30835_rip90_nut.eps
%figure(3);print -depsc 30835_rip90_torque.eps

