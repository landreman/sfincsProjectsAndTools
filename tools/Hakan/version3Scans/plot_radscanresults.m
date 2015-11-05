function [runs,miss]=plot_radscanresults(directory)

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

%runs

[tauP,tauE]=calcAlltau(runs,151,35);


%the_first_two_should_be_equal=[tauP,-pbar*runs.NTV(:,1),tauE]

fig(1)
if Nspec==2
  plot(runs.rN,-pbar*runs.NTV(:,1),'r',runs.rN,-pbar*runs.NTV(:,2),'b',...
       runs.rN,-pbar*runs.NTVfromFlux(:,1),'r--',runs.rN,-pbar*runs.NTVfromFlux(:,2),'b--')
  legend('\tau e','\tau i','\tau fromFlux e','\tau fromFlux i')
elseif Nspec==1
    plot(runs.rN,-pbar*runs.NTV(:,1),'b',...
       runs.rN,-pbar*runs.NTVfromFlux(:,1),'b--')
  legend('\tau','\tau fromFlux')
end
title('\tau ( = -NTV )')
xlabel('r / a')


fig(2)

s=runs.rN.^2;

if Nspec==2
  NTVtot=pbar*(runs.NTV(:,1)+runs.NTV(:,2));
  NTVfromFluxtot=pbar*(runs.NTVfromFlux(:,1)+runs.NTVfromFlux(:,2));
  plot(runs.rN,-NTVtot,'g',...
       runs.rN,-NTVfromFluxtot,'g--')
  title('\tau total (all species together)')
  xlabel('r /a')
elseif Nspec==1
  NTVtot=pbar*runs.NTV(:,1);
  NTVfromFluxtot=pbar*runs.NTVfromFlux(:,1);
  plot(runs.rN,-NTVtot,'g',...
       runs.rN,-NTVfromFluxtot,'g--')
  title('\tau total')
  xlabel('r /a')
end


integr=NTVtot'.*abs(4*pi^2./runs.FSABHat2.*(runs.GHat+runs.iota.*runs.IHat)*...
       psiAHat*Rbar^3);
%abs is taken because previously psiAHat had the wrong sign!

if 0
  fig(3)
  plot(s,integr)
  title('integrand for total NTV')
  xlabel('s')
end

%NTVtot0=NTVtot(1)-s(1)*(NTVtot(2)-NTVtot(1))/(s(2)-s(1));
%NTVtot1=NTVtot(end)+(1-s(end))*(NTVtot(end)-NTVtot(end-1))/(s(end)-s(end-1));



if 0 %EXTRAPOLATE
  integr0=integr(1)-s(1)*(integr(2)-integr(1))/(s(2)-s(1));
  integr1=integr(end)+(1-s(end))*(integr(end)-integr(end-1))/(s(end)-s(end-1));
  s=[0,s,1];
  integr=[integr0,integr,integr1];
end

NTV_Nm=trapz(s,integr)

A1=(runs.dnHatdpsiN(:,ion)./runs.nHats(:,ion)...
    +(runs.dPhiHatdpsiN'*Phibar).*runs.Zs(:,ion).*e./(runs.THats(:,ion)*Tbar)...
    -3/2*runs.dTHatdpsiN(:,ion)./runs.THats(:,ion))/psiAHat;
A2=runs.dTHatdpsiN(:,ion)./runs.THats(:,ion)/psiAHat;
kappaiFSAB2=vbar*Bbar*runs.FSABFlow(:,ion)./runs.nHats(:,ion)+...
    runs.THats(:,ion)*Tbar.*G./runs.Zs(:,ion)/e./iota.*...
    (A1+5/2*A2);
FSAomegai=1./(runs.GHat'+runs.iota'.*runs.IHat').*...
    (vbar*Bbar*runs.FSABFlow(:,ion)./runs.nHats(:,ion)...
     -runs.THats(:,ion)*Tbar.*runs.IHat'./runs.Zs(:,ion)/e.*...
     (A1+5/2*A2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Approximate analytical calculations in the ripple plateau
% for the AUG equilibrium rip_n16
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


FSAB2=zeros(length(runs.rN),1);
FSAg_phiphi=zeros(length(runs.rN),1);
integr=zeros(length(runs.rN),1);
FSAdelta2overB=zeros(length(runs.rN),1);

for ind=1:length(runs.rN)
  fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b%3i/%3i',ind,length(runs.rN))
  rnorm=runs.rN(ind);
  rind=findnearest(Geom.rnorm,rnorm);
  Ntheta=151;
  Nzeta=55;%Ntheta;
  [Ham,Booz]=makeHamada(Geom,rind,Ntheta,Nzeta);
  FSAg_phiphi(ind)=Booz.FSAg_phiphi;
  FSAB2(ind)=Booz.FSAB2;
  Btormean=mean(Booz.B,2)*ones(1,Nzeta);
  deltatw=(Booz.B-Btormean)./Btormean;
  %integr(ind)=4*pi^2/(Ntheta*Nzeta)*sum(sum(deltatw.^2./Btormean.^3))*2;
  FSAdelta2overB(ind)=2*sum(sum(deltatw.^2./Btormean.^3))/sum(sum(1./Booz.B.^2));
end
fprintf(1,'\n')

plat.tau=Geom.Nperiods*mi.^2.*vTi.^3.*(ni20*1e20)./(Z*e)./iota*...
          sqrt(pi)/4.*(G+iota.*I).*FSAdelta2overB.*(A1+A2*3);




%%%%%%%%%% This is for step 1 %%%%%%%%%%
fig(5)
if Nspec==2
  plot(runs.rN,runs.FSABFlow(:,1),runs.rN,runs.FSABFlow(:,2))
  legend('e','i')
else
  plot(runs.rN,runs.FSABFlow(:,1))
end
title('FSABFlow')
xlabel('r / a')

fig(6)
plot(runs.rN,kappaiFSAB2)
title('\kappa_i <B^2>')
xlabel('r / a')

tau_anal_over_tau_calc=(-plat.tau./-NTVtot)'
fig(7)
plot(runs.rN,-NTVtot,runs.rN,plat.tau)%,runs.rN,tauP)%*sqrt(2))%testing
%figure(1)
%hold on
%plot(runs.rN,plat.tau2)
%hold off
title('\tau')
legend('calc','anal')
xlabel('r / a')
