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
pbar=nbar*Tbar;
psiAHat=runs.psiAHat(1);
vbar=sqrt(1e3*e*2/mp);
iota=runs.iota';
G=runs.GHat'*Bbar*Rbar;
psiAHat=runs.psiAHat(1);

Nspec=size(runs.NTV,2);

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

fig(3)
plot(s,integr)
title('integrand for total NTV')
xlabel('s')


%NTVtot0=NTVtot(1)-s(1)*(NTVtot(2)-NTVtot(1))/(s(2)-s(1));
%NTVtot1=NTVtot(end)+(1-s(end))*(NTVtot(end)-NTVtot(end-1))/(s(end)-s(end-1));



if 0 %EXTRAPOLATE
  integr0=integr(1)-s(1)*(integr(2)-integr(1))/(s(2)-s(1));
  integr1=integr(end)+(1-s(end))*(integr(end)-integr(end-1))/(s(end)-s(end-1));
  s=[0,s,1];
  integr=[integr0,integr,integr1];
end

NTV_Nm=trapz(s,integr)

ion=find(runs.Zs(1,:)~=-1);
A1=(runs.dnHatdpsiN(:,ion)./runs.nHats(:,ion)...
    +runs.dPhiHatdpsiN'.*runs.Zs(:,ion).*e./(runs.THats(:,ion)*Tbar)...
    -3/2*runs.dTHatdpsiN(:,ion)./runs.THats(:,ion))/psiAHat;
A2=runs.dTHatdpsiN(:,ion)./runs.THats(:,ion)/psiAHat;
kappaiFSAB2=vbar*Bbar*runs.FSABFlow(:,ion)./runs.nHats(:,ion)+...
    runs.THats(:,ion)*Tbar.*G./runs.Zs(:,ion)/e./iota.*...
    (A1+5/2*A2);
FSAomegai=1./(runs.GHat'+runs.iota'.*runs.IHat').*...
    (vbar*Bbar*runs.FSABFlow(:,ion)./runs.nHats(:,ion)...
     -runs.THats(:,ion)*Tbar.*runs.IHat'./runs.Zs(:,ion)/e.*...
     (A1+5/2*A2));


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

