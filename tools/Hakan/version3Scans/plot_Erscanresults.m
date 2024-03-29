function [runs,miss,ErVarroots]=plot_Erscanresults(directory,varargin)

isdirlistfile=0;
if length(directory)>4
  if directory(end-3:end)=='.dat'
    isdirlistfile=1;
  end
end

if isdirlistfile
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
  [runs,miss]=getresults(dirlist,'dPhiHatdpsiN');
else
  [runs,miss]=getresults(directory,'dPhiHatdpsiN');
end

if runs.NumElements==0
  error('Nothing found!')
end

if nargin>1
  ErVarStr=varargin{1};
else 
  ErVarStr='dPhiHatdpsiN';ErVar
end
ErVar=getfield(runs,ErVarStr);

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

plotNTV=0;

if plotNTV
  fig(11)
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
  
  
  fig(12)
  
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
  
  fig(13)
  plot(s,integr)
  title('integrand for total NTV')
  xlabel('s')
  
  NTV_Nm=trapz(s,integr);
  tau_Nm=-NTV_Nm
end


if 0
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
end


if 0
  %SATAKE FIGURE
  aw7x=0.51092;
  fig(7)
  dpsidr_w7x=psiAHat*2*runs.rN'/aw7x;
  Gamma=(runs.GHat'+runs.iota'.*runs.IHat')*4*pi./runs.FSABHat2'...
        *psiAHat*nbar*vbar.*runs.particleFlux_vm_psiN(:,1);
  Gamma_m2s=1./dpsidr_w7x ...
            *psiAHat*nbar*vbar.*runs.particleFlux_vm_psiN(:,1);
  
  Er_w7x=-runs.dPhiHatdpsiN'/psiAHat.*dpsidr_w7x;
  plot(Er_w7x,Gamma_m2s)
  ylabel('\Gamma [1/m^2s] ')
  xlabel('E_r [keV/m]')
  title('Satake case')
end

fig(1)
if Nspec==2
  semilogy(ErVar,runs.particleFlux_vm_psiN(:,1),...
           ErVar,runs.particleFlux_vm_psiN(:,2))
  legend('e','i')
elseif Nspec==3
  semilogy(ErVar,runs.particleFlux_vm_psiN(:,1),...
           ErVar,runs.particleFlux_vm_psiN(:,2),...
           ErVar,runs.particleFlux_vm_psiN(:,3))
  legend('e','i','Z')
else
  plot(ErVar,runs.particleFlux_vm_psiN(:,1))
end
title('particleFlux')
if strcmp(ErVarStr,'dPhiHatdpsiN')
  xlabel('d\Phi/d\psi_N')
else
  xlabel(ErVarStr)
end

if Nspec>1
   jnorm=sum(runs.particleFlux_vm_psiN.*runs.Zs,2);
   if any(diff(sign(jnorm)))
     ErVarroots=interp1(jnorm,ErVar,0,'pchip')
   else
     ErVarroots=NaN;
     ErVarroots_extrap=interp1(jnorm,ErVar,0,'pchip','extrap')
   end
   
   fig(7)
   plot(ErVar,jnorm,'b-',...
        ErVarroots,zeros(size(ErVarroots)),'r+')
   title('sum_{spec} (particleFlux*Z)')
   if strcmp(ErVarStr,'dPhiHatdpsiN')
     xlabel('d\Phi/d\psi_N')
   else
     xlabel(ErVarStr)
   end
   
end


fig(2)
if Nspec==2
  plot(ErVar,runs.heatFlux_vm_psiN(:,1),...
       ErVar,runs.heatFlux_vm_psiN(:,2))
  legend('e','i')
elseif Nspec==3
  plot(ErVar,runs.heatFlux_vm_psiN(:,1),...
       ErVar,runs.heatFlux_vm_psiN(:,2),...
       ErVar,runs.heatFlux_vm_psiN(:,3))
  legend('e','i','Z')
else
  plot(ErVar,runs.heatFlux_vm_psiN(:,1))
end
title('heatFlux')
if strcmp(ErVarStr,'dPhiHatdpsiN')
  xlabel('d\Phi/d\psi_N')
else
  xlabel(ErVarStr)
end


fig(3)
if Nspec==2
  plot(ErVar,runs.FSABFlow(:,1),...
       ErVar,runs.FSABFlow(:,2))
  legend('e','i')
elseif Nspec==3
  plot(ErVar,runs.FSABFlow(:,1),...
       ErVar,runs.FSABFlow(:,2),...
       ErVar,runs.FSABFlow(:,3))
  legend('e','i','Z')
else
  plot(ErVar,runs.FSABFlow(:,1))
end
title('FSABFlow')
if strcmp(ErVarStr,'dPhiHatdpsiN')
  xlabel('d\Phi/d\psi_N')
else
  xlabel(ErVarStr)
end




