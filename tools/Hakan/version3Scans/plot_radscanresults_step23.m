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
    if tline(1)~='%' && tline(1)~='!'
      dirlist={dirlist{:},tline};
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
pbar=nbar*Tbar;
mi=runs2.mHats(:,2)*mp;
ni=runs2.nHats(:,2)*nbar;
alpha=1;
Deltahalf=runs2.Delta' / 2;
vbar=sqrt(1e3*e*2/mp);

tau1=runs2.NTV*pbar;
tauin=runs3.NTV*pbar;

if runs2.NumElements~=runs3.NumElements
  error('The number of finished runs in step 2 and 3 are different!')
end
if any(runs2.rN~=runs3.rN)
  error('The radii in step 2 and three do not match')
end

%We can extract omega_meas from the input data of step 2.
%the density gradient was calucated using
%dni20dPsiN_step2= -Z  *ni20./TikeV.*...
%      (alpha*dPotentialkVdPsiN + iota.*halfDelta.*omega_torrot/vbar);

omega_torrot=vbar.*Deltahalf./runs2.iota' .* ( ...
    -runs2.dnHatdpsiN(:,2)./runs2.Zs(:,2)./runs2.nHats(:,2).*runs2.THats(:,2)...
    -alpha*runs2.dPhiHatdpsiN');

dPhiHatdpsiN=runs2.dPhiHatdpsiN
runs2.dnHatdpsiN(:,2)
runs2.iota

if 0 %approximately
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

nut = -tau1(:,2)./mi./ni.*FSAg_phiphi./omega_torrot;
omegain = -tau1(:,2)./tauin(:,2).*omega_torrot;

NTVtot=(tau1(:,1)+tauin(:,1)+tau1(:,2)+tauin(:,2))';
integr=NTVtot.*4*pi./runs2.FSABHat2.*(runs2.GHat+runs2.iota.*runs2.IHat);
s=sqrt(runs2.rN);
NTV_Nm=trapz(s,integr)


fig(1)
plot(runs2.rN,omega_torrot,runs2.rN,omegain)
xlabel('r / a')
ylabel('\omega')
legend('\omega_{meas}','\omega_{in}')


fig(2)
plot(runs2.rN,1./nut)
xlabel('r / a')
ylabel('1 / \nu_t')

fig(3)
plot(runs2.rN,tau1(:,1)+tauin(:,1),runs2.rN,tau1(:,2)+tauin(:,2))
title('torque')
xlabel('r / a')
ylabel('Nm')

