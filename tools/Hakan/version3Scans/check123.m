function check123(step1dir,step2dir,step3dir)

if step1dir(end-3:end)=='.dat'
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
  [runs1,miss1]=getresults(dirlist);
else
  [runs1,miss1]=getresults(step1dir);
end

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
ion=1;
mi=runs2.mHats(:,ion)*mp;
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

omega_torrot=vbar.*Deltahalf./iota/psiAHat .* ( ...
    -runs2.dnHatdpsiN(:,ion)./runs2.Zs(:,ion)...
    ./runs2.nHats(:,ion).*runs2.THats(:,ion)...
    -alpha*runs2.dPhiHatdpsiN');


dTikeVdPsiN=runs3.dTHatdpsiN(:,ion);
dPotentialkVdPsiNLHS=runs2.dPhiHatdpsiN';
if any(runs3.dPhiHatdpsiN'~=runs2.dPhiHatdpsiN')
  error('not the same LHS Er in steps 2 and 3!')
end
dpdpsi=(ni20.*dTikeVdPsiN+TikeV.*runs1.dnHatdpsiN(1:8,ion))...
       *1e20*1e3*e/psiAHat;
FSAVBistep1=vbar*Bbar*runs1.FSABFlow(:,ion)./runs1.nHats(:,ion);

ELHSshould=-iota.*(omega_torrot+1./(G+iota.*I).*(I./ni./Z/e.*dpdpsi-FSAVBistep1(1:8)));

%fig(20)
%plot(runs1.rN,runs1.FSABFlow(:,ion))


omega_test=-dPotentialkVdPsiNLHS./iota/psiAHat...
    -1./(G+iota.*I).*(I./ni./Z/e.*dpdpsi-FSAVBistep1(1:8));

% Calculate thermodyn. forces
A1istep1=(runs1.dnHatdpsiN(:,ion)./runs1.nHats(:,ion) ...
          -3/2*runs1.dTHatdpsiN(:,ion)./runs1.THats(:,ion))/psiAHat;
A2istep1=(runs1.dTHatdpsiN(:,ion)./runs1.THats(:,ion))/psiAHat;

A1istep2=(runs2.dnHatdpsiN(:,ion)./ni20 ...
          +Z./TikeV.*alpha.*dPotentialkVdPsiNLHS)/psiAHat;
A2istep2=zeros(size(A1istep2));

A1istep3=(runs3.dnHatdpsiN(:,ion)./ni20 ...
          +Z./TikeV.*alpha.*dPotentialkVdPsiNLHS ...
          -3/2*dTikeVdPsiN./TikeV)/psiAHat;
A2istep3=(dTikeVdPsiN./TikeV)/psiAHat;


kappaiFSAB2_step1=vbar*Bbar*runs1.FSABFlow(:,ion)./runs1.nHats(:,ion)+...
      runs1.THats(:,ion)*Tbar.*runs1.GHat'./runs1.Zs(:,ion)/e./runs1.iota'.*...
      (A1istep1+5/2*A2istep1);
kappaiFSAB2_step2=vbar*Bbar*runs2.FSABFlow(:,ion)./runs2.nHats(:,ion)+...
      runs2.THats(:,ion)*Tbar.*G./runs2.Zs(:,ion)/e./iota.*...
      (A1istep2+5/2*A2istep2);
kappaiFSAB2_step3=vbar*Bbar*runs3.FSABFlow(:,ion)./runs3.nHats(:,ion)+...
      runs3.THats(:,ion)*Tbar.*G./runs3.Zs(:,ion)/e./iota.*...
      (A1istep3+5/2*A2istep3);


FSAomegai_step1=1./(runs1.GHat'+runs1.iota'.*runs1.IHat').*...
    (vbar*Bbar*runs1.FSABFlow(:,ion)./runs1.nHats(:,ion)...
     -runs1.THats(:,ion)*Tbar.*runs1.IHat'./runs1.Zs(:,ion)/e.*...
     (A1istep1+5/2*A2istep1));
FSAomegai_step2=1./(G+iota.*I).*...
    (vbar*Bbar*runs2.FSABFlow(:,ion)./runs2.nHats(:,ion)...
     -runs2.THats(:,ion)*Tbar.*I./runs2.Zs(:,ion)/e.*...
     (A1istep2+5/2*A2istep2));
FSAomegai_step3=1./(G+iota.*I).*...
    (vbar*Bbar*runs3.FSABFlow(:,ion)./runs3.nHats(:,ion)...
     -runs3.THats(:,ion)*Tbar.*I./runs3.Zs(:,ion)/e.*...
     (A1istep3+5/2*A2istep3));


fig(7)
plot(runs1.rN,kappaiFSAB2_step1,...
     runs2.rN,kappaiFSAB2_step2,...
     runs3.rN,kappaiFSAB2_step3,...
     runs2.rN,kappaiFSAB2_step2+kappaiFSAB2_step3)
title('\kappa_i <B^2>')
xlabel('r / a')
legend('step 1','step 2','step 3','2 + 3')

fig(8)
plot(runs1.rN,FSAomegai_step1,'b',...
     runs2.rN,FSAomegai_step2,'m',...
     runs3.rN,FSAomegai_step3,'c',...
     runs2.rN,FSAomegai_step2+FSAomegai_step3,'r--',...
     runs2.rN,omega_torrot,'k',...
     runs2.rN,omega_test,'y')
title('<\omega_i>')
xlabel('r / a')
legend('step 1','step 2','step 3','2 + 3','input')

fig(9)
plot(runs2.rN,ELHSshould*psiAHat/1e3,...
     runs2.rN,dPotentialkVdPsiNLHS)