function L11=calcL11fromsimpleruns(runs)
%Use [runs,miss]=getconvergencescanresults(dir,1) first

if any(runs.dTHatdpsiN)
  error('runs.dTHatdpsiN elements not equal to zero exist! ')
end
if any(runs.EParallelHat)
  error('runs.EParallelHat elements not equal to zero exist! ')
end

Zs=runs.Zs';
mHats=runs.mHats';
nHats=runs.nHats';
THats=runs.THats';
dnHatdpsiN=runs.dnHatdpsiN';
dTHatdpsiN=runs.dTHatdpsiN';
particleFlux_vm_psiN=runs.particleFlux_vm_psiN';

A1tmp=dnHatdpsiN./nHats ...
      +(1e3*runs.dPhiHatdpsiN).*Zs/(THats*1e3) ...
      -3/2*dTHatdpsiN./THats;

e=1.6022e-19;
mp=1.6726e-27;
vT=sqrt(2*(1e3*THats*e)./(mHats*mp));
vbar=sqrt(2*1e3*e/mp);
vTbar=sqrt(THats./mHats);
nbar=1e20;

dpsiNdpsi=1./runs.psiAHat;

A1=A1tmp.*runs.GHat.*(THats*1e3)./Zs./runs.B0OverBBar./vT.*dpsiNdpsi;


normflux=nbar*vbar*particleFlux_vm_psiN./dpsiNdpsi.*Zs ...
         .*(runs.GHat+runs.iota.*runs.IHat) ...
         ./(nHats*nbar)./(THats*1e3)./runs.GHat;

%L11=normfluxvbar./A1vbar;
L11=normflux./A1;



%%%%%%%% nuPrime calculation
loglambdadef='17.3'
if loglambdadef=='17.3'
  loglambda=17.3
elseif loglambdadef=='HMCB'
  ind1=find(THats<0.050);
  ind2=find(THats>=0.050);
  if not(isempty(ind1))
    loglambda(ind1)=23.4-1.15*log10(nHats(ind1)*1e14)+3.45*log10(THats(ind1)*1e3);
  end
  if not(isempty(ind2))
    loglambda(ind2)=25.3-1.15*log10(nHats(ind2)*1e14)+2.30*log10(THats(ind2)*1e3);
  end
end

epsilon_0=8.8542e-12;
nu_ii = 4*sqrt(2*pi)*(nHats*1e20).*Zs.^2*e^4.*loglambda./ ...
        (3*(4*pi*epsilon_0)^2.*sqrt(mp*mHats).*(THats*1e3*e).^(3/2))

nuPrime=nu_ii./sqrt(2*e*1e3/mp./mHats)...
        .*(runs.GHat+runs.iota.*runs.IHat)./runs.B0OverBBar./sqrt(THats)