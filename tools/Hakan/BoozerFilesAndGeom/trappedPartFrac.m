function [ft,FSAB,FSAB2]=trappedPartFrac(Geom,rind,Ntheta,Nzeta)

%Boozer discretisation
Dzeta=2*pi/Nzeta/Geom.Nperiods;
Dtheta=2*pi/Ntheta;
zetavec=(0:Nzeta-1)*Dzeta;
thetavec=(0:Ntheta-1)'*Dtheta;
[theta, zeta] = ndgrid(thetavec,zetavec); %the same as meshgrid but 
                                          %more logical argument order

NPeriods=Geom.Nperiods;
Bmn=Geom.Bmn{rind}; %NOTE: Not normalised to B00
NHarmonics=Geom.nmodes(rind);

if Geom.StelSym
  parity=ones(size(Bmn));
else
  parity=Geom.parity{rind};
end

B = zeros(Ntheta,Nzeta);

for i=1:NHarmonics
  if parity(i) %The cosine components of B
    B = B + Bmn(i) *...
           cos(Geom.m{rind}(i) * theta - Geom.n{rind}(i) * NPeriods * zeta);
  else  %The sine components of B
    B = B + Bmn(i) *...
           sin(Geom.m{rind}(i) * theta - Geom.n{rind}(i) * NPeriods * zeta);
  end
end

Bmax=max(max(B));
sumsumBm2=sum(sum(B.^(-2)));
FSAB2=4*pi^2/ (NPeriods*Dzeta*Dtheta*sumsumBm2);
FSAB=sum(sum(B.^(-1)))/sumsumBm2;


%VHatprime=Geom.dVdsoverNper(rind)*NPeriods/Geom.torfluxtot/(G+iota*I)
%FSABHat2=4*pi^2/VHatprime

Nx=100;
x=linspace(0,1,Nx);
FSAg=zeros(1,Nx);
for xind=1:Nx
  FSAg(xind)=sum(sum(sqrt(1-x(xind)*B/Bmax)./B.^2))/sumsumBm2;
end

%fig(1)
%plot(x,FSAg)

ft=1-3/4*FSAB2/Bmax^2*trapz(x,x./FSAg);