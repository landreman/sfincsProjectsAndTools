function CoordPlots(Discr)
%Discr should be Ham or Booz as created by 
%[Ham,Booz]=makeHamada(Geom,rind,Ntheta,Nzeta);

Nperiods=Discr.Nperiods;
R=Discr.cylR;
Z=Discr.cylZ;
B=Discr.B;
zeta=Discr.zeta;
theta=Discr.theta;
phi=Discr.phi;
vthet=Discr.vthet;
cylphi=Discr.cylphi;
if Discr.cylPossible
  cylth=Discr.cylth;
  cylr=Discr.cylr;
end
cylR=Discr.cylR;
cylZ=Discr.cylZ;


zetaf=makefulltordata(zeta,Nperiods,'toroidal coordinate');
phif=makefulltordata(phi,Nperiods,'toroidal coordinate');
cylphif=makefulltordata(cylphi,Nperiods,'toroidal coordinate');
thetaf=makefulltordata(theta,Nperiods,'poloidal coordinate');
vthetf=makefulltordata(vthet,Nperiods,'poloidal coordinate');
if Discr.cylPossible
  cylthf=makefulltordata(cylth,Nperiods,'poloidal coordinate');
end

%make versions of coords which are within (0,2*pi)
zetafmod=mod(zetaf,2*pi);
thetafmod=mod(thetaf,2*pi);
phifmod=mod(phif,2*pi);
vthetfmod=mod(vthetf,2*pi);
cylphifmod=mod(cylphif,2*pi);
if Discr.cylPossible
  cylthfmod=mod(cylthf,2*pi);
end

geomang=-cylphi; %cylphi is minus the geometrical angle
X=R.*cos(geomang);
Y=R.*sin(geomang);

geomangf=-cylphif;
Rf=makefulltordata(R,Nperiods);
Xf=Rf.*cos(geomangf);
Yf=Rf.*sin(geomangf);
Zf=makefulltordata(Z,Nperiods);
Bf=makefulltordata(B,Nperiods);
if Discr.cylPossible
  cylrf=makefulltordata(cylr,Nperiods);
end
cylRf=makefulltordata(cylR,Nperiods);
cylZf=makefulltordata(cylZ,Nperiods);


fig(3)
surf(Xf,Yf,Zf,Bf)
%surf(X,Y,Z,B)
title('Magnetic field')
shading flat
axis equal

fig(1)
surf(Xf,Yf,Zf,phifmod)
%surf(X,Y,Z,phi)
title('Hamada toroidal coordinate')
shading flat
axis equal

fig(2)
surf(Xf,Yf,Zf,vthetfmod)
%surf(X,Y,Z,vthet)
title('Hamada poloidal coordinate')
shading flat
axis equal

fig(5)
surf(Xf,Yf,Zf,cylphifmod)
%surf(X,Y,Z,phi)
title('Cylindrical toroidal coordinate')
shading flat
axis equal

if Discr.cylPossible
  fig(6)
  surf(Xf,Yf,Zf,mod(cylthfmod,pi/2))
  %surf(X,Y,Z,vthet)
  title('Cylindrical poloidal coordinate')
  shading flat
  axis equal
  view(0,90);colorbar
end

fig(7)
surf(Xf,Yf,Zf,cylRf)
%surf(X,Y,Z,vthet)
title('Cylindrical major radius coordinate')
shading flat
axis equal
view(0,90);colorbar

if Discr.cylPossible
  fig(8)
  surf(Xf,Yf,Zf,cylrf)
  %surf(X,Y,Z,vthet)
  title('Cylindrical minor radius coordinate')
  shading flat
  axis equal
  colorbar
end

fig(12)
surf(Xf,Yf,Zf,cylZf)
%surf(X,Y,Z,vthet)
title('vertical distance')
shading flat
axis equal
colorbar

fig(9)
surf(Xf,Yf,Zf,zetafmod)
%surf(X,Y,Z,phi)
title('Boozer toroidal coordinate')
shading flat
axis equal

fig(10)
surf(Xf,Yf,Zf,thetafmod)
%surf(X,Y,Z,vthet)
title('Boozer poloidal coordinate')
shading flat
axis equal
view(0,90);colorbar