function CoordPlots(Discr,varargin)
% Discr should be the discretisation Ham, Booz or Pest as created by 
% [Ham,Booz,Cyl,Pest]=makeHamada(Geom,rind,Ntheta,Nzeta);
%
% whichtoplot=varargin{1} is a string telling which coordinates to plot
% If not provided, all will be plotted
% Otherwise the options are
% whichtoplot='Hamada'
% whichtoplot='Boozer'
% whichtoplot='Cylindrical'
% whichtoplot='Pest'
% whichtoplot='Aligned'
% whichtoplot='all'

if nargin==1
  whichtoplot='all';
else
  whichtoplot=varargin{1};
end


Nperiods=Discr.Nperiods;
R=Discr.cylR;
Z=Discr.cylZ;
B=Discr.B;
zeta=Discr.zeta;
theta=Discr.theta;
phi=Discr.phi;
vthet=Discr.vthet;
pzeta=Discr.pzeta;
ptheta=Discr.ptheta;
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
pzetaf=makefulltordata(pzeta,Nperiods,'toroidal coordinate');
thetaf=makefulltordata(theta,Nperiods,'poloidal coordinate');
vthetf=makefulltordata(vthet,Nperiods,'poloidal coordinate');
pthetaf=makefulltordata(ptheta,Nperiods,'poloidal coordinate');
if Discr.cylPossible
  cylthf=makefulltordata(cylth,Nperiods,'poloidal coordinate');
end

%make versions of coords which are within (0,2*pi)
zetafmod=mod(zetaf,2*pi);
thetafmod=mod(thetaf,2*pi);
phifmod=mod(phif,2*pi);
vthetfmod=mod(vthetf,2*pi);
pzetafmod=mod(pzetaf,2*pi);
pthetafmod=mod(pthetaf,2*pi);
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



fig(8)
surf(Xf,Yf,Zf,Bf)
%surf(X,Y,Z,B)
title('Magnetic field')
shading flat
axis equal


if strcmp(whichtoplot,'Hamada') || strcmp(whichtoplot,'all')
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
end

if strcmp(whichtoplot,'Cylindrical') || strcmp(whichtoplot,'all')
  fig(9)
  surf(Xf,Yf,Zf,cylphifmod)
  %surf(X,Y,Z,phi)
  title('Cylindrical toroidal coordinate')
  shading flat
  axis equal

  if Discr.cylPossible
    fig(10)
    surf(Xf,Yf,Zf,mod(cylthfmod,pi/2))
    %surf(X,Y,Z,vthet)
    title('Cylindrical poloidal coordinate')
    shading flat
    axis equal
    view(0,90);colorbar
  end

  fig(11)
  surf(Xf,Yf,Zf,cylRf)
  %surf(X,Y,Z,vthet)
  title('Cylindrical major radius coordinate')
  shading flat
  axis equal
  view(0,90);colorbar

  if Discr.cylPossible
    fig(12)
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
end

if strcmp(whichtoplot,'Boozer') || strcmp(whichtoplot,'all')
  fig(1+2*strcmp(whichtoplot,'all'))
  surf(Xf,Yf,Zf,zetafmod)
  %surf(X,Y,Z,phi)
  title('Boozer toroidal coordinate')
  shading flat
  axis equal

  fig(2+2*strcmp(whichtoplot,'all'))
  surf(Xf,Yf,Zf,thetafmod)
  %surf(X,Y,Z,vthet)
  title('Boozer poloidal coordinate')
  shading flat
  axis equal
  %view(0,90);colorbar
end

if strcmp(whichtoplot,'Pest') || strcmp(whichtoplot,'all')
  fig(1+4*strcmp(whichtoplot,'all'))
  surf(Xf,Yf,Zf,pzetafmod)
  %surf(X,Y,Z,phi)
  title('Pest toroidal coordinate')
  shading flat
  axis equal

  fig(2+4*strcmp(whichtoplot,'all'))
  surf(Xf,Yf,Zf,pthetafmod)
  %surf(X,Y,Z,vthet)
  title('Pest poloidal coordinate')
  shading flat
  axis equal
  %view(0,90);colorbar
end


if Discr.Nperiods==1 %Aligned coordinates only work fo this case
  if strcmp(whichtoplot,'Align') || strcmp(whichtoplot,'all')
    if strcmp(Discr.name,'Cylindrical')
      error('Cannot make field aligned coordinates from non-straight-field-line coordinates!')
    elseif strcmp(Discr.name,'Boozer') || strcmp(Discr.name,'Hamada') || strcmp(Discr.name,'Pest')
      [Align,Discr]=makeAlign(Discr);
      if strcmp(Discr.name,'Boozer')
        atorf=makefulltordata(Discr.alzeta,Nperiods,'aligned toroidal coordinate',Discr.iota);
        apolf=makefulltordata(Discr.altheta,Nperiods,'poloidal coordinate');
      elseif strcmp(Discr.name,'Hamada')
        atorf=makefulltordata(Discr.alphi,Nperiods,'aligned toroidal coordinate',Discr.iota);
        apolf=makefulltordata(Discr.alvthet,Nperiods,'poloidal coordinate');
      elseif strcmp(Discr.name,'Pest')
        atorf=makefulltordata(Discr.alpzeta,Nperiods,'aligned toroidal coordinate',Discr.iota);
        apolf=makefulltordata(Discr.alptheta,Nperiods,'poloidal coordinate');
      end     
      fig(1+6*strcmp(whichtoplot,'all'))
      surf(Xf,Yf,Zf,atorf)
      %surf(X,Y,Z,phi)
      title(['Aligned ',Discr.name,' toroidal coordinate'])
      shading flat
      axis equal      
      
      fig(2+6*strcmp(whichtoplot,'all'))
      surf(Xf,Yf,Zf,apolf)
      %surf(X,Y,Z,vthet)
      title(['Aligned ',Discr.name,' poloidal coordinate'])
      shading flat
      axis equal
      
    else
      error('Non-recognised Discr.name!')
    end
  end
end

%atorf/2/pi