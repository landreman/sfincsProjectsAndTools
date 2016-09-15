function [Ham,Booz,Cyl,Pest]=makeHamada(Geom,rind,Ntheta,Nzeta,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function transforms Boozer coordinate data to other coordinate systems,
% Hamada, Pest and (toroidally) Cylindrical.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The name of the output struct signifies the uniform grid coordinates upon which
% included quantities are discretised. The struct also contains values of the 
% other coordinates on that grid.
% 
% The input struct Geom which contains the Boozer coordinate data 
% can be generated with readBoozerfile.m, readBoozerXformfile.m
% or makeBcfromVmec.m
%
% rind is the index of the chosen flux surface
%
% Ntheta and Nzeta is the resolution of the spatial grid used in the
% transformation to other coordinates. Ntheta and Nzeta must both be odd!
% Typical values are 101 - 151. Larger values make the calculation heavy.
%
% Ntheta and Nzeta should be chosen to resolve the necessary toroidal
% and poloidal mode numbers of the magnetic field. If you set Ntheta to Nzeta to
% less than the numbers reuired to resolve the modes of Bmn in Geom, an error message
% is produced. If you still want these low Ntheta and Nzeta, you can set the optional
% argument varargin to be the string 'forceSize'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The name of the output struct denotes the uniform grid coordinates upon which
% included quantities are discretised. The struct also contains values of the 
% other coordinates on that grid.
%
% The basic output is the same in all four discretisation output structs. 
% Let us call one of the output structs Discr in general, i.e. 
% Discr = Ham, Booz, Cyl or Pest. The Fourier representation of B, R and Z in the respective 
% coordinates is stored as
% Discr.m, Discr.n, Discr.parity, Discr.Bmn, Discr.Rmn and Discr.Zmn
% where parity=1 signifies a cosinus term for Bmn and Rmn, but a sinus term for Zmn.
%
% Names of coordinates (right-handed):
%
% Boozer:      psi,  theta,  zeta
% Hamada:      psi,  vthet,  phi   (Note: psi is the radial coordinate, not V)
% Cylindrical: cylr, cylth,  cylphi (Not always possible to construct)
% Pest:        psi,  ptheta, pzeta
%
% Jacobians are given in the respective discretisation structs as
% Ham.Jacob_psi_vthet_phi
% Booz.Jacob_psi_theta_zeta
%
% Some differences between the coordinates are also given. They are denoted with a D
% followed by the names of the two coordinates. For instance, Booz.Dzetaphi is the difference
% between zeta and phi at the discretisation points of the uniform (theta,zeta) grid. 
%
% Some metric elements are also given. The notation is that, e.g.,
% Booz.g_thetazeta is g_{\theta\zeta} = \mathbf{e}_\theta \cdot \mathbf{e}_\zeta
% Booz.gpsipsi     is g^{\psi\psi} = \nabla\psi\cdot\nabla\psi
%
% The quantity u satisfies
% (B dot grad) u = 2 iota |B|^-3 B x nabla psi dot nabla |B|
%
% G and I are defined by the expression for the magnetic field
% \mathbf{B} = I\nabla\theta + G\nabla\zeta + B_\psi(\theta,\zeta)\nabla\psi
% \mathbf{B} = I\nabla\vthet + G\nabla\phi + \nabla H(\psi,\vthet,\phi)
%
% where Booz.Bpsitilde is B_\psi up to an additive constant. 
% (The m=0,n=0 comp. of Booz.Bpsitilde is set to zero.)
%
% Other quantities:
% h = 1/B^2
% FSAB2       = <B^2>   (FSA stands for flux surface average)
% FSAu2B2     = <u^2 B^2>
% FSAgpsipsi  = <g^{\psi\psi}>
% FSAg_phiphi = <g_{\phi\phi}>
% (CylR,CylZ,cylphi) or (CylR,-cylphi,CylZ): Coordinates of discretisation points in
%                                            the "major radius cylinder"
% (X,Y,CylZ) : Cartesian coordinates of discretisation points
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

useFFT = 1; %to compare results of 0 and 1 here, one cannot use the ifftOpt='forceSize' option

if isempty(rind)
  error('rind is empty!')
end
ifftOpt='';
if nargin==5
  ifftOpt=varargin{1};
end

%identifying names for the coordinate systems:
Booz.name='Boozer';
Ham.name='Hamada';
Cyl.name='Cylindrical';
Pest.name='Pest';

NPeriods=Geom.Nperiods;
I=Geom.Btheta(rind);
G=Geom.Bphi(rind);
iota=Geom.iota(rind);
Bmn=Geom.Bmn{rind}; %NOTE: Not normalised to B00
NHarmonics=Geom.nmodes(rind);
mu0=1.256637061e-06;
mu0dpdpsi=Geom.dpds(rind)/Geom.psi_a*mu0; %[T/m^2]

%Boozer discretisation
Dzeta=2*pi/Nzeta/Geom.Nperiods;
Dtheta=2*pi/Ntheta;
zetavec=(0:Nzeta-1)*Dzeta;
thetavec=(0:Ntheta-1)'*Dtheta;
[theta, zeta] = ndgrid(thetavec,zetavec); %the same as meshgrid but 
                                          %more logical argument order

use4Dshortvecs=0;
if use4Dshortvecs %This approach does not save any time
  %short vectors S
  mmaxS=max(Geom.m{rind});
  mminS=0;
  mvecS=mminS:mmaxS;
  nmaxS=max(abs(Geom.n{rind}));
  nminS=-nmaxS;
  nvecS=nminS:nmaxS;
  tic
    [theta4D, zeta4D, m4D, n4D] = ndgrid(thetavec,zetavec,mvecS,nvecS);
    c4D=cos(m4D .* theta4D - n4D * NPeriods .* zeta4D);
    s4D=sin(m4D .* theta4D - n4D * NPeriods .* zeta4D);
  toc
end

%This approach does in deed save time because nvecL is a long vector
%nvecL=-floor(Nzeta/2):(floor(Nzeta/2)-1);
%nvecL0=1:floor(Nzeta/2)-1;
%[theta3D, zeta3D, n3D] = ndgrid(thetavec,zetavec,nvecL);
%[theta3D0, zeta3D0, n3D0] = ndgrid(thetavec,zetavec,nvecL0);
Ntheta_even=not(mod(Ntheta,2));
Nzeta_even=not(mod(Nzeta,2));
nvecL=-floor(Nzeta/2)+Nzeta_even:floor(Nzeta/2);
nvecL0c=0:floor(Nzeta/2);
[theta3D, zeta3D, n3D] = ndgrid(thetavec,zetavec,nvecL);
[theta3D0c, zeta3D0c, n3D0c] = ndgrid(thetavec,zetavec,nvecL0c);


if Geom.StelSym
  parity=ones(size(Bmn));
else
  parity=Geom.parity{rind};
end

if useFFT 
  %tic
  if 0
    Blist.cosparity=parity;
    Blist.m=Geom.m{rind};
    Blist.n=Geom.n{rind};
    Blist.data=Bmn;
    Rlist.cosparity=parity;
    Rlist.m=Geom.m{rind};
    Rlist.n=Geom.n{rind};
    Rlist.data=Geom.R{rind};
    Zlist.cosparity=not(parity);
    Zlist.m=Geom.m{rind};
    Zlist.n=Geom.n{rind};
    Zlist.data=Geom.Z{rind};
    Dzetacylphilist.cosparity=not(parity);
    Dzetacylphilist.m=Geom.m{rind};
    Dzetacylphilist.n=Geom.n{rind};
    Dzetacylphilist.data=2*pi/NPeriods*Geom.Dphi{rind};
    %[B,R,Z,Dzetacylphi]=ifftmn(mnmat([Blist,Rlist,Zlist,Dzetacylphilist],Ntheta,Nzeta,ifftOpt),NPeriods,Ntheta,Nzeta,ifftOpt);
    [B,R,Z,Dzetacylphi]=ifftmn(mnmat([Blist,Rlist,Zlist,Dzetacylphilist],Ntheta,Nzeta,ifftOpt),NPeriods);
    [dBdtheta,dBdzeta,dRdtheta,dRdzeta,dZdtheta,dZdzeta,...
     dDzetacylphi_dtheta,dDzetacylphi_dzeta]=...
        ifftmn(mnmat([mnlistgrad(Blist,NPeriods),mnlistgrad(Rlist,NPeriods),...
                      mnlistgrad(Zlist,NPeriods),mnlistgrad(Dzetacylphilist,NPeriods)],...
                     Ntheta,Nzeta,ifftOpt),NPeriods);
  else
    mnmats.B=mnmat(Geom,rind,'B',Ntheta,Nzeta,ifftOpt);
    mnmats.R=mnmat(Geom,rind,'R',Ntheta,Nzeta,ifftOpt);
    mnmats.Z=mnmat(Geom,rind,'Z',Ntheta,Nzeta,ifftOpt);
    mnmats.Dzetacylphi=mnmat(Geom,rind,'Dzetacylphi',Ntheta,Nzeta,ifftOpt);
    [B,R,Z,Dzetacylphi]=ifftmn([mnmats.B,mnmats.R,mnmats.Z,mnmats.Dzetacylphi]);
    [mnmats.dBdtheta,mnmats.dBdzeta]=grad(mnmats.B,NPeriods);
    [mnmats.dRdtheta,mnmats.dRdzeta]=grad(mnmats.R,NPeriods);
    [mnmats.dZdtheta,mnmats.dZdzeta]=grad(mnmats.Z,NPeriods);
    [mnmats.dDzetacylphidtheta,mnmats.dDzetacylphidzeta]=...
        grad(mnmats.Dzetacylphi,NPeriods);
    [dBdtheta,dBdzeta]=ifftmn([mnmats.dBdtheta,mnmats.dBdzeta]);
    [dRdtheta,dRdzeta]=ifftmn([mnmats.dRdtheta,mnmats.dRdzeta]);
    [dZdtheta,dZdzeta]=ifftmn([mnmats.dZdtheta,mnmats.dZdzeta]);
    [dDzetacylphi_dtheta,dDzetacylphi_dzeta]=...
        ifftmn([mnmats.dDzetacylphidtheta,mnmats.dDzetacylphidzeta]);

    [d2Rdtheta2,d2Rdthetadzeta]=ifftmn(grad(mnmats.dRdtheta,NPeriods));
    [d2Rdthetadzeta,d2Rdzeta2] =ifftmn(grad(mnmats.dRdzeta,NPeriods));
    [d2Zdtheta2,d2Zdthetadzeta]=ifftmn(grad(mnmats.dZdtheta,NPeriods));
    [d2Zdthetadzeta,d2Zdzeta2] =ifftmn(grad(mnmats.dZdzeta,NPeriods));
    [d2Dzetacylphi_dtheta2,d2Dzetacylphi_dthetadzeta]=ifftmn(grad(mnmats.dDzetacylphidtheta,NPeriods));
    [d2Dzetacylphi_dthetadzeta,d2Dzetacylphi_dzeta2] = ...
        ifftmn(grad(mnmats.dDzetacylphidzeta,NPeriods));
  end
  %toc
else %not(useFFT)
  %tic
  B           = zeros(Ntheta,Nzeta);
  dBdtheta    = zeros(Ntheta,Nzeta);
  dBdzeta     = zeros(Ntheta,Nzeta);
  R           = zeros(Ntheta,Nzeta);
  dRdtheta    = zeros(Ntheta,Nzeta);
  dRdzeta     = zeros(Ntheta,Nzeta);
  Z           = zeros(Ntheta,Nzeta);
  dZdtheta    = zeros(Ntheta,Nzeta);
  dZdzeta     = zeros(Ntheta,Nzeta);
  Dzetacylphi = zeros(Ntheta,Nzeta);
  dDzetacylphi_dtheta = zeros(Ntheta,Nzeta);
  dDzetacylphi_dzeta  = zeros(Ntheta,Nzeta);

  d2Rdtheta2     = zeros(Ntheta,Nzeta);
  d2Rdzeta2      = zeros(Ntheta,Nzeta);
  d2Rdthetadzeta = zeros(Ntheta,Nzeta);
  d2Zdtheta2     = zeros(Ntheta,Nzeta);
  d2Zdzeta2      = zeros(Ntheta,Nzeta);
  d2Zdthetadzeta = zeros(Ntheta,Nzeta);
  d2Dzetacylphi_dtheta2     = zeros(Ntheta,Nzeta);
  d2Dzetacylphi_dzeta2      = zeros(Ntheta,Nzeta);
  d2Dzetacylphi_dthetadzeta = zeros(Ntheta,Nzeta);
  
  for i=1:NHarmonics
    m=Geom.m{rind}(i);
    n=Geom.n{rind}(i);
    
    if use4Dshortvecs %preparation of c4D,s4D takes too much time
      mind=find(mvecS==m);
      nind=find(nvecS==n);  
      c=squeeze(c4D(:,:,mind,nind));
      s=squeeze(s4D(:,:,mind,nind));
    else
      c=cos(m * theta - n * NPeriods * zeta);
      s=sin(m * theta - n * NPeriods * zeta);
    end
    if parity(i) %The cosine components of B
      B = B + Bmn(i) * c;
      dBdtheta = dBdtheta - Bmn(i) * m * s;
      dBdzeta  = dBdzeta  + Bmn(i) * n * NPeriods * s;
      
      R    = R + Geom.R{rind}(i) * c;
      dRdtheta = dRdtheta - Geom.R{rind}(i) * m * s;
      dRdzeta  = dRdzeta  + Geom.R{rind}(i) * n * NPeriods * s;
      d2Rdtheta2 = d2Rdtheta2 - Geom.R{rind}(i) * m^2 * c;
      d2Rdzeta2  = d2Rdzeta2 - Geom.R{rind}(i) * (n * NPeriods)^2 * c;
      d2Rdthetadzeta=d2Rdthetadzeta + Geom.R{rind}(i) * m * n * NPeriods * c;
      
      Z    = Z + Geom.Z{rind}(i) * s;
      dZdtheta = dZdtheta + Geom.Z{rind}(i) * m * c;
      dZdzeta  = dZdzeta  - Geom.Z{rind}(i) * n * NPeriods * c;
      d2Zdtheta2 = d2Zdtheta2 - Geom.Z{rind}(i) * m^2 * s;
      d2Zdzeta2  = d2Zdzeta2 - Geom.Z{rind}(i) * (n * NPeriods)^2 * s;
      d2Zdthetadzeta=d2Zdthetadzeta + Geom.Z{rind}(i) * m * n * NPeriods * s;
      
      Dzetacylphi = Dzetacylphi + 2*pi/NPeriods*Geom.Dphi{rind}(i) * s;
      dDzetacylphi_dtheta = dDzetacylphi_dtheta + 2*pi/NPeriods*Geom.Dphi{rind}(i) * m * c;
      dDzetacylphi_dzeta  = dDzetacylphi_dzeta - 2*pi/NPeriods*Geom.Dphi{rind}(i) * n * NPeriods * c;
      d2Dzetacylphi_dtheta2 = d2Dzetacylphi_dtheta2 - 2*pi/NPeriods*Geom.Dphi{rind}(i) * m^2 * s;
      d2Dzetacylphi_dzeta2  = d2Dzetacylphi_dzeta2 - 2*pi/NPeriods*Geom.Dphi{rind}(i) * (n * NPeriods)^2 * s;
      d2Dzetacylphi_dthetadzeta=d2Dzetacylphi_dthetadzeta + 2*pi/NPeriods*Geom.Dphi{rind}(i) * m * n * NPeriods * s;
      
    else  %The sine components of B
      B = B + Bmn(i) * s;
      dBdtheta = dBdtheta + Bmn(i) * m * c;
      dBdzeta  = dBdzeta  - Bmn(i) * n * NPeriods * c; 
      
      R    = R + Geom.R{rind}(i) * s;
      dRdtheta = dRdtheta + Geom.R{rind}(i) * m * c;
      dRdzeta  = dRdzeta  - Geom.R{rind}(i) * n * NPeriods * c; 
      d2Rdtheta2 = d2Rdtheta2 - Geom.R{rind}(i) * m^2 * s;
      d2Rdzeta2  = d2Rdzeta2 - Geom.R{rind}(i) * (n * NPeriods)^2 * s;
      d2Rdthetadzeta=d2Rdthetadzeta + Geom.R{rind}(i) * m * n * NPeriods * s;      
      
      Z    = Z + Geom.Z{rind}(i) * c;
      dZdtheta = dZdtheta - Geom.Z{rind}(i) * m * s;
      dZdzeta  = dZdzeta  + Geom.Z{rind}(i) * n * NPeriods * s;
      d2Zdtheta2 = d2Zdtheta2 - Geom.Z{rind}(i) * m^2 * c;
      d2Zdzeta2  = d2Zdzeta2 - Geom.Z{rind}(i) * (n * NPeriods)^2 * c;
      d2Zdthetadzeta=d2Zdthetadzeta + Geom.Z{rind}(i) * m * n * NPeriods * c;
      
      Dzetacylphi = Dzetacylphi + 2*pi/NPeriods*Geom.Dphi{rind}(i) * c;
      dDzetacylphi_dtheta = dDzetacylphi_dtheta - 2*pi/NPeriods*Geom.Dphi{rind}(i) * m * s;
      dDzetacylphi_dzeta  = dDzetacylphi_dzeta + 2*pi/NPeriods*Geom.Dphi{rind}(i) * n * NPeriods * s;
      d2Dzetacylphi_dtheta2 = d2Dzetacylphi_dtheta2 - 2*pi/NPeriods*Geom.Dphi{rind}(i) * m^2 * c;
      d2Dzetacylphi_dzeta2  = d2Dzetacylphi_dzeta2 - 2*pi/NPeriods*Geom.Dphi{rind}(i) * (n * NPeriods)^2 * c;
      d2Dzetacylphi_dthetadzeta=d2Dzetacylphi_dthetadzeta + 2*pi/NPeriods*Geom.Dphi{rind}(i) * m * n * NPeriods * c;
      
    end
  end
  %toc
end

if any(size(zeta)~=size(Dzetacylphi)) && not(strcmp(ifftOpt,'forceSize'))
  error('Size mismatch! consider using the forceSize option.')
end

cylphi=zeta-Dzetacylphi;     %cylphi is minus the geometrical toroidal angle.
geomang=-cylphi;             %(R,Z,cylphi) and (R,geomang,Z) are right handed systems. 
dgeomangdtheta=dDzetacylphi_dtheta;
dgeomangdzeta =dDzetacylphi_dzeta - 1;
d2geomangdtheta2=d2Dzetacylphi_dtheta2;
d2geomangdzeta2=d2Dzetacylphi_dzeta2;
d2geomangdthetadzeta=d2Dzetacylphi_dthetadzeta;
X=R.*cos(geomang);
Y=R.*sin(geomang);

dXdtheta=dRdtheta.*cos(geomang)-R.*dgeomangdtheta.*sin(geomang);
dXdzeta =dRdzeta .*cos(geomang)-R.*dgeomangdzeta .*sin(geomang);
dYdtheta=dRdtheta.*sin(geomang)+R.*dgeomangdtheta.*cos(geomang);
dYdzeta =dRdzeta .*sin(geomang)+R.*dgeomangdzeta .*cos(geomang);

d2Xdtheta2=d2Rdtheta2.*cos(geomang) ...
    -2*dRdtheta.*dgeomangdtheta.*sin(geomang)...
    -R.*d2geomangdtheta2.*sin(geomang) ...
    -R.*dgeomangdtheta.^2.*cos(geomang);
d2Xdthetadzeta=d2Rdthetadzeta.*cos(geomang) ...
    -(dRdtheta.*dgeomangdzeta+dRdzeta.*dgeomangdtheta).*sin(geomang) ...
    -R.*d2geomangdthetadzeta.*sin(geomang) ...
    -R.*dgeomangdtheta.*dgeomangdzeta.*cos(geomang);
d2Xdzeta2=d2Rdzeta2.*cos(geomang) ...
    -2*dRdzeta.*dgeomangdzeta.*sin(geomang)...
    -R.*d2geomangdzeta2.*sin(geomang) ...
    -R.*dgeomangdzeta.^2.*cos(geomang);

d2Ydtheta2=d2Rdtheta2.*sin(geomang) ...
    +2*dRdtheta.*dgeomangdtheta.*cos(geomang)...
    +R.*d2geomangdtheta2.*cos(geomang) ...
    -R.*dgeomangdtheta.^2.*sin(geomang);
d2Ydthetadzeta=d2Rdthetadzeta.*sin(geomang) ...
    +(dRdtheta.*dgeomangdzeta+dRdzeta.*dgeomangdtheta).*cos(geomang) ...
    +R.*d2geomangdthetadzeta.*cos(geomang) ...
    -R.*dgeomangdtheta.*dgeomangdzeta.*sin(geomang);
d2Ydzeta2=d2Rdzeta2.*sin(geomang) ...
    +2*dRdzeta.*dgeomangdzeta.*cos(geomang)...
    +R.*d2geomangdzeta2.*cos(geomang) ...
    -R.*dgeomangdzeta.^2.*sin(geomang);

g_thetatheta=dXdtheta.^2+dYdtheta.^2+dZdtheta.^2;
g_zetazeta  =dXdzeta.^2 +dYdzeta.^2 +dZdzeta.^2;
g_thetazeta =dXdtheta.*dXdzeta+dYdtheta.*dYdzeta+dZdtheta.*dZdzeta;

gradpsi.X=B.^2./(G+iota*I).*(dYdtheta.*dZdzeta-dZdtheta.*dYdzeta);
gradpsi.Y=B.^2./(G+iota*I).*(dZdtheta.*dXdzeta-dXdtheta.*dZdzeta);
gradpsi.Z=B.^2./(G+iota*I).*(dXdtheta.*dYdzeta-dYdtheta.*dXdzeta);
gpsipsi=gradpsi.X.^2+gradpsi.Y.^2+gradpsi.Z.^2;

Booz.XYZ.r=zeros(3,Ntheta,Nzeta);
Booz.XYZ.r(1,:,:)=X;
Booz.XYZ.r(2,:,:)=Y;
Booz.XYZ.r(3,:,:)=Z;

Booz.XYZ.gradpsi=zeros(3,Ntheta,Nzeta);
Booz.XYZ.gradpsi(1,:,:)=gradpsi.X;
Booz.XYZ.gradpsi(2,:,:)=gradpsi.Y;
Booz.XYZ.gradpsi(3,:,:)=gradpsi.Z;

Booz.XYZ.e_theta=zeros(3,Ntheta,Nzeta);
Booz.XYZ.e_theta(1,:,:)=dXdtheta;
Booz.XYZ.e_theta(2,:,:)=dYdtheta;
Booz.XYZ.e_theta(3,:,:)=dZdtheta;

Booz.XYZ.e_zeta=zeros(3,Ntheta,Nzeta);
Booz.XYZ.e_zeta(1,:,:)=dXdzeta;
Booz.XYZ.e_zeta(2,:,:)=dYdzeta;
Booz.XYZ.e_zeta(3,:,:)=dZdzeta;

Booz.XYZ.d2rdtheta2=zeros(3,Ntheta,Nzeta);
Booz.XYZ.d2rdtheta2(1,:,:)=d2Xdtheta2;
Booz.XYZ.d2rdtheta2(2,:,:)=d2Ydtheta2;
Booz.XYZ.d2rdtheta2(3,:,:)=d2Zdtheta2;

Booz.XYZ.d2rdzeta2=zeros(3,Ntheta,Nzeta);
Booz.XYZ.d2rdzeta2(1,:,:)=d2Xdzeta2;
Booz.XYZ.d2rdzeta2(2,:,:)=d2Ydzeta2;
Booz.XYZ.d2rdzeta2(3,:,:)=d2Zdzeta2;

Booz.XYZ.d2rdthetadzeta=zeros(3,Ntheta,Nzeta);
Booz.XYZ.d2rdthetadzeta(1,:,:)=d2Xdthetadzeta;
Booz.XYZ.d2rdthetadzeta(2,:,:)=d2Ydthetadzeta;
Booz.XYZ.d2rdthetadzeta(3,:,:)=d2Zdthetadzeta;

CX=(d2Xdzeta2+2*iota*d2Xdthetadzeta+iota^2*d2Xdtheta2).*(B.^2/(G+iota*I)).^2;
CY=(d2Ydzeta2+2*iota*d2Ydthetadzeta+iota^2*d2Ydtheta2).*(B.^2/(G+iota*I)).^2;
CZ=(d2Zdzeta2+2*iota*d2Zdthetadzeta+iota^2*d2Zdtheta2).*(B.^2/(G+iota*I)).^2;

BdotgradabsB=B.^2/(G+iota*I).*(iota*dBdtheta+dBdzeta);
BxgradpsidotgradabsB=B.^2/(G+iota*I).*(G.*dBdtheta-I.*dBdzeta);
Booz.BdotgradabsB=BdotgradabsB;
Booz.BxgradpsidotgradabsB=BxgradpsidotgradabsB;

if 0 %Double-check the FFT results against the slow method
  fig(1+useFFT*4)
  surf(dYdzeta);shading flat;view(0,90)
  fig(2+useFFT*4)
  surf(dXdzeta);shading flat;view(0,90)
  fig(3+useFFT*4)
  surf(dYdtheta);shading flat;view(0,90)
  fig(4+useFFT*4)
  surf(dXdtheta);shading flat;view(0,90)
  error('stop stop!')
end

%%%%% Direction of B
if 0
  alt=1;
  %I thought this was wrong, but it is right
  BR=(dRdzeta+iota*dRdtheta).*B.^2/(G+iota*I);
  BZ=(dZdzeta+iota*dZdtheta).*B.^2/(G+iota*I);
  Bgeomang=(R.*dgeomangdzeta+iota*R.*dgeomangdtheta).*B.^2/(G+iota*I);
  BX=BR.*cos(geomang)-Bgeomang.*sin(geomang);
  BY=BR.*sin(geomang)+Bgeomang.*cos(geomang);
  Booz.XYZ.B=zeros(3,Ntheta,Nzeta);
  Booz.XYZ.B(1,:,:)=BX;
  Booz.XYZ.B(2,:,:)=BY;
  Booz.XYZ.B(3,:,:)=BZ;
else 
  alt=2;
  Booz.XYZ.B=zeros(3,Ntheta,Nzeta);
  Booz.XYZ.B(1,:,:)=...
      squeeze(Booz.XYZ.e_zeta(1,:,:)+iota*Booz.XYZ.e_theta(1,:,:))...
      .*B.^2/(G+iota*I);
  Booz.XYZ.B(2,:,:)=...
      squeeze(Booz.XYZ.e_zeta(2,:,:)+iota*Booz.XYZ.e_theta(2,:,:))...
      .*B.^2/(G+iota*I);
  Booz.XYZ.B(3,:,:)=...
      squeeze(Booz.XYZ.e_zeta(3,:,:)+iota*Booz.XYZ.e_theta(3,:,:))...
      .*B.^2/(G+iota*I);
  
  BX=squeeze(Booz.XYZ.B(1,:,:));
  BY=squeeze(Booz.XYZ.B(2,:,:));
  BZ=squeeze(Booz.XYZ.B(3,:,:));
  BR=BX.*cos(geomang)+BY.*sin(geomang);
  Bgeomang=-BX.*sin(geomang)+BY.*cos(geomang);
end

if 0 %testing
  BRopt=(dRdzeta+iota*dRdtheta).*B.^2/(G+iota*I);
  BZopt=(dZdzeta+iota*dZdtheta).*B.^2/(G+iota*I);
  Bgeomangopt=(R.*dgeomangdzeta+iota*R.*dgeomangdtheta).*B.^2/(G+iota*I);
  BXopt=BRopt.*cos(geomang)-Bgeomangopt.*sin(geomang);
  BYopt=BRopt.*sin(geomang)+Bgeomangopt.*cos(geomang);

  fig(50)
  surf(BXopt.^2+BYopt.^2+BZopt.^2-B.^2)
  %surf(BYopt-BY)
  %surf((dRdzeta+iota*dRdtheta).*B.^2/(G+iota*I).*cos(geomang)...
  %     -(R.*dgeomangdzeta+iota*R.*dgeomangdtheta).*B.^2/(G+iota*I).*sin(geomang)-BX)
  fig(51)
  %surf(BZopt-BZ)
  %surf((dZdzeta+iota*dZdtheta).*B.^2/(G+iota*I)-BZ)
  %surf((dRdzeta+iota*dRdtheta).*B.^2/(G+iota*I).*sin(geomang)...
  %     +(R.*dgeomangdzeta+iota*R.*dgeomangdtheta).*B.^2/(G+iota*I).*cos(geomang)-BY)
  surf(BX.^2+BY.^2+BZ.^2-B.^2)
end


Booz.XYZ.Bxgradpsi=cross(Booz.XYZ.B,Booz.XYZ.gradpsi);

%% Calculate the curvature
Booz.XYZ.curv=zeros(3,Ntheta,Nzeta);
Booz.XYZ.curv(1,:,:)=(CX+BX./B.*BdotgradabsB)./B.^2;
Booz.XYZ.curv(2,:,:)=(CY+BY./B.*BdotgradabsB)./B.^2;
Booz.XYZ.curv(3,:,:)=(CZ+BZ./B.*BdotgradabsB)./B.^2;

Booz.curv_normal=squeeze(dot(Booz.XYZ.gradpsi,Booz.XYZ.curv))./sqrt(gpsipsi);
Booz.curv_geodes1=BxgradpsidotgradabsB./sqrt(gpsipsi)./B.^2;
Booz.curv_geodes2=squeeze(dot(cross(Booz.XYZ.B,Booz.XYZ.gradpsi),Booz.XYZ.curv))./sqrt(gpsipsi)./B;
Booz.gradpsidotgradB=Booz.curv_normal.*sqrt(gpsipsi).*B.^2-mu0dpdpsi.*gpsipsi;

% The following is for equilibria close to axi-symmetry, where we want to calculate 
% the deviation from axi-symmetry in terms of dBinsurf and dBperp.
% Approximating that dr/dzeta is roughly in the geometrical toroidal direction I 
% assume that poloidal components of B in the "closest" axisymmetric equilibrium can
% be obtained by averaging BR and BZ over the Boozer toroidal angle. Then we project
% the deviation from axi-symmetry on the two directions perpedicular to the toroidal direction.
BRmean=mean(BR')'*ones(1,size(dRdtheta,2)); %mean over the Boozer toroidal angle
BZmean=mean(BZ')'*ones(1,size(dRdtheta,2)); %mean over the Boozer toroidal angle

bRmean=BRmean./sqrt(BRmean.^2+BZmean.^2);  %unit vector tangent to flux surface
bZmean=BZmean./sqrt(BRmean.^2+BZmean.^2);  %in the poloidal plane

dBR=BR-BRmean;
dBZ=BZ-BZmean;
Booz.dBinsurf=dBR.*bRmean+dBZ.*bZmean; %dot product 
Booz.dBperp=dBR.*bZmean-dBZ.*bRmean;   %cross product 
%sqrt((dBR-Booz.dBinsurf*bRmean).^2+(dBZ-Booz.dBinsurf*bZmean).^2);


Booz.Dzetacylphi=Dzetacylphi;     %for later use
Booz.dBdtheta=dBdtheta;
Booz.dBdzeta=dBdzeta;
Booz.dRdtheta=dRdtheta;
Booz.dRdzeta=dRdzeta;
Booz.dZdtheta=dZdtheta;
Booz.dZdzeta=dZdzeta;
Booz.BR=BR;%wrong:(dRdzeta+iota*dRdtheta).*B.^2/(G+iota*I);
Booz.BZ=BZ;%(dZdzeta+iota*dZdtheta).*B.^2/(G+iota*I);
Booz.Bgeomang=Bgeomang;%wrong:(R.*dgeomangdzeta+iota*R.*dgeomangdtheta).*B.^2/(G+iota*I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if cylindrical coordinates are possible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rinboard=zeros(1,Nzeta);
Routboard=zeros(1,Nzeta);
cylPossible=1;
for zind=1:Nzeta
  Zv=Z(:,zind);
  inds=find(diff(sign(Zv+eps))~=0);
  if length(inds)==1
    inds=[inds;Ntheta];
  elseif length(inds)==3
    Zv
    inds
    error('not implemented case')
  elseif isempty(inds)
    cylPossible=0;
  end
  if cylPossible
    Zvf=[Zv;Zv(1)];
    Rvf=[R(:,zind);R(1,zind)];
    R0s(1)=interp1(Zvf(inds(1):inds(1)+1),Rvf(inds(1):inds(1)+1),0);
    R0s(2)=interp1(Zvf(inds(2):inds(2)+1),Rvf(inds(2):inds(2)+1),0);
    Rinboard(zind)=min(R0s);
    Routboard(zind)=max(R0s);
  end
end  
if cylPossible
  minRoutboard=min(Routboard);
  maxRinboard=max(Rinboard);
  cylPossible=(maxRinboard<minRoutboard);
end
Booz.cylPossible=cylPossible;
Ham.cylPossible=cylPossible;
Pest.cylPossible=cylPossible;

if cylPossible
  cylR00=(maxRinboard+minRoutboard)/2;
  cylth=mod(atan2(Z,R-cylR00),2*pi); %cylth is the geometrical toroidal angle.
  Booz.Dthetacylth=theta-cylth;               %for later use
  high=find(Booz.Dthetacylth>pi);
  Booz.Dthetacylth(high)=Booz.Dthetacylth(high)-2*pi;
  low=find(Booz.Dthetacylth<-pi);
  Booz.Dthetacylth(low)=Booz.Dthetacylth(low)+2*pi;
  cylth=theta-Booz.Dthetacylth;
  
  cylr=sqrt(Z.^2+(R-cylR00).^2);
end

% ---------------------------------------------------------------------------------------
% Calculate parallel current u from harmonics of 1/B^2. Used in NTV calculation.
% \nabla_\parallel u = (2/B^4) \nabla B \times \vector{B} \cdot \iota \nabla \psi 
% ---------------------------------------------------------------------------------------
u = zeros(Ntheta,Nzeta); 
Dzetaphi        =zeros(Ntheta,Nzeta); %difference between zeta and phi (tor Hamada coord) 
dDzetaphi_dtheta=zeros(Ntheta,Nzeta);
dDzetaphi_dzeta =zeros(Ntheta,Nzeta);
%dudtheta = zeros(Ntheta,Nzeta);
%dudzeta = zeros(Ntheta,Nzeta);
h=1./(B.^2);
VPrimeHat=sum(sum(h))*4*pi^2/(Nzeta*Ntheta); %Note: VPrime=VPrimeHat*(G+iota*I)
FSAB2=4*pi^2/VPrimeHat;
%h00=VPrimeHat/(4*pi^2);
if useFFT
  %tic
  if 1
    [u,mnmats.u]=calcu(fftmn(h),G,I,iota,NPeriods);
  else
    hmn=unmnmat(fftmn(h));
    umn.c=iota*(G*hmn.m + I*hmn.n * NPeriods)./(hmn.n * NPeriods - iota*hmn.m) .* hmn.c;
    umn.c(hmn.m0ind,hmn.n0ind)=0;
    umn.s=iota*(G*hmn.m + I*hmn.n * NPeriods)./(hmn.n * NPeriods - iota*hmn.m) .* hmn.s;
    umn.s(hmn.m0ind,hmn.n0ind)=NaN; %Important to set. Indicates M odd/even
    umn.m=hmn.m;umn.n=hmn.n;
    umn.m0ind=hmn.m0ind;umn.n0ind=hmn.n0ind;
    mnmats.u=mnmat(umn);
    u=ifftmn(mnmats.u);
  end
  
  if 1
    mnmats.Dzetaphi=invJacBdotgrad(fftmn(1-h*FSAB2),iota,NPeriods);
    Dzetaphi=ifftmn(mnmats.Dzetaphi);
  else
    hmn=unmnmat(fftmn(h));
    umn=unmnmat(fftmn(u));
    Dzetaphi_mn.s  =  FSAB2/(G+iota*I)./hmn.m.*(umn.c/iota-I*hmn.c);     %for m~=0
    Dzetaphi_mn.c  = -FSAB2/(G+iota*I)./hmn.m.*(umn.s/iota-I*hmn.s);    %for m~=0
    Dzetaphi2_mn.s =  FSAB2/(G+iota*I)./hmn.n/NPeriods.*(umn.c+G*hmn.c); %for n~=0
    Dzetaphi2_mn.c = -FSAB2/(G+iota*I)./hmn.n/NPeriods.*(umn.s+G*hmn.s); %for n~=0
    Dzetaphi_mn.s(hmn.m0ind,:)=Dzetaphi2_mn.s(hmn.m0ind,:);
    Dzetaphi_mn.c(hmn.m0ind,:)=Dzetaphi2_mn.c(hmn.m0ind,:);
    
    Dzetaphi_mn.s(hmn.m0ind,hmn.n0ind)=NaN;
    Dzetaphi_mn.c(hmn.m0ind,hmn.n0ind)=0;
    if not(mod(Nzeta,2))
      Dzetaphi_mn.c(hmn.m0ind,end)=0; %info was lost here
      Dzetaphi_mn.s(hmn.m0ind,end)=NaN; %Important to set. Indicates M odd/even
      if not(mod(Ntheta,2))
        Dzetaphi_mn.c(end,end)=0;       %info was lost here
        Dzetaphi_mn.s(end,end)=NaN;       %info was lost here
      end
    end
    if not(mod(Ntheta,2))
      Dzetaphi_mn.c(end,hmn.n0ind)=0; %info was lost here
      Dzetaphi_mn.s(end,hmn.n0ind)=NaN; %info was lost here
    end
    Dzetaphi_mn.m=hmn.m;Dzetaphi_mn.n=hmn.n;
    Dzetaphi_mn.m0ind=hmn.m0ind;Dzetaphi_mn.n0ind=hmn.n0ind;
    Dzetaphi_mn=mnmat(Dzetaphi_mn);
    Dzetaphi=ifftmn(Dzetaphi_mn);
  end
  [dDzetaphi_dtheta,dDzetaphi_dzeta]=ifftmn(grad(mnmats.Dzetaphi,NPeriods));
  %toc
else %not(useFFT)
  %tic
    if not(Geom.StelSym) %sine components exist
      for m=0:floor(Ntheta/2) %Nyquist max freq.
        if m==0 || m==Ntheta/2
          nrange=nvecL0c;%=0:floor(Nzeta/2);
          c3D=cos(m * theta3D0c - n3D0c * NPeriods .* zeta3D0c);
          s3D=sin(m * theta3D0c - n3D0c * NPeriods .* zeta3D0c);
        else
          nrange=nvecL;% =-floor(Nzeta/2)+Nzeta_even:floor(Nzeta/2);
          c3D=cos(m * theta3D - n3D * NPeriods .* zeta3D);
          s3D=sin(m * theta3D - n3D * NPeriods .* zeta3D);
        end
        for n=nrange
          nind=find(nrange==n);  
          c=squeeze(c3D(:,:,nind)); %FAST
          s=squeeze(s3D(:,:,nind));
          %c=cos(m * theta - n * NPeriods * zeta); %SLOW
          %s=sin(m * theta - n * NPeriods * zeta);
          
          %cos
          if (m==0 && n==0) 
            %We know that u_00=0 by definition!
            hmnc=0;umnc=0;Dzetaphi1_mns=0;Dzetaphi2_mns=0;
          else
            iscorner=((m==0 && n==0) || ...
                      (m==Ntheta/2 && n==0) || ...
                      (m==Ntheta/2 && n==Nzeta/2) || ...
                      (m==0 && n==Nzeta/2));
            hmnc = (2-iscorner)/(Ntheta*Nzeta) *...
                   sum(sum(c.*h));
            umnc = ...
                iota*(G*m + I*n * NPeriods)/(n * NPeriods - iota*m) * hmnc;
            %There are two ways of calculating Dzetaphi. We need both.
            Dzetaphi1_mns = ...
                FSAB2/(G+iota*I)/m*(umnc/iota-I*hmnc);     %for m~=0
            Dzetaphi2_mns = ...
                FSAB2/(G+iota*I)/n/NPeriods*(umnc+G*hmnc); %for n~=0
          end
          
          %sin
          if (m==0 && n==0) || (m==0 && n==Nzeta/2) ...
                || (m==Ntheta/2 && n==0) || (m==Ntheta/2 && n==Nzeta/2)
            hmns=0;umns=0;Dzetaphi1_mnc=0;Dzetaphi2_mnc=0;
          else
            hmns = 2/(Ntheta*Nzeta) *...
                   sum(sum(s.*h));
            umns = ...
                iota*(G*m + I*n * NPeriods)/(n * NPeriods - iota*m) * hmns;
            Dzetaphi1_mnc = ...  
                -FSAB2/(G+iota*I)/m*(umns/iota-I*hmns);     %for m~=0
            Dzetaphi2_mnc = ...
                -FSAB2/(G+iota*I)/n/NPeriods*(umns+G*hmns); %for n~=0
          end
          
          %assemble
          u = u + umnc * c + umns * s;
          if m>0
            Dzetaphi = Dzetaphi + Dzetaphi1_mnc * c + Dzetaphi1_mns * s; 
            dDzetaphi_dtheta = dDzetaphi_dtheta + m*(-Dzetaphi1_mnc*s+Dzetaphi1_mns*c);
            dDzetaphi_dzeta  = dDzetaphi_dzeta  - n*NPeriods*(-Dzetaphi1_mnc*s+Dzetaphi1_mns*c);
          else
            Dzetaphi = Dzetaphi + Dzetaphi2_mnc * c + Dzetaphi2_mns * s; 
            dDzetaphi_dzeta  = dDzetaphi_dzeta  - n*NPeriods*(-Dzetaphi2_mnc*s+Dzetaphi2_mns*c);            
          end
          %dudtheta = dudtheta ...
          %    + umns * m * c - umnc * m * s;
          %dudzeta = dudzeta ...
          %    + umnc * n * NPeriods * s - umns * n * NPeriods * c; 
        end
      end
    else %only cosinus components
      for m=0:floor(Ntheta/2) %Nyquist max freq.
        if m==0 || m==Ntheta/2
          nrange=nvecL0c;%=0:floor(Nzeta/2);
          c3D=cos(m * theta3D0c - n3D0c * NPeriods .* zeta3D0c);
          s3D=sin(m * theta3D0c - n3D0c * NPeriods .* zeta3D0c);
        else
          nrange=nvecL;% =-floor(Nzeta/2)+Nzeta_even:floor(Nzeta/2);
          c3D=cos(m * theta3D - n3D * NPeriods .* zeta3D);
          s3D=sin(m * theta3D - n3D * NPeriods .* zeta3D);
        end
        for n=nrange
          nind=find(nrange==n);  
          c=squeeze(c3D(:,:,nind)); %FAST
          s=squeeze(s3D(:,:,nind));
          %c=cos(m * theta - n * NPeriods * zeta); %SLOW
          %s=sin(m * theta - n * NPeriods * zeta);

          if (m==0 && n==0) 
            %We know that u_00=0 by definition!
            hmnc=0;umnc=0;Dzetaphi1_mns=0;Dzetaphi2_mns=0;
          else
            iscorner=((m==0 && n==0) || ...
                      (m==Ntheta/2 && n==0) || ...
                      (m==Ntheta/2 && n==Nzeta/2) || ...
                      (m==0 && n==Nzeta/2));
                      
            hmnc = (2-iscorner)/(Ntheta*Nzeta) *...
                   sum(sum(c.*h));
            umnc = ...
                iota*(G*m + I*n * NPeriods)/(n * NPeriods - iota*m) * hmnc;
            Dzetaphi1_mns = ...
                FSAB2/(G+iota*I)/m*(umnc/iota-I*hmnc);     %for m~=0
            Dzetaphi2_mns = ...
                FSAB2/(G+iota*I)/n/NPeriods*(umnc+G*hmnc); %for n~=0
          end
          
          %disp([num2str(m),', ',num2str(n),', ',num2str(umnc)]) %control output
          u = u + umnc * c;
          if m>0
            Dzetaphi = Dzetaphi + Dzetaphi1_mns * s;    
            dDzetaphi_dtheta = dDzetaphi_dtheta + m*Dzetaphi1_mns*c;
            dDzetaphi_dzeta  = dDzetaphi_dzeta  - n*NPeriods*Dzetaphi1_mns*c;
          else
            Dzetaphi = Dzetaphi + Dzetaphi2_mns * s; 
            dDzetaphi_dzeta  = dDzetaphi_dzeta  - n*NPeriods*Dzetaphi2_mns*c;            
          end
          %dudtheta = dudtheta ...
          %    - umnc * m * s;
          %dudzeta = dudzeta ...
          %          + umnc * n * NPeriods * s;   
        end              
      end
    end
  %toc
end

%[dDzetaphi_dtheta,dDzetaphi_dzeta]=ifftmn(grad(Dzetaphi_mn,NPeriods));
%old way:
%dDzetaphi_dtheta_old= FSAB2/(G+iota*I)/iota*(u-iota*I*(h-1/FSAB2));
%dDzetaphi_dzeta_old =-FSAB2/(G+iota*I) *    (u +    G*(h-1/FSAB2));
%dDzetaphi_dtheta-dDzetaphi_dtheta_old
%dDzetaphi_dzeta-dDzetaphi_dzeta_old

%B_psi_tilde is defined by B_psi=B_psi_00+B_psi_tilde
B_psi_tilde=-mu0dpdpsi/FSAB2*(G+iota*I)*Dzetaphi;
%B_psi_tilde_mn=-mu0dpdpsi/FSAB2*(G+iota*I)*Dzetaphi_mn;


XtozmXzot=dXdtheta.*dDzetaphi_dzeta-dXdzeta.*dDzetaphi_dtheta;
YtozmYzot=dYdtheta.*dDzetaphi_dzeta-dYdzeta.*dDzetaphi_dtheta;
ZtozmZzot=dZdtheta.*dDzetaphi_dzeta-dZdzeta.*dDzetaphi_dtheta;

dXdphi=B.^2/FSAB2.*(dXdzeta+iota*XtozmXzot);
dYdphi=B.^2/FSAB2.*(dYdzeta+iota*YtozmYzot);
dZdphi=B.^2/FSAB2.*(dZdzeta+iota*ZtozmZzot);

dXdvthet=B.^2/FSAB2.*(dXdtheta-XtozmXzot);
dYdvthet=B.^2/FSAB2.*(dYdtheta-YtozmYzot);
dZdvthet=B.^2/FSAB2.*(dZdtheta-ZtozmZzot);

g_phiphi    =dXdphi.^2  +dYdphi.^2  +dZdphi.^2;
g_vthetvthet=dXdvthet.^2+dYdvthet.^2+dZdvthet.^2;
g_vthetphi  =dXdphi.*dXdvthet+dYdphi.*dYdvthet+dZdphi.*dZdvthet;

%%%%%%%%%%%%%%%%% calculate B_psi_00, grad(theta) and grad(zeta) %%%%%%%%%%%%%%%%%%%%%
% Is not possible, don't do it!

if 0
d1=1./(dXdtheta.*dYdzeta-dXdzeta.*dYdtheta)

atheta_Y0=-dXdzeta.*d1;
atheta_YZ=-(dXdtheta.*dZdzeta-dXdzeta.*dZdtheta).*d1;
atheta_X0=dYdzeta.*d1;
atheta_XZ=(dYdtheta.*dZdzeta-dYdzeta.*dZdtheta).*d1;

azeta_Y0=dXdtheta.*d1;
azeta_YZ=-(dXdtheta.*dZdzeta-dXdzeta.*dZdtheta).*d1;
azeta_X0=-dYdtheta.*d1;
azeta_XZ=(dYdtheta.*dZdzeta-dYdzeta.*dZdtheta).*d1;

DX=BX-B_psi_tilde.*gradpsi.X
DY=BY-B_psi_tilde.*gradpsi.Y
DZ=BZ-B_psi_tilde.*gradpsi.Z

atheta_YZ_c = atheta_YZ - gradpsi.Y./gradpsi.Z
azeta_YZ_c  = azeta_YZ  - gradpsi.Y./gradpsi.Z
atheta_XZ_c = atheta_XZ - gradpsi.X./gradpsi.Z
azeta_XZ_c  = azeta_XZ  - gradpsi.X./gradpsi.Z

azeta_XZ_c .* atheta_YZ_c 
azeta_YZ_c .* atheta_XZ_c

d2=1/G/I./(azeta_XZ_c .* atheta_YZ_c - azeta_YZ_c .* atheta_XZ_c);

EX=DX-DZ.*gradpsi.X./gradpsi.Z-G*azeta_X0-I*atheta_X0;
EY=DY-DZ.*gradpsi.Y./gradpsi.Z-G*azeta_Y0-I*atheta_Y0;

EX./atheta_XZ_c
EY./atheta_YZ_c

gradzeta.Z =I*d2.*( EX.*atheta_YZ_c -EY.*atheta_XZ_c);
gradtheta.Z=G*d2.*(-EX.*azeta_YZ_c  +EY.*azeta_XZ_c );

gradzeta.X = azeta_X0 + azeta_XZ.*gradzeta.Z;
gradzeta.Y = azeta_Y0 + azeta_YZ.*gradzeta.Z;

gradtheta.X = atheta_X0 + atheta_XZ.*gradtheta.Z;
gradtheta.Y = atheta_Y0 + atheta_YZ.*gradtheta.Z;

B_psi_00 = (DZ-G*gradzeta.Z-I*gradtheta.Z)./gradpsi.Z;
error('stoppppp')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store variables in the Booz struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Booz.Nperiods=Geom.Nperiods;
Booz.Nzeta=Nzeta;
Booz.Ntheta=Ntheta;
Booz.Dzeta=Dzeta;
Booz.Dtheta=Dtheta;
Booz.zeta=zeta;
Booz.theta=theta;

Booz.cylR=R;
Booz.cylZ=Z;
Booz.cylphi=cylphi;
Booz.R00=mean(mean(R));
if cylPossible
  Booz.cylth=cylth;
  Booz.cylr=cylr;
  Booz.cylR00=cylR00;
end
Booz.X=X;
Booz.Y=Y;
Booz.iota=Geom.iota(rind);
Booz.G=Geom.Bphi(rind);
Booz.I=Geom.Btheta(rind);
Booz.mu0dpdpsi=mu0dpdpsi;
Booz.B=B;
Booz.u=u;
Booz.h=h;
Booz.Bpsitilde=B_psi_tilde;
Booz.dBpsidtheta=-Booz.mu0dpdpsi/iota.*(u-iota*I*(B.^-2 - 1/FSAB2));
Booz.dBpsidzeta = Booz.mu0dpdpsi/iota.*(u+G*(B.^-2 - 1/FSAB2));
Booz.FSAB2=FSAB2;
Booz.FSAu2B2=sum(sum(u.^2.*B.^2.*h))/sum(sum(h));
Booz.Jacob_psi_theta_zeta=Booz.h*(G+iota*I);
Booz.g_thetatheta=g_thetatheta;
Booz.g_thetazeta =g_thetazeta;
Booz.g_zetazeta  =g_zetazeta;
Booz.g_phiphi=g_phiphi;
Booz.g_vthetvthet=g_vthetvthet;
Booz.g_vthetphi=g_vthetphi;
Booz.gpsipsi=gpsipsi;
Booz.FSAg_phiphi=sum(sum(g_phiphi.*h))/sum(sum(h));
Booz.FSAgpsipsi=sum(sum(gpsipsi.*h))/sum(sum(h));


Booz.dpsidr_eff=sqrt(Booz.FSAgpsipsi);
Booz.r_eff=sqrt(Booz.FSAg_phiphi*Booz.FSAgpsipsi)/abs(Booz.G);
Booz.r_eff_appr1=sqrt(sum(sum(gpsipsi.*h./B.^2))/sum(sum(h)));
Booz.r_eff_appr2=sqrt((Booz.FSAg_phiphi-Booz.FSAu2B2-G^2/FSAB2)/iota^2);
%r_eff_appr3=sqrt(mean(mean(cylr.^2)))

%Copy data from Geom struct
Booz.B00=Geom.B00(rind);
Booz.m=Geom.m{rind};
Booz.n=Geom.n{rind};
Booz.Bmn=Geom.Bmn{rind};
Booz.StelSym=Geom.StelSym;
if Geom.StelSym %sine components exist
  Booz.parity=parity;
else
  Booz.parity=Geom.parity{rind};
end

if useFFT
  mnmats.gradpsidotgradB=fftmn(Booz.gradpsidotgradB);
  mnmats.dBpsidtheta=fftmn(Booz.dBpsidtheta);
  mnmats.dBpsidzeta=fftmn(Booz.dBpsidzeta);
  mnmats.h=fftmn(Booz.h);
  Booz.mnmat=mnmats;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate to a phi vthet grid 
% Here, phi is the Hamada tor. coord. and vthet (\vartheta) is the Hamada poloidal coordinate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
  [Booz,Ham] =interp2straightfieldlinecoords(Booz,Ham,Dzetaphi,'vthet','phi');
  [Booz,Pest]=interp2straightfieldlinecoords(Booz,Pest,Dzetacylphi,'ptheta','pzeta');
  Ham.Dzetapzeta=interp2_cyclic(Booz.theta,Booz.zeta,Booz.Dzetapzeta,Ham.theta,Ham.zeta,NPeriods);
  Ham.Dphipzeta=Ham.Dzetapzeta-Ham.Dzetaphi;
  Ham.pzeta=Ham.phi+Ham.Dzetaphi;
  Ham.ptheta=Ham.vthet+iota*Ham.Dzetaphi;
  Pest.Dzetaphi=interp2_cyclic(Booz.theta,Booz.zeta,Booz.Dzetaphi,Pest.theta,Pest.zeta,NPeriods);
  Pest.Dpzetaphi=Pest.Dzetaphi-Pest.Dzetapzeta;
  Pest.phi=Pest.pzeta-Pest.Dpzetaphi;
  Pest.vthet=Pest.ptheta-iota*Pest.Dpzetaphi;
else %This whole chunk has been moved to the above function interp2straightfieldlinecoords
  Booz.Dzetaphi=Dzetaphi;
  Booz.phi=zeta-Booz.Dzetaphi;
  Booz.vthet=theta-iota*Booz.Dzetaphi;

  %Choose the same resolution. Do not change this. It means 3D tensors can be used again 
  Nphi=Nzeta;
  Nvthet=Ntheta;
  Ham.Nperiods=Geom.Nperiods;
  Ham.Nphi=Nphi;
  Ham.Nvthet=Nvthet;

  %Hamada discretisation
  Ham.Dphi=2*pi/Nphi/NPeriods;
  Ham.Dvthet=2*pi/Nvthet;
  phivec=(0:Nphi-1)*Ham.Dphi;
  vthetvec=(0:Nvthet-1)'*Ham.Dvthet;
  [vthet, phi] = ndgrid(vthetvec,phivec);
  Ham.phi=phi;
  Ham.vthet=vthet;

  if 1 %new method
    Nleftadd=ceil(max(Booz.phi(:,1))/(2*pi/NPeriods));
    Nrightadd=ceil(1-min(Booz.phi(:,end))/(2*pi/NPeriods));
    BoozphiBig      =NaN*zeros(Ntheta*3,Nphi*(Nleftadd+1+Nrightadd));
    BoozvthetBig    =NaN*zeros(Ntheta*3,Nphi*(Nleftadd+1+Nrightadd));
    BoozDzetaphiBig=NaN*zeros(Ntheta*3,Nphi*(Nleftadd+1+Nrightadd));

    for addi=1:Nleftadd+1+Nrightadd
      addci=addi-(Nleftadd+1);
      toaddphi=2*pi/NPeriods*addci;
      BoozphiBig(1:Ntheta*3,(addi-1)*Nphi+1:addi*Nphi)=...
          [Booz.phi+toaddphi;
           Booz.phi+toaddphi;
           Booz.phi+toaddphi];
      BoozvthetBig(1:Ntheta*3,(addi-1)*Nphi+1:addi*Nphi)=...
          [Booz.vthet-2*pi;
           Booz.vthet;
           Booz.vthet+2*pi];
      BoozDzetaphiBig(1:Ntheta*3,(addi-1)*Nphi+1:addi*Nphi)=...
          [Booz.Dzetaphi;
           Booz.Dzetaphi;
           Booz.Dzetaphi];
    end
  else %old method
    BoozphiBig=[Booz.phi-2*pi/NPeriods, Booz.phi, Booz.phi+2*pi/NPeriods;
                Booz.phi-2*pi/NPeriods, Booz.phi, Booz.phi+2*pi/NPeriods;
                Booz.phi-2*pi/NPeriods, Booz.phi, Booz.phi+2*pi/NPeriods];
    
    BoozvthetBig=[Booz.vthet-2*pi,Booz.vthet-2*pi,Booz.vthet-2*pi;
                  Booz.vthet,     Booz.vthet,     Booz.vthet;
                  Booz.vthet+2*pi,Booz.vthet+2*pi,Booz.vthet+2*pi];

    BoozDzetaphiBig=[Booz.Dzetaphi,Booz.Dzetaphi,Booz.Dzetaphi;
                     Booz.Dzetaphi,Booz.Dzetaphi,Booz.Dzetaphi;
                     Booz.Dzetaphi,Booz.Dzetaphi,Booz.Dzetaphi];
  end

  %The following is the most time consuming step
  %tic
  Ham.Dzetaphi=griddata(BoozphiBig,BoozvthetBig,BoozDzetaphiBig,phi,vthet);
  %toc 

  Ham.zeta =Ham.phi+Ham.Dzetaphi;
  Ham.theta=Ham.vthet+iota*Ham.Dzetaphi;

  zetaDirection=sign(Booz.zeta(1,end)-Booz.zeta(1,1));
  thetaDirection=sign(Booz.theta(end,1)-Booz.theta(1,1));

  Boozzetabig=[Booz.zeta,Booz.zeta(:,1)+zetaDirection*2*pi/NPeriods];
  Boozzetabig=[Boozzetabig;Boozzetabig(1,:)];
  Boozthetabig=[Booz.theta;Booz.theta(1,:)+thetaDirection*2*pi];
  Boozthetabig=[Boozthetabig,Boozthetabig(:,1)];
  BoozBbig=[Booz.B,Booz.B(:,1)];
  BoozBbig=[BoozBbig;BoozBbig(1,:)];

  %We have two options on how to interpolate. 
  %We use interp2, because it is much faster than griddata

  % major radius R
  BoozRbig=[Booz.cylR,Booz.cylR(:,1)];
  BoozRbig=[BoozRbig;BoozRbig(1,:)];
  Ham.cylR=interp2(Boozzetabig,Boozthetabig,BoozRbig,...
                   mod(Ham.zeta,2*pi/NPeriods),mod(Ham.theta,2*pi)); 
  % vertical coordinate Z
  BoozZbig=[Booz.cylZ,Booz.cylZ(:,1)];
  BoozZbig=[BoozZbig;BoozZbig(1,:)];
  Ham.cylZ=interp2(Boozzetabig,Boozthetabig,BoozZbig,...
                   mod(Ham.zeta,2*pi/NPeriods),mod(Ham.theta,2*pi)); 
  % Geometrical toroidal angle (cylphi)
  Boozcylphibig=[Booz.cylphi,Booz.cylphi(:,1)+zetaDirection*2*pi/NPeriods];
  Boozcylphibig=[Boozcylphibig;Boozcylphibig(1,:)];
  Ham.cylphi=interp2(Boozzetabig,Boozthetabig,Boozcylphibig,...
                     mod(Ham.zeta,2*pi/NPeriods),mod(Ham.theta,2*pi))+...
      Ham.zeta-mod(Ham.zeta,2*pi/NPeriods); 
  Ham.R00=mean(mean(Ham.cylR));
  if cylPossible
    % Geometrical poloidal angle (cylth)
    Boozcylthbig=[Booz.cylth;Booz.cylth(1,:)+thetaDirection*2*pi];
    Boozcylthbig=[Boozcylthbig,Boozcylthbig(:,1)];
    Ham.cylth=interp2(Boozzetabig,Boozthetabig,Boozcylthbig,...
                      mod(Ham.zeta,2*pi/NPeriods),mod(Ham.theta,2*pi))+...
              Ham.theta-mod(Ham.theta,2*pi); 

    Ham.cylr=sqrt(Ham.cylZ.^2+(Ham.cylR-cylR00).^2);

    Ham.cylR00=cylR00;
  end
  Ham.X=Ham.cylR.*cos(-Ham.cylphi);
  Ham.Y=Ham.cylR.*sin(-Ham.cylphi);

  %Option 1:
  %Ham.B=griddata(BoozphiBig,BoozvthetBig,BoozBBig,phi,vthet); %time=1 second
  %Option 2:
  Ham.B=interp2(Boozzetabig,Boozthetabig,BoozBbig,...
                mod(Ham.zeta,2*pi/NPeriods),mod(Ham.theta,2*pi)); %time=3 milliseconds

  % Interpolate the rest of the interesting stored data:
  % The parallel curent u:
  Boozubig=[Booz.u,Booz.u(:,1)];
  Boozubig=[Boozubig;Boozubig(1,:)];
  Ham.u=interp2(Boozzetabig,Boozthetabig,Boozubig,...
                mod(Ham.zeta,2*pi/NPeriods),mod(Ham.theta,2*pi));
  Ham.h=1./Ham.B.^2;
  Ham.iota=Booz.iota;
  Ham.G=Booz.G;
  Ham.I=Booz.I;
  Ham.FSAB2=Booz.FSAB2;
  Ham.FSAg_phiphi=Booz.FSAg_phiphi;
  Ham.FSAgpsipsi =Booz.FSAgpsipsi;
  Ham.FSAu2B2=Booz.FSAu2B2;
end
Ham.Jacob_psi_vthet_phi=(G+iota*I)/Ham.FSAB2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If possible, create a cylindrical discretization cylth,cylphi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(cylPossible)
  Cyl='Not possible to find R00!';
else
  Cyl.Nperiods=Geom.Nperiods;
  Cyl.Ncylphi=Nzeta;
  Cyl.Ncylth=Ntheta;
  Cyl.Dcylphi=Dzeta;
  Cyl.Dcylth=Dtheta;

  %Cyl discretisation
  cylphivec=(0:Cyl.Ncylphi-1)*Cyl.Dcylphi;
  cylthvec=(0:Cyl.Ncylth-1)'*Cyl.Dcylth;
  [cylth, cylphi] = ndgrid(cylthvec,cylphivec);
  Cyl.cylphi=cylphi;
  Cyl.cylth=cylth;

  BoozcylphiBig=[Booz.cylphi-2*pi/NPeriods, Booz.cylphi, Booz.cylphi+2*pi/NPeriods;
                 Booz.cylphi-2*pi/NPeriods, Booz.cylphi, Booz.cylphi+2*pi/NPeriods;
                 Booz.cylphi-2*pi/NPeriods, Booz.cylphi, Booz.cylphi+2*pi/NPeriods];        
  BoozcylthBig=[Booz.cylth-2*pi,Booz.cylth-2*pi,Booz.cylth-2*pi;
                Booz.cylth,     Booz.cylth,     Booz.cylth;
                Booz.cylth+2*pi,Booz.cylth+2*pi,Booz.cylth+2*pi];
  Booz.DthetacylthBig=[Booz.Dthetacylth,Booz.Dthetacylth,Booz.Dthetacylth;
                 Booz.Dthetacylth,Booz.Dthetacylth,Booz.Dthetacylth;
                 Booz.Dthetacylth,Booz.Dthetacylth,Booz.Dthetacylth];
  Booz.DzetacylphiBig=[Booz.Dzetacylphi,Booz.Dzetacylphi,Booz.Dzetacylphi;
                 Booz.Dzetacylphi,Booz.Dzetacylphi,Booz.Dzetacylphi;
                 Booz.Dzetacylphi,Booz.Dzetacylphi,Booz.Dzetacylphi];

  Cyl.Dthetacylth=griddata(BoozcylphiBig,BoozcylthBig,Booz.DthetacylthBig,Cyl.cylphi,Cyl.cylth);
  Cyl.Dzetacylphi=griddata(BoozcylphiBig,BoozcylthBig,Booz.DzetacylphiBig,Cyl.cylphi,Cyl.cylth);
  Cyl.theta = Cyl.cylth+Cyl.Dthetacylth;
  Cyl.zeta = Cyl.cylphi+Cyl.Dzetacylphi;

  %BoozzetaBig=[Booz.zeta-2*pi/NPeriods, Booz.zeta, Booz.zeta+2*pi/NPeriods;
  %             Booz.zeta-2*pi/NPeriods, Booz.zeta, Booz.zeta+2*pi/NPeriods;
  %             Booz.zeta-2*pi/NPeriods, Booz.zeta, Booz.zeta+2*pi/NPeriods];
  %BoozthetaBig=[Booz.theta-2*pi,Booz.theta-2*pi,Booz.theta-2*pi;
  %              Booz.theta,     Booz.theta,     Booz.theta;
  %              Booz.theta+2*pi,Booz.theta+2*pi,Booz.theta+2*pi];
  %Cyl.zeta =griddata(BoozcylphiBig,BoozcylthBig,BoozzetaBig,Cyl.cylphi,Cyl.cylth);
  %Cyl.theta=griddata(BoozcylphiBig,BoozcylthBig,BoozthetaBig,Cyl.cylphi,Cyl.cylth);

  zetaDirection=sign(Booz.zeta(1,end)-Booz.zeta(1,1));
  thetaDirection=sign(Booz.theta(end,1)-Booz.theta(1,1));
  % Cyl toroidal angle
  Boozzetabig=[Booz.zeta,Booz.zeta(:,1)+zetaDirection*2*pi/NPeriods];
  Boozzetabig=[Boozzetabig;Boozzetabig(1,:)];
  Boozthetabig=[Booz.theta;Booz.theta(1,:)+thetaDirection*2*pi];
  Boozthetabig=[Boozthetabig,Boozthetabig(:,1)];
  Boozphibig=[Booz.phi,Booz.phi(:,1)+zetaDirection*2*pi/NPeriods];
  Boozphibig=[Boozphibig;Boozphibig(1,:)];
  Cyl.phi=interp2(Boozzetabig,Boozthetabig,Boozphibig,...
                  mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi))+...
          Cyl.zeta-mod(Cyl.zeta,2*pi/NPeriods); 
  % Cyl poloidal angle
  Boozvthetbig=[Booz.vthet;Booz.vthet(1,:)+thetaDirection*2*pi];
  Boozvthetbig=[Boozvthetbig,Boozvthetbig(:,1)];
  Cyl.vthet=interp2(Boozzetabig,Boozthetabig,Boozvthetbig,...
                    mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi))+...
            Cyl.theta-mod(Cyl.theta,2*pi); 

  % major radius R
  BoozRbig=[Booz.cylR,Booz.cylR(:,1)];
  BoozRbig=[BoozRbig;BoozRbig(1,:)];
  Cyl.cylR=interp2(Boozzetabig,Boozthetabig,BoozRbig,...
                   mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi)); 
  % vertical coordinate Z
  BoozZbig=[Booz.cylZ,Booz.cylZ(:,1)];
  BoozZbig=[BoozZbig;BoozZbig(1,:)];
  Cyl.cylZ=interp2(Boozzetabig,Boozthetabig,BoozZbig,...
                   mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi)); 
  Cyl.cylr=sqrt(Cyl.cylZ.^2+(Cyl.cylR-cylR00).^2);
  Cyl.cylR00=cylR00;
  
  Cyl.X=Cyl.cylR.*cos(-Cyl.cylphi);
  Cyl.Y=Cyl.cylR.*sin(-Cyl.cylphi);
  
  BoozBbig=[Booz.B,Booz.B(:,1)];
  BoozBbig=[BoozBbig;BoozBbig(1,:)];
  Cyl.B=interp2(Boozzetabig,Boozthetabig,BoozBbig,...
                mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi)); 
  Cyl.B00=sum(sum(Cyl.B))/(Cyl.Ncylth*Cyl.Ncylphi);
  
  %B vector
  BoozBRbig=[Booz.BR,Booz.BR(:,1)];
  BoozBRbig=[BoozBRbig;BoozBRbig(1,:)];
  BoozBZbig=[Booz.BZ,Booz.BZ(:,1)];
  BoozBZbig=[BoozBZbig;BoozBZbig(1,:)];
  BoozBgeomangbig=[Booz.Bgeomang,Booz.Bgeomang(:,1)];
  BoozBgeomangbig=[BoozBgeomangbig;BoozBgeomangbig(1,:)];
  Cyl.BR=interp2(Boozzetabig,Boozthetabig,BoozBRbig,...
                mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi)); 
  Cyl.BZ=interp2(Boozzetabig,Boozthetabig,BoozBZbig,...
                mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi)); 
  Cyl.Bgeomang=interp2(Boozzetabig,Boozthetabig,BoozBgeomangbig,...
                mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi)); 
  BRmean=mean(Cyl.BR')'*ones(1,size(dRdtheta,2));
  BZmean=mean(Cyl.BZ')'*ones(1,size(dRdtheta,2));
  bRmean=BRmean./sqrt(BRmean.^2+BZmean.^2);
  bZmean=BZmean./sqrt(BRmean.^2+BZmean.^2);

  dBR=Cyl.BR-BRmean;
  dBZ=Cyl.BZ-BZmean;
  Cyl.dBinsurf=dBR.*bRmean+dBZ.*bZmean;
  Cyl.dBperp=dBR.*bZmean-dBZ.*bRmean;
  
  
  if 0
    if useFFT %calculate Bmn in Cyl coordinates
      precision=eps*1e3;
      min_Bmn=max(Geom.Bfilter.min_Bmn*Booz.B00/Cyl.B00,precision);
      CylBmn=mnlist(fftmn(Cyl.B),min_Bmn,'relative');
      Cyl.parity=CylBmn.cosparity;
      Cyl.Bmn=CylBmn.data;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculate Bmn in Hamada coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ham.B00=sum(sum(Ham.B))/(Ham.Nvthet*Ham.Nphi);

precision=eps*1e3;
believedAccuracyLoss=1;
min_Bmn=max(Geom.Bfilter.min_Bmn*Booz.B00/Ham.B00*believedAccuracyLoss,precision);

if useFFT
  %tic
  HamBmn=mnlist(fftmn(Ham.B),min_Bmn,'relative');
  HamRmn=mnlist(fftmn(Ham.cylR),min_Bmn,'relative');
  HamZmn=mnlist(fftmn(Ham.cylZ),min_Bmn,'relative');
  Ham.m=HamBmn.m;
  Ham.n=HamBmn.n;
  if Geom.StelSym
    if any(HamBmn.cosparity~=1)
      error(['Something has gone wrong. The calculated Hamada Bmn is not stellarator ', ...
             'symmetric!'])
    end
  else
    Ham.parity=HamBmn.cosparity;
  end  
  Ham.Bmn=HamBmn.data;
  Ham.Rmn=HamRmn.data;
  Ham.Zmn=HamZmn.data;
  
  %toc
else
  vthet3D  = theta3D;    %Because Nphi=Nzeta;Nvthet=Ntheta;
  phi3D    = zeta3D;     %Because Nphi=Nzeta;Nvthet=Ntheta;
  vthet3D0c= theta3D0c;   %Because Nphi=Nzeta;Nvthet=Ntheta;
  phi3D0c  = zeta3D0c;    %Because Nphi=Nzeta;Nvthet=Ntheta;

  maxvthet=floor(Nvthet/2)-1; %Nyquist max freq.
  maxphi=floor(Nphi/2)-1;
  Ham.m=0;
  Ham.n=0;
  Ham.Bmn=Ham.B00; %Add this here to make sure it comes first in the list
  ind=1;
  %tic
  if not(Geom.StelSym) %sine components exist
    Ham.parity=1; %this is for the 00 component
    for m=0:floor(Nvthet/2) %Nyquist max freq.
      if m==0 || m==Ntheta/2
        nrange=nvecL0c;%=0:floor(Nphi/2);
        c3D=cos(m * vthet3D0c - n3D0c * NPeriods .* phi3D0c);
        s3D=sin(m * vthet3D0c - n3D0c * NPeriods .* phi3D0c);
      else
        nrange=nvecL;% =-floor(Nphi/2)+Nzeta_even:floor(Nphi/2);
        c3D=cos(m * vthet3D - n3D * NPeriods .* phi3D);
        s3D=sin(m * vthet3D - n3D * NPeriods .* phi3D);
      end
      for n=nrange
        nind=find(nrange==n);  
        c=squeeze(c3D(:,:,nind));
        s=squeeze(s3D(:,:,nind));
        %cos
        if m~=0 && n~=0 %m=n=0 case was already added separately above
          iscorner=((m==0 && n==0) || ...
                    (m==Ntheta/2 && n==0) || ...
                    (m==Ntheta/2 && n==Nzeta/2) || ...
                    (m==0 && n==Nzeta/2));
          Bmnc = (2-iscorner)/(Ntheta*Nzeta) *...
                 sum(sum(c.*Ham.B)); %FAST
                                     %Bmnc = 2/(Ntheta*Nzeta) *...
                                     %       sum(sum(cos(m * Ham.vthet - n * NPeriods * Ham.phi).*Ham.B)); %SLOW
          if abs(Bmnc)>min_Bmn
            ind=ind+1;
            Ham.m(ind)=m;
            Ham.n(ind)=n;
            Ham.Bmn(ind)=Bmnc;
            Ham.parity(ind)=1;
          end
        end
        %sin
        if not(m==0 && n==0) && not(m==0 && n==Nzeta/2) ...
                 && not(m==Ntheta/2 && n==0) && not(m==Ntheta/2 && n==Nzeta/2)
          Bmns = 2/(Ntheta*Nzeta) *...
                 sum(sum(s.*Ham.B)); %FAST
                                     %Bmns = 2/(Ntheta*Nzeta) *...
                                     %       sum(sum(sin(m * Ham.vthet - n * NPeriods * Ham.phi).*Ham.B)); %SLOW
          if abs(Bmns)>min_Bmn
            ind=ind+1;
            Ham.m(ind)=m;
            Ham.n(ind)=n;
            Ham.Bmn(ind)=Bmns;
            Ham.parity(ind)=0;
          end
        end
      end
    end
  else %only cosinus components
    for m=0:floor(Ntheta/2) %Nyquist max freq.
      if m==0 || m==Ntheta/2
        nrange=nvecL0c;%=0:floor(Nphi/2);
        c3D=cos(m * vthet3D0c - n3D0c * NPeriods .* phi3D0c);
      else
        nrange=nvecL;% =-floor(Nphi/2)+Nzeta_even:floor(Nphi/2);
        c3D=cos(m * vthet3D - n3D * NPeriods .* phi3D);
      end
      for n=nrange
        nind=find(nrange==n);
        if m~=0 && n~=0 %m=n=0 case was already added separately above
          c=squeeze(c3D(:,:,nind));
          iscorner=((m==0 && n==0) || ...
                    (m==Ntheta/2 && n==0) || ...
                    (m==Ntheta/2 && n==Nzeta/2) || ...
                    (m==0 && n==Nzeta/2));
          Bmnc = (2-iscorner)/(Ntheta*Nzeta) *...
                 sum(sum(c.*Ham.B)); %FAST
                                     %Bmnc = 2/(Ntheta*Nzeta) *...
                                     %       sum(sum(cos(m * Ham.vthet - n * NPeriods * Ham.phi).*Ham.B)); %SLOW
          if abs(Bmnc)>min_Bmn
            ind=ind+1;
            Ham.m(ind)=m;
            Ham.n(ind)=n;
            Ham.Bmn(ind)=Bmnc;
            Ham.parity(ind)=1;
          end
        end
      end              
    end
  end
  %toc
end
Ham.StelSym=Geom.StelSym;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculate Bmn in Pest coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pest.StelSym=Geom.StelSym;
Pest.B00=sum(sum(Pest.B))/(Pest.Nptheta*Pest.Npzeta);

precision=eps*1e3;
believedAccuracyLoss=1;
min_Bmn=max(Geom.Bfilter.min_Bmn*Booz.B00/Pest.B00*believedAccuracyLoss,precision);

PestBmn=mnlist(fftmn(Pest.B),min_Bmn,'relative');
%PestRmn=mnlist(fftmn(Pest.cylR),min_Bmn,'relative');
%PestZmn=mnlist(fftmn(Pest.cylZ),min_Bmn,'relative');
Pest.m=PestBmn.m;
Pest.n=PestBmn.n;
if Geom.StelSym
  if any(PestBmn.cosparity~=1)
    error(['Something has gone wrong. The calculated Pest Bmn is not stellarator ', ...
           'symmetric!'])
  end
else
  Pest.parity=PestBmn.cosparity;
end  
Pest.Bmn=PestBmn.data;
Pest.StelSym=Geom.StelSym;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Test if the Fourier decomposition was done correctly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0 %Not necessary. I already tested it.
  HamB = zeros(Ntheta,Nzeta);
  for i=1:length(Ham.Bmn)
    if Geom.StelSym
      HamB = HamB + Ham.Bmn(i) *...
             cos(Ham.m(i) * Ham.vthet - Ham.n(i) * NPeriods * Ham.phi);
    else
      if Ham.parity(i) %The cosine components of B
        HamB = HamB + Ham.Bmn(i) *...
               cos(Ham.m(i) * Ham.vthet - Ham.n(i) * NPeriods * Ham.phi);
      else  %The sine components of B
        HamB = HamB + Ham.Bmn(i) *...
               sin(Ham.m(i) * Ham.vthet - Ham.n(i) * NPeriods * Ham.phi);
      end
    end
  end
  fig(1);surf(Ham.vthet,Ham.phi,HamB);view(0,90);shading flat
  xlabel('\vartheta');ylabel('\phi')
  fig(2);surf(Ham.vthet,Ham.phi,Ham.B);view(0,90);shading flat
  xlabel('\vartheta');ylabel('\phi')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some test plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0  
  
  if 0
    fig(3);surf(Booz.theta,Booz.zeta,Booz.B);view(0,90);shading flat
    xlabel('\theta');ylabel('\zeta')
    fig(4);surf(Ham.theta,Ham.zeta,Ham.B);view(0,90);shading flat
    xlabel('\theta');ylabel('\zeta')


    fig(7);surf(Booz.vthet,Booz.phi,Booz.B);view(0,90);shading flat
    xlabel('\vartheta');ylabel('\phi')
    fig(8);surf(Ham.vthet,Ham.phi,Ham.B);view(0,90);shading flat
    xlabel('\vartheta');ylabel('\phi')
  end

  if 0
    fig(1);surf(Booz.theta,Booz.zeta,Booz.vthet);view(0,90);shading flat
    xlabel('\theta');ylabel('\zeta')
    fig(2);surf(Ham.theta,Ham.zeta,Ham.vthet);view(0,90);shading flat
    xlabel('\theta');ylabel('\zeta')


    fig(5);surf(Booz.vthet,Booz.phi,Booz.vthet);view(0,90);shading flat
    xlabel('\vartheta');ylabel('\phi')
    fig(6);surf(Ham.vthet,Ham.phi,Ham.vthet);view(0,90);shading flat
    xlabel('\vartheta');ylabel('\phi')
  end

  if 0
    fig(1);surf(Booz.theta,Booz.zeta,Booz.phi);view(0,90);shading flat
    xlabel('\theta');ylabel('\zeta')
    fig(2);surf(Ham.theta,Ham.zeta,Ham.phi);view(0,90);shading flat
    xlabel('\theta');ylabel('\zeta')


    fig(5);surf(Booz.vthet,Booz.phi,Booz.phi);view(0,90);shading flat
    xlabel('\vartheta');ylabel('\phi')
    fig(6);surf(Ham.vthet,Ham.phi,Ham.phi);view(0,90);shading flat
    xlabel('\vartheta');ylabel('\phi')
  end

  if 0
    fig(1);surf(Booz.theta,Booz.zeta,Booz.Dzetaphi);view(0,90);shading flat;colorbar
    xlabel('\theta');ylabel('\zeta')
    fig(2);surf(Ham.theta,Ham.zeta,Ham.Dzetaphi);view(0,90);shading flat;colorbar
    xlabel('\theta');ylabel('\zeta')


    fig(5);surf(Booz.vthet,Booz.phi,Booz.Dzetaphi);view(0,90);shading flat;colorbar
    xlabel('\vartheta');ylabel('\phi')
    fig(6);surf(Ham.vthet,Ham.phi,Ham.Dzetaphi);view(0,90);shading flat;colorbar
    xlabel('\vartheta');ylabel('\phi')
  end

end