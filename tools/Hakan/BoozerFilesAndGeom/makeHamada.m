function [Ham,Booz,Cyl]=makeHamada(Geom,rind,Ntheta,Nzeta)

useFFT=1;

if isempty(rind)
  error('rind is empty!')
end

NPeriods=Geom.Nperiods;
I=Geom.Btheta(rind);
G=Geom.Bphi(rind);
iota=Geom.iota(rind);
Bmn=Geom.Bmn{rind}; %NOTE: Not normalised to B00
NHarmonics=Geom.nmodes(rind);

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
nvecL=-floor(Nzeta/2):(floor(Nzeta/2)-1);
nvecL0=1:floor(Nzeta/2)-1;
[theta3D, zeta3D, n3D] = ndgrid(thetavec,zetavec,nvecL);
[theta3D0, zeta3D0, n3D0] = ndgrid(thetavec,zetavec,nvecL0);

if Geom.StelSym
  parity=ones(size(Bmn));
else
  parity=Geom.parity{rind};
end

if useFFT 
  %tic
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
  tmp=mnmatrix(Blist);
  [B,R,Z,Dzetacylphi]=ifftmn([Blist,Rlist,Zlist,Dzetacylphilist],NPeriods,Ntheta,Nzeta);
  [dBdtheta,dBdzeta,dRdtheta,dRdzeta,dZdtheta,dZdzeta,...
   dDzetacylphi_dtheta,dDzetacylphi_dzeta]=...
      ifftmn([mngrad(Blist,NPeriods),mngrad(Rlist,NPeriods),...
              mngrad(Zlist,NPeriods),mngrad(Dzetacylphilist,NPeriods)],NPeriods,Ntheta,Nzeta);
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
      
      Z    = Z + Geom.Z{rind}(i) * s;
      dZdtheta = dZdtheta + Geom.Z{rind}(i) * m * c;
      dZdzeta  = dZdzeta  - Geom.Z{rind}(i) * n * NPeriods * c;
      
      Dzetacylphi = Dzetacylphi + 2*pi/NPeriods*Geom.Dphi{rind}(i) * s;
      dDzetacylphi_dtheta = dDzetacylphi_dtheta + 2*pi/NPeriods*Geom.Dphi{rind}(i) * m * c;
      dDzetacylphi_dzeta  = dDzetacylphi_dzeta - 2*pi/NPeriods*Geom.Dphi{rind}(i) * n * NPeriods * c;
      
    else  %The sine components of B
      B = B + Bmn(i) * s;
      dBdtheta = dBdtheta + Bmn(i) * m * c;
      dBdzeta  = dBdzeta  - Bmn(i) * n * NPeriods * c; 
      
      R    = R + Geom.R{rind}(i) * s;
      dRdtheta = dRdtheta + Geom.R{rind}(i) * m * c;
      dRdzeta  = dRdzeta  - Geom.R{rind}(i) * n * NPeriods * c; 
      
      Z    = Z + Geom.Z{rind}(i) * c;
      dZdtheta = dZdtheta - Geom.Z{rind}(i) * m * s;
      dZdzeta  = dZdzeta  + Geom.Z{rind}(i) * n * NPeriods * s;
      
      Dzetacylphi = Dzetacylphi + 2*pi/NPeriods*Geom.Dphi{rind}(i) * c;
      dDzetacylphi_dtheta = dDzetacylphi_dtheta - 2*pi/NPeriods*Geom.Dphi{rind}(i) * m * s;
      dDzetacylphi_dzeta  = dDzetacylphi_dzeta + 2*pi/NPeriods*Geom.Dphi{rind}(i) * n * NPeriods * s;
      
    end
  end
  %toc
end

%Double-check the FFT results against the slow method
%fig(1+useFFT*4)
%surf(B);shading flat;view(0,90)
%fig(2+useFFT*4)
%surf(R);shading flat;view(0,90)
%fig(3+useFFT*4)
%surf(Z);shading flat;view(0,90)
%fig(4+useFFT*4)
%surf(Dzetacylphi);shading flat;view(0,90)
%error('stop stop!')

cylphi=zeta-Dzetacylphi;     %cylphi is minus the geometrical toroidal angle.
geomang=-cylphi;             %(R,Z,cylphi) and (R,geomang,Z) are right handed systems. 
dgeomangdtheta=dDzetacylphi_dtheta;
dgeomangdzeta =dDzetacylphi_dzeta - 1;
X=R.*cos(geomang);
Y=R.*sin(geomang);
dXdtheta=dRdtheta.*cos(geomang)-R.*dgeomangdtheta.*sin(geomang);
dXdzeta =dRdzeta .*cos(geomang)-R.*dgeomangdzeta .*sin(geomang);
dYdtheta=dRdtheta.*sin(geomang)+R.*dgeomangdtheta.*cos(geomang);
dYdzeta =dRdzeta .*sin(geomang)+R.*dgeomangdzeta .*cos(geomang);
g_thetatheta=dXdtheta.^2+dYdtheta.^2+dZdtheta.^2;
g_zetazeta  =dXdzeta.^2 +dYdzeta.^2 +dZdzeta.^2;
g_thetazeta =dXdtheta.*dXdzeta+dYdtheta.*dYdzeta+dZdtheta.*dZdzeta;

gradpsi.X=B.^2./(G+iota*I).*(dYdtheta.*dZdzeta-dZdtheta.*dYdzeta);
gradpsi.Y=B.^2./(G+iota*I).*(dZdtheta.*dXdzeta-dXdtheta.*dZdzeta);
gradpsi.Z=B.^2./(G+iota*I).*(dXdtheta.*dYdzeta-dYdtheta.*dZdzeta);
gpsipsi=gradpsi.X.^2+gradpsi.Y.^2+gradpsi.Z.^2;

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

Booz.Dzetacylphi=Dzetacylphi;     %for later use
Booz.dBdtheta=dBdtheta;
Booz.dBdzeta=dBdzeta;
Booz.dRdtheta=dRdtheta;
Booz.dRdzeta=dRdzeta;
Booz.dZdtheta=dZdtheta;
Booz.dZdzeta=dZdzeta;

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
%dudtheta = zeros(Ntheta,Nzeta);
%dudzeta = zeros(Ntheta,Nzeta);
h=1./(B.^2);
VPrimeHat=sum(sum(h))*4*pi^2/(Nzeta*Ntheta); %Note: VPrime=VPrimeHat*(G+iota*I)
FSAB2=4*pi^2/VPrimeHat;
%h00=VPrimeHat/(4*pi^2);
if useFFT
  %tic
  hmn=fftmn(h);
  umn.c=iota*(G*hmn.m + I*hmn.n * NPeriods)./(hmn.n * NPeriods - iota*hmn.m) .* hmn.c;
  umn.c(hmn.m0ind,hmn.n0ind)=0;
  umn.s=iota*(G*hmn.m + I*hmn.n * NPeriods)./(hmn.n * NPeriods - iota*hmn.m) .* hmn.s;
  umn.s(hmn.m0ind,hmn.n0ind)=NaN;
  Dzetaphi_mn.s  =  FSAB2/(G+iota*I)./hmn.m.*(umn.c/iota-I*hmn.c);     %for m~=0
  Dzetaphi_mn.c  = -FSAB2/(G+iota*I)./hmn.m.*(umn.s/iota-I*hmn.s);    %for m~=0
  Dzetaphi2_mn.s =  FSAB2/(G+iota*I)./hmn.n/NPeriods.*(umn.c+G*hmn.c); %for n~=0
  Dzetaphi2_mn.c = -FSAB2/(G+iota*I)./hmn.n/NPeriods.*(umn.s+G*hmn.s); %for n~=0
  Dzetaphi_mn.s(hmn.m0ind,:)=Dzetaphi2_mn.s(hmn.m0ind,:);
  Dzetaphi_mn.c(hmn.m0ind,:)=Dzetaphi2_mn.c(hmn.m0ind,:);
  Dzetaphi_mn.s(hmn.m0ind,hmn.n0ind)=NaN;
  Dzetaphi_mn.c(hmn.m0ind,hmn.n0ind)=0;
  Dzetaphi_mn.c(hmn.m0ind,end)=0; %info was lost here
  Dzetaphi_mn.c(end,end)=0;       %info was lost here
  Dzetaphi_mn.c(end,hmn.n0ind)=0; %info was lost here
  [u,Dzetaphi]=ifftmn([umn,Dzetaphi_mn]);
  %toc
else %not(useFFT)
  %tic
    if not(Geom.StelSym) %sine components exist
      for m=0:floor(Ntheta/2)-1 %Nyquist max freq.
        if m==0
          nrange=nvecL0;%=1:floor(Nzeta/2)-1;
          c3D=cos(-n3D0 * NPeriods .* zeta3D0);
          s3D=sin(-n3D0 * NPeriods .* zeta3D0);
        else
          nrange=nvecL;% =-floor(Nzeta/2):(floor(Nzeta/2)-1);
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
          hmnc = 2/(Ntheta*Nzeta) *...
                 sum(sum(c.*h));
          umnc = ...
              iota*(G*m + I*n * NPeriods)/(n * NPeriods - iota*m) * hmnc;
          %There are two ways of calculating Dzetaphi. We need both.
          Dzetaphi1_mns = ...
              FSAB2/(G+iota*I)/m*(umnc/iota-I*hmnc);     %for m~=0
          Dzetaphi2_mns = ...
              FSAB2/(G+iota*I)/n/NPeriods*(umnc+G*hmnc); %for n~=0
          
          %sin
          hmns = 2/(Ntheta*Nzeta) *...
                 sum(sum(s.*h));
          umns = ...
              iota*(G*m + I*n * NPeriods)/(n * NPeriods - iota*m) * hmns;
          Dzetaphi1_mnc = ...  
              -FSAB2/(G+iota*I)/m*(umns/iota-I*hmns);     %for m~=0
          Dzetaphi2_mnc = ...
              -FSAB2/(G+iota*I)/n/NPeriods*(umns+G*hmns); %for n~=0

          %assemble
          u = u + umnc * c + umns * s;
          if m>0
            Dzetaphi = Dzetaphi + Dzetaphi1_mnc * c + Dzetaphi1_mns * s;    
          else
            Dzetaphi = Dzetaphi + Dzetaphi2_mnc * c + Dzetaphi2_mns * s; 
          end
          %dudtheta = dudtheta ...
          %    + umns * m * c - umnc * m * s;
          %dudzeta = dudzeta ...
          %    + umnc * n * NPeriods * s - umns * n * NPeriods * c; 
        end
      end
    else %only cosinus components
      for m=0:floor(Ntheta/2)-1 %Nyquist max freq.
        if m==0
          nrange=nvecL0;%=1:floor(Nzeta/2)-1;
          c3D=cos(-n3D0 * NPeriods .* zeta3D0);
          s3D=sin(-n3D0 * NPeriods .* zeta3D0);
        else
          nrange=nvecL;% =-floor(Nzeta/2):(floor(Nzeta/2)-1);
          c3D=cos(m * theta3D - n3D * NPeriods .* zeta3D);
          s3D=sin(m * theta3D - n3D * NPeriods .* zeta3D);
        end
        for n=nrange
          nind=find(nrange==n);  
          c=squeeze(c3D(:,:,nind)); %FAST
          s=squeeze(s3D(:,:,nind));
          %c=cos(m * theta - n * NPeriods * zeta); %SLOW
          %s=sin(m * theta - n * NPeriods * zeta);

          hmnc = 2/(Ntheta*Nzeta) *...
                 sum(sum(c.*h));
          umnc = ...
              iota*(G*m + I*n * NPeriods)/(n * NPeriods - iota*m) * hmnc;
          Dzetaphi1_mns = ...
              FSAB2/(G+iota*I)/m*(umnc/iota-I*hmnc);     %for m~=0
          Dzetaphi2_mns = ...
              FSAB2/(G+iota*I)/n/NPeriods*(umnc+G*hmnc); %for n~=0

          u = u + umnc * c;
          if m>0
            Dzetaphi = Dzetaphi + Dzetaphi1_mns * s;    
          else
            Dzetaphi = Dzetaphi + Dzetaphi2_mns * s; 
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


dDzetaphi_dtheta= FSAB2/(G+iota*I)/iota*(u-iota*I*(h-1/FSAB2));
dDzetaphi_dzeta =-FSAB2/(G+iota*I) *    (u +    G*(h-1/FSAB2));

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


if 0 %Double-check the FFT results against the slow method
  fig(1+useFFT*4)
  surf(g_phiphi);shading flat;view(0,90)
  fig(2+useFFT*4)
  surf(dXdphi);shading flat;view(0,90)
  fig(3+useFFT*4)
  surf(dYdphi);shading flat;view(0,90)
  fig(4+useFFT*4)
  surf(dZdphi);shading flat;view(0,90)
end

Booz.Nperiods=Geom.Nperiods;
Booz.Nzeta=Nzeta;
Booz.Ntheta=Ntheta;
Booz.Dzeta=Dzeta;
Booz.Dtheta=Dtheta;
Booz.zeta=zeta;
Booz.theta=theta;

Booz.Dzetaphi=Dzetaphi;
Booz.phi=zeta-Booz.Dzetaphi;
Booz.vthet=theta-iota*Booz.Dzetaphi;
Booz.cylR=R;
Booz.cylZ=Z;
Booz.cylphi=cylphi;
if cylPossible
  Booz.cylth=cylth;
  Booz.cylr=cylr;
  Booz.cylR00=cylR00;
end
Booz.iota=Geom.iota(rind);
Booz.G=Geom.Bphi(rind);
Booz.I=Geom.Btheta(rind);
Booz.B=B;
Booz.u=u;
Booz.h=h;
Booz.FSAB2=FSAB2;
Booz.FSAu2B2=sum(sum(u.^2.*B.^2.*h))/sum(sum(h));
Booz.Jacob=Booz.h*(G+iota*I);
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
if not(Geom.StelSym) %sine components exist
  Booz.parity=Geom.parity{rind};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate to a phi vthet grid 
% Here, phi is the Hamada tor. coord. and vthet (\vartheta) is the Hamada poloidal coordinate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

vthet3D = theta3D;    %Because Nphi=Nzeta;Nvthet=Ntheta;
phi3D   = zeta3D;     %Because Nphi=Nzeta;Nvthet=Ntheta;
vthet3D0= theta3D0;   %Because Nphi=Nzeta;Nvthet=Ntheta;
phi3D0  = zeta3D0;    %Because Nphi=Nzeta;Nvthet=Ntheta;


BoozphiBig=[Booz.phi-2*pi/NPeriods, Booz.phi, Booz.phi+2*pi/NPeriods;
            Booz.phi-2*pi/NPeriods, Booz.phi, Booz.phi+2*pi/NPeriods;
            Booz.phi-2*pi/NPeriods, Booz.phi, Booz.phi+2*pi/NPeriods];
            
BoozvthetBig=[Booz.vthet-2*pi,Booz.vthet-2*pi,Booz.vthet-2*pi;
              Booz.vthet,     Booz.vthet,     Booz.vthet;
              Booz.vthet+2*pi,Booz.vthet+2*pi,Booz.vthet+2*pi];
%BoozBBig=[Booz.B,Booz.B,Booz.B;
%          Booz.B,Booz.B,Booz.B;
%          Booz.B,Booz.B,Booz.B];
Booz.DzetaphiBig=[Booz.Dzetaphi,Booz.Dzetaphi,Booz.Dzetaphi;
                 Booz.Dzetaphi,Booz.Dzetaphi,Booz.Dzetaphi;
                 Booz.Dzetaphi,Booz.Dzetaphi,Booz.Dzetaphi];

%The following is the most time consuming step
%tic
HamDzetaphi=griddata(BoozphiBig,BoozvthetBig,Booz.DzetaphiBig,phi,vthet);
%toc 

Ham.zeta =Ham.phi+HamDzetaphi;
Ham.theta=Ham.vthet+iota*HamDzetaphi;

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

Ham.iota=Booz.iota;
Ham.G=Booz.G;
Ham.I=Booz.I;

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

% The quantity h=1/B^2:
if 0 %interpolate
  Boozhbig=[Booz.h,Booz.h(:,1)];
  Boozhbig=[Boozhbig;Boozhbig(1,:)];
  Ham.h=interp2(Boozzetabig,Boozthetabig,Boozhbig,...
                mod(Ham.zeta,2*pi/NPeriods),mod(Ham.theta,2*pi)); 
else
  Ham.h=1./Ham.B.^2;
end

%Ham.FSAB2=(Nvthet*Nphi)/sum(sum(Ham.h)); %Wrong, only holds for Boozer
Ham.FSAB2=Booz.FSAB2;
Ham.Jacob=1/Ham.FSAB2;
Ham.FSAg_phiphi=Booz.FSAg_phiphi;
Ham.FSAgpsipsi =Booz.FSAgpsipsi;
Ham.FSAu2B2=Booz.FSAu2B2;

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

  % Hamada toroidal angle
  Boozphibig=[Booz.phi,Booz.phi(:,1)+zetaDirection*2*pi/NPeriods];
  Boozphibig=[Boozphibig;Boozphibig(1,:)];
  Cyl.phi=interp2(Boozzetabig,Boozthetabig,Boozphibig,...
                  mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi))+...
          Cyl.zeta-mod(Cyl.zeta,2*pi/NPeriods); 
  % Hamada poloidal angle
  Boozvthetbig=[Booz.vthet;Booz.vthet(1,:)+thetaDirection*2*pi];
  Boozvthetbig=[Boozvthetbig,Boozvthetbig(:,1)];
  Cyl.vthet=interp2(Boozzetabig,Boozthetabig,Boozvthetbig,...
                    mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi))+...
            Cyl.theta-mod(Cyl.theta,2*pi); 

  Cyl.cylR=interp2(Boozzetabig,Boozthetabig,BoozRbig,...
                   mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi)); 
  Cyl.cylZ=interp2(Boozzetabig,Boozthetabig,BoozZbig,...
                   mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi)); 
  Cyl.cylr=sqrt(Cyl.cylZ.^2+(Cyl.cylR-cylR00).^2);
  Cyl.cylR00=cylR00;
  
  Cyl.B=interp2(Boozzetabig,Boozthetabig,BoozBbig,...
                mod(Cyl.zeta,2*pi/NPeriods),mod(Cyl.theta,2*pi)); 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculate Bmn in Hamada coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ham.B00=sum(sum(Ham.B))/(Nvthet*Nphi);

precision=eps*1e3;
believedAccuracyLoss=1;
min_Bmn=max(Geom.Bfilter.min_Bmn*Booz.B00/Ham.B00*believedAccuracyLoss,precision);

if useFFT
  %tic
  HamBmn=mnlist(fftmn(Ham.B),min_Bmn,'relative');
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
  %toc
else
  maxvthet=floor(Nvthet/2)-1; %Nyquist max freq.
  maxphi=floor(Nphi/2)-1;
  Ham.m=0;
  Ham.n=0;
  Ham.Bmn=Ham.B00;
  ind=1;
  %tic
  if not(Geom.StelSym) %sine components exist
    Ham.parity=1; %this is for the 00 component
    for m=0:floor(Nvthet/2)-1 %Nyquist max freq.
      if m==0
        nrange=nvecL0;%=1:floor(Nphi/2)-1;
        c3D=cos(-n3D0 * NPeriods .* phi3D0);
        s3D=sin(-n3D0 * NPeriods .* phi3D0);
      else
        nrange=nvecL;% =-floor(Nphi/2):(floor(Nphi/2)-1);
        c3D=cos(m * vthet3D - n3D * NPeriods .* phi3D);
        s3D=sin(m * vthet3D - n3D * NPeriods .* phi3D);
      end
      for n=nrange
        nind=find(nrange==n);  
        c=squeeze(c3D(:,:,nind));
        s=squeeze(s3D(:,:,nind));
        %cos
        Bmnc = 2/(Ntheta*Nzeta) *...
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
        %sin
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
  else %only cosinus components
    for m=0:floor(Ntheta/2)-1 %Nyquist max freq.
      if m==0
        nrange=nvecL0;%=1:floor(Nphi/2)-1;
        c3D=cos(-n3D0 * NPeriods .* phi3D0);
      else
        nrange=nvecL;% =-floor(Nphi/2):(floor(Nphi/2)-1);
        c3D=cos(m * vthet3D - n3D * NPeriods .* phi3D);
      end
      for n=nrange
        nind=find(nrange==n);
        c=squeeze(c3D(:,:,nind));
        Bmnc = 2/(Ntheta*Nzeta) *...
               sum(sum(c.*Ham.B)); %FAST
                                   %Bmnc = 2/(Ntheta*Nzeta) *...
                                   %       sum(sum(cos(m * Ham.vthet - n * NPeriods * Ham.phi).*Ham.B)); %SLOW
        if abs(Bmnc)>min_Bmn
          ind=ind+1;
          Ham.m(ind)=m;
          Ham.n(ind)=n;
          Ham.Bmn(ind)=Bmnc;
        end
      end              
    end
  end
  %toc
end
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