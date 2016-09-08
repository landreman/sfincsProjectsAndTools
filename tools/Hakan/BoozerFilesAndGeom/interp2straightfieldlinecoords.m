function [Booz,New]=interp2straightfieldlinecoords(Booz,New,Dzeta_Newtor,polname,torname)

%This is used to interpolate Boozer coordinates to other field aligned coordinates,
%such as Hamada or Pest coordinates.
%The names vthet and phi are the ones I usually use to refer to Hamada, but when the
%input Dzeta_Newtor is something else than Dzetadphi other coordinates come out.
%the names in the struct will be given by polname,torname

Booz=setfield(Booz,['Dzeta',torname],Dzeta_Newtor);
Booz=setfield(Booz,torname,Booz.zeta-Dzeta_Newtor);
Booz=setfield(Booz,polname,Booz.theta-Booz.iota*Dzeta_Newtor);
%Locally here, use the names phi and vthet.
Boozphi=Booz.zeta-Dzeta_Newtor;
Boozvthet=Booz.theta-Booz.iota*Dzeta_Newtor;
BoozDzetaphi=Dzeta_Newtor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate to a phi vthet grid 
% Here, phi is the Hamada tor. coord. and vthet (\vartheta) is the Hamada poloidal coordinate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Choose the same resolution. Do not change this. It means 3D tensors can be used again 

Nphi=Booz.Nzeta;
Nvthet=Booz.Ntheta;
NPeriods=Booz.Nperiods;
iota=Booz.iota;
New.Nperiods=NPeriods;
New=setfield(New,['N',torname],Nphi);
New=setfield(New,['N',polname],Nvthet);
%New.Nphi=Nphi;
%New.Nvthet=Nvthet;

%New discretisation

Dphi=2*pi/Nphi/NPeriods;
Dvthet=2*pi/Nvthet;
phivec=(0:Nphi-1)*Dphi;
vthetvec=(0:Nvthet-1)'*Dvthet;
[vthet, phi] = ndgrid(vthetvec,phivec);

New=setfield(New,['D',torname],Dphi);
New=setfield(New,['D',polname],Dvthet);
New=setfield(New,torname,phi);
New=setfield(New,polname,vthet);


Nleftadd=ceil(max(Booz.phi(:,1))/(2*pi/NPeriods));
Nrightadd=ceil(1-min(Booz.phi(:,end))/(2*pi/NPeriods));
BoozphiBig      =NaN*zeros(Booz.Ntheta*3,Nphi*(Nleftadd+1+Nrightadd));
BoozvthetBig    =NaN*zeros(Booz.Ntheta*3,Nphi*(Nleftadd+1+Nrightadd));
BoozDzetaphiBig=NaN*zeros(Booz.Ntheta*3,Nphi*(Nleftadd+1+Nrightadd));

for addi=1:Nleftadd+1+Nrightadd
  addci=addi-(Nleftadd+1);
  toaddphi=2*pi/NPeriods*addci;
  BoozphiBig(1:Booz.Ntheta*3,(addi-1)*Nphi+1:addi*Nphi)=...
      [Boozphi+toaddphi;
       Boozphi+toaddphi;
       Boozphi+toaddphi];
  BoozvthetBig(1:Booz.Ntheta*3,(addi-1)*Nphi+1:addi*Nphi)=...
      [Boozvthet-2*pi;
       Boozvthet;
       Boozvthet+2*pi];
  BoozDzetaphiBig(1:Booz.Ntheta*3,(addi-1)*Nphi+1:addi*Nphi)=...
      [BoozDzetaphi;
       BoozDzetaphi;
       BoozDzetaphi];
end

%The following is the most time consuming step
%tic
Dzetaphi=griddata(BoozphiBig,BoozvthetBig,BoozDzetaphiBig,phi,vthet);
%toc 
New.zeta =phi+Dzetaphi;
New.theta=vthet+iota*Dzetaphi;
New=setfield(New,['Dzeta',torname],Dzetaphi);

zetaDirection=sign(Booz.zeta(1,end)-Booz.zeta(1,1));
thetaDirection=sign(Booz.theta(end,1)-Booz.theta(1,1));

Boozzetabig=[Booz.zeta,Booz.zeta(:,1)+zetaDirection*2*pi/NPeriods];
Boozzetabig=[Boozzetabig;Boozzetabig(1,:)];
Boozthetabig=[Booz.theta;Booz.theta(1,:)+thetaDirection*2*pi];
Boozthetabig=[Boozthetabig,Boozthetabig(:,1)];
BoozBbig=[Booz.B,Booz.B(:,1)];
BoozBbig=[BoozBbig;BoozBbig(1,:)];

% major radius R
BoozRbig=[Booz.cylR,Booz.cylR(:,1)];
BoozRbig=[BoozRbig;BoozRbig(1,:)];
New.cylR=interp2(Boozzetabig,Boozthetabig,BoozRbig,...
              mod(New.zeta,2*pi/NPeriods),mod(New.theta,2*pi)); 
% vertical coordinate Z
BoozZbig=[Booz.cylZ,Booz.cylZ(:,1)];
BoozZbig=[BoozZbig;BoozZbig(1,:)];
New.cylZ=interp2(Boozzetabig,Boozthetabig,BoozZbig,...
              mod(New.zeta,2*pi/NPeriods),mod(New.theta,2*pi)); 
% Geometrical toroidal angle (cylphi)
Boozcylphibig=[Booz.cylphi,Booz.cylphi(:,1)+zetaDirection*2*pi/NPeriods];
Boozcylphibig=[Boozcylphibig;Boozcylphibig(1,:)];
New.cylphi=interp2(Boozzetabig,Boozthetabig,Boozcylphibig,...
              mod(New.zeta,2*pi/NPeriods),mod(New.theta,2*pi))+...
              New.zeta-mod(New.zeta,2*pi/NPeriods); 
New.R00=mean(mean(New.cylR));

if Booz.cylPossible
  % Geometrical poloidal angle (cylth)
  Boozcylthbig=[Booz.cylth;Booz.cylth(1,:)+thetaDirection*2*pi];
  Boozcylthbig=[Boozcylthbig,Boozcylthbig(:,1)];
  New.cylth=interp2(Boozzetabig,Boozthetabig,Boozcylthbig,...
                    mod(New.zeta,2*pi/NPeriods),mod(New.theta,2*pi))+...
            New.theta-mod(New.theta,2*pi); 

  New.cylr=sqrt(New.cylZ.^2+(New.cylR-Booz.cylR00).^2);

  New.cylR00=Booz.cylR00;
end
New.X=New.cylR.*cos(-New.cylphi);
New.Y=New.cylR.*sin(-New.cylphi);

New.iota=Booz.iota;
New.G=Booz.G;
New.I=Booz.I;

%We have two options on how to interpolate. 
%We use interp2, because it is much faster than griddata
%Option 1:
%New.B=griddata(BoozphiBig,BoozvthetBig,BoozBBig,phi,vthet); %time=1 second
%Option 2:
New.B=interp2(Boozzetabig,Boozthetabig,BoozBbig,...
              mod(New.zeta,2*pi/NPeriods),mod(New.theta,2*pi)); %time=3 milliseconds

% Interpolate the rest of the interesting stored data:
% The parallel curent u:
Boozubig=[Booz.u,Booz.u(:,1)];
Boozubig=[Boozubig;Boozubig(1,:)];
New.u=interp2(Boozzetabig,Boozthetabig,Boozubig,...
              mod(New.zeta,2*pi/NPeriods),mod(New.theta,2*pi)); 

% The quantity h=1/B^2:
if 0 %interpolate
  Boozhbig=[Booz.h,Booz.h(:,1)];
  Boozhbig=[Boozhbig;Boozhbig(1,:)];
  New.h=interp2(Boozzetabig,Boozthetabig,Boozhbig,...
                mod(New.zeta,2*pi/NPeriods),mod(New.theta,2*pi)); 
else
  New.h=1./New.B.^2;
end

New.FSAB2=Booz.FSAB2;
New.FSAg_phiphi=Booz.FSAg_phiphi;
New.FSAgpsipsi =Booz.FSAgpsipsi;
New.FSAu2B2=Booz.FSAu2B2;
 