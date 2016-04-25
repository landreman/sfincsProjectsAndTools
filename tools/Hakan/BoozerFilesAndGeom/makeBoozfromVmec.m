function Booz=makeBoozfromVmec(woutin,s_wish,Nu,Nw,min_Bmn)

%This routine has been tested and compared to jmc and/or booz2xform to find that B00
%and R00 coincides.

%woutin is the wout file name or just the netcdf variables from the wout file
%(Too old matlab versions do not have the necessary netcdf routines.)

if not(isstruct(woutin))%if not already loaded, assume woutin is a string with the file name
  wout=struct();
  if isstr(woutin)
    tmp=ncinfo(woutin);
    for vi=1:length(tmp.Variables)
      data=setfield(wout,tmp.Variables(vi).Name,...
                         ncread(woutin,tmp.Variables(vi).Name));
    end
  else
    error('Unknown input. wout must be a file name or the netcdf variables from the wout file!')
  end
end



% let w = -v!

signchange=double(wout.signgs); % is -1, because vmec is left handed

Geom.headertext.input_extension = ...
    wout.input_extension; %!<  suffix of the vmec-input file: input.<input_extension>

Geom.m0b      = double(wout.mpol);       %!< number of poloidal fourier modes
Geom.n0b      = double(wout.ntor);       %!< upper bound toroidal fourier modes: -ntor <= n <= ntor
Geom.nsurf    = double(wout.ns);         %!< number of radial surfaces
Geom.Nperiods = double(wout.nfp);        %!< number of field periods
if not(isfield(wout,'iasym'))
  wout.iasym=isfield(wout,'bmns');
end
Geom.StelSym  = not(wout.iasym); %!<  defines stellarator symmetry for iasym=0, otherwise =
Geom.torfluxtot = ...
   wout.phi(wout.ns)*signchange;%!<  total toroidal flux within the boundary (s=1)
Geom.psi_a=Geom.torfluxtot/2/pi;

Geom.minorradiusVMEC= wout.Aminor_p;  %!<  minor plasma radius
Geom.majorradiusLastbcR00 = NaN; %Calculate this below (not necessary)
Geom.minorradiusW7AS= NaN; %Calculate this below (not necessary)
Geom.majorradiusVMEC= wout.Rmajor_p;  %!<  major plasma radius

Geom.fullgrid.s     = wout.phi/wout.phi(wout.ns);     %full grid
Geom.fullgrid.rnorm = sqrt(Geom.fullgrid.s);     %full grid

skip=1; %this is how many elements are skipped at low radii when going to half grid

Geom.s=(Geom.fullgrid.s(skip:end-1)+Geom.fullgrid.s(skip+1:end))/2; %half grid
Geom.rnorm=sqrt(Geom.s); %half grid

%Geom.dVdsoverNper=dVdsoverNper*signchange;
Geom.Bphi  = wout.bvco*signchange;%direction switch sign
Geom.Btheta= wout.buco;%*signchange;

Geom.Bphi=Geom.Bphi(skip+1:end);
Geom.Btheta=Geom.Btheta(skip+1:end);

Geom.fullgrid.iota = wout.iotaf*signchange; 
Geom.iota = wout.iotas*signchange; %half mesh
Geom.iota = Geom.iota(skip+1:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specific for this radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dummy,rindh]=min(abs(Geom.s-s_wish));
s=Geom.s(rindh);

Booz.s=Geom.s(rindh);
Booz.rnorm=Geom.rnorm(rindh);
Booz.iota=Geom.iota(rindh);
Booz.Bphi=Geom.Bphi(rindh);
Booz.Btheta=Geom.Btheta(rindh);

G=Geom.Bphi(rindh);   %=wout.bsubvmnc(n0ind,m0ind,rind)
I=Geom.Btheta(rindh); %=wout.bsubumnc(n0ind,m0ind,rind)
iota=Geom.iota(rindh);

skrindh=skip+rindh; %shorthand notation.
rindf_plus =skrindh;
rindf_minus=skrindh-1;

%now use a RH coordinate system (u,w,s) with w=-v
if not(Geom.StelSym)
  error('Non-stelllarator symmmetric case not implemented!')
end

w.Nu=Nu;
w.Nw=Nw;

w.Rmnlist.m=double(wout.xm);
w.Rmnlist.n=signchange*double(wout.xn)/Geom.Nperiods;
w.Rmnlist.cosparity=ones(size(w.Rmnlist.m));
w.Rmnlist.data=(wout.rmnc(:,rindf_minus)+wout.rmnc(:,rindf_plus))/2;
w.Rmn=mnmat(w.Rmnlist,w.Nu,w.Nw,'forceSize');
w.R=ifftmn(w.Rmn,Geom.Nperiods,w.Nu,w.Nw,'forceSize');

w.Zmnlist.m=double(wout.xm);
w.Zmnlist.n=signchange*double(wout.xn)/Geom.Nperiods;
w.Zmnlist.cosparity = 0*ones(size(w.Zmnlist.m));
w.Zmnlist.data=(wout.zmns(:,rindf_minus)+wout.zmns(:,rindf_plus))/2;
w.Zmn=mnmat(w.Zmnlist,w.Nu,w.Nw,'forceSize');
w.Z=ifftmn(w.Zmn,Geom.Nperiods,w.Nu,w.Nw,'forceSize');

w.Bmnlist.m=double(wout.xm_nyq);
w.Bmnlist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
w.Bmnlist.cosparity=ones(size(w.Bmnlist.m));
w.Bmnlist.data=wout.bmnc(:,skrindh);
w.Bmn=mnmat(w.Bmnlist,w.Nu,w.Nw,'forceSize');
w.B=ifftmn(w.Bmn,Geom.Nperiods,w.Nu,w.Nw,'forceSize');

w.Jmnlist.m=double(wout.xm_nyq);
w.Jmnlist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
w.Jmnlist.cosparity=ones(size(w.Jmnlist.m));
w.Jmnlist.data=wout.gmnc(:,skrindh) * signchange;
w.Jmn=mnmat(w.Jmnlist,w.Nu,w.Nw,'forceSize');
w.J=ifftmn(w.Jmn,Geom.Nperiods,w.Nu,w.Nw,'forceSize');

w.B_ulist.m=double(wout.xm_nyq);
w.B_ulist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
w.B_ulist.cosparity=ones(size(w.B_ulist.m));
w.B_ulist.data=wout.bsubumnc(:,skrindh);
w.B_umn=mnmat(w.B_ulist,w.Nu,w.Nw,'forceSize');
w.B_umntilde=remove00(w.B_umn);

w.B_wlist.m=double(wout.xm_nyq);
w.B_wlist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
w.B_wlist.cosparity=ones(size(w.B_wlist.m)); 
w.B_wlist.data=wout.bsubvmnc(:,skrindh) * signchange;
w.B_wmn=mnmat(w.B_wlist,w.Nu,w.Nw,'forceSize');
w.B_wmntilde=remove00(w.B_wmn);

w.llist.m=double(wout.xm);
w.llist.n=signchange*double(wout.xn)/Geom.Nperiods;
w.llist.cosparity = 0*ones(size(w.llist.m));
w.llist.data=wout.lmns(:,skrindh);
w.lmn=mnmat(w.llist,w.Nu,w.Nw,'forceSize');
w.l=ifftmn(w.lmn,Geom.Nperiods,w.Nu,w.Nw,'forceSize');

%%%%% This is the transformation. See th notes by Hirshman
%%%%% 'Transformation from VMEC to Boozer Coordinates', April 1995
alphamn1=invgrad(w.B_umntilde,w.B_wmntilde,Geom.Nperiods,1);
alphamn2=invgrad(w.B_umntilde,w.B_wmntilde,Geom.Nperiods,2);

w.pmn=(alphamn2-I*w.lmn)*(1/(G+iota*I));

w.Dw=2*pi/w.Nw/Geom.Nperiods;
w.Du=2*pi/w.Nu;
w.uvec=(0:w.Nu-1)*w.Du;
w.wvec=(0:w.Nw-1)'*w.Dw;
[w.u,w.w] = ndgrid(w.uvec,w.wvec);

w.Dthetau=ifftmn(w.lmn+iota*w.pmn,Geom.Nperiods,w.Nu,w.Nw,'forceSize');
w.Dzetaw =ifftmn(w.pmn,Geom.Nperiods,w.Nu,w.Nw,'forceSize');
w.theta=w.Dthetau+w.u;
w.zeta=w.Dzetaw+w.w;

Booz.Ntheta=w.Nu;
Booz.Nzeta=w.Nw;
Booz.Dtheta=w.Du;
Booz.Dzeta=w.Dw;
Booz.zeta=w.w;
Booz.theta=w.u;

Booz.Dthetau=griddatacyclic(w.theta,w.zeta,w.Dthetau,Geom.Nperiods);
Booz.Dzetaw=griddatacyclic(w.theta,w.zeta,w.Dzetaw,Geom.Nperiods);
Booz_u=Booz.theta-Booz.Dthetau;
Booz_w=Booz.zeta-Booz.Dzetaw;

Booz.B=interp2_cyclic(w.u,w.w,w.B,Booz_u,Booz_w,Geom.Nperiods);
Booz.R=interp2_cyclic(w.u,w.w,w.R,Booz_u,Booz_w,Geom.Nperiods);
Booz.Z=interp2_cyclic(w.u,w.w,w.Z,Booz_u,Booz_w,Geom.Nperiods);

Booz.mnmat.B=fftmn(Booz.B);
Booz.mnmat.R=fftmn(Booz.R);
Booz.mnmat.Z=fftmn(Booz.Z);

Booz.R00=mean(mean(Booz.R));
Booz.B00=mean(mean(Booz.B));

%[mean(mean(Booz.R)),5.5170]
%[mean(mean(Booz.B)),2.8234]



if 0
fig(1)
%surf(w.u,w.w,w.Dzetaw);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
%surf(w.u,w.w,w.zeta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
surf(w.u,w.w,w.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])

fig(2)
%surf(w.u,w.w,w.Dthetau);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
surf(w.u,w.w,w.theta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
%surf(w.u,w.w,w.u);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])

fig(5)
%surf(Booz_u,Booz_w,Booz.Dzetaw);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
%surf(Booz_u,Booz_w,Booz.zeta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
surf(Booz_u,Booz_w,Booz.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])

fig(6)
%surf(Booz_u,Booz_w,Booz.Dthetau);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
surf(Booz_u,Booz_w,Booz.theta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
%surf(Booz_u,Booz_w,Booz.interpu);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])

fig(3)
surf(w.theta,w.zeta,w.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
fig(7)
surf(Booz.theta,Booz.zeta,Booz.B);shading flat;view(0,90);colorbar;axis( ...
    [-0.5,6.5,-0.2,1.4])

end