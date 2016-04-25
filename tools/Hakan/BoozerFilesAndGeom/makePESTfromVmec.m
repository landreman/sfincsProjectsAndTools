function [Pest,Vmecrh]=makePESTfromVmec(woutin,s_wish,Nu,Nw)
%Both sets of coordinates are left handed.

%woutin is the wout file name or just the netcdf variables from the wout file
%(Too old matlab versions do not have the necessary netcdf routines.)

if not(isstruct(woutin))%if not already loaded, assume woutin is a string with the file name
  wout=struct();
  if isstr(woutin)
    tmp=ncinfo(woutin);
    for vi=1:length(tmp.Variables)
      wout=setfield(wout,tmp.Variables(vi).Name,...
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

Pest.s=Geom.s(rindh);
Pest.rnorm=Geom.rnorm(rindh);
Pest.iota=Geom.iota(rindh);
Pest.Bphi=Geom.Bphi(rindh);
Pest.Btheta=Geom.Btheta(rindh);

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

Vmecrh.Nu=Nu;
Vmecrh.Nw=Nw;

Vmecrh.Rmnlist.m=double(wout.xm);
Vmecrh.Rmnlist.n=signchange*double(wout.xn)/Geom.Nperiods;
Vmecrh.Rmnlist.cosparity=ones(size(Vmecrh.Rmnlist.m));
Vmecrh.Rmnlist.data=(wout.rmnc(:,rindf_minus)+wout.rmnc(:,rindf_plus))/2;
Vmecrh.Rmn=mnmat(Vmecrh.Rmnlist,Vmecrh.Nu,Vmecrh.Nw,'forceSize');
Vmecrh.R=ifftmn(Vmecrh.Rmn,Geom.Nperiods,Vmecrh.Nu,Vmecrh.Nw,'forceSize');

Vmecrh.Zmnlist.m=double(wout.xm);
Vmecrh.Zmnlist.n=signchange*double(wout.xn)/Geom.Nperiods;
Vmecrh.Zmnlist.cosparity = 0*ones(size(Vmecrh.Zmnlist.m));
Vmecrh.Zmnlist.data=(wout.zmns(:,rindf_minus)+wout.zmns(:,rindf_plus))/2;
Vmecrh.Zmn=mnmat(Vmecrh.Zmnlist,Vmecrh.Nu,Vmecrh.Nw,'forceSize');
Vmecrh.Z=ifftmn(Vmecrh.Zmn,Geom.Nperiods,Vmecrh.Nu,Vmecrh.Nw,'forceSize');

Vmecrh.Bmnlist.m=double(wout.xm_nyq);
Vmecrh.Bmnlist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
Vmecrh.Bmnlist.cosparity=ones(size(Vmecrh.Bmnlist.m));
Vmecrh.Bmnlist.data=wout.bmnc(:,skrindh);
Vmecrh.Bmn=mnmat(Vmecrh.Bmnlist,Vmecrh.Nu,Vmecrh.Nw,'forceSize');
Vmecrh.B=ifftmn(Vmecrh.Bmn,Geom.Nperiods,Vmecrh.Nu,Vmecrh.Nw,'forceSize');

Vmecrh.Jmnlist.m=double(wout.xm_nyq);
Vmecrh.Jmnlist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
Vmecrh.Jmnlist.cosparity=ones(size(Vmecrh.Jmnlist.m));
Vmecrh.Jmnlist.data=wout.gmnc(:,skrindh) * signchange;
Vmecrh.Jmn=mnmat(Vmecrh.Jmnlist,Vmecrh.Nu,Vmecrh.Nw,'forceSize');
Vmecrh.J=ifftmn(Vmecrh.Jmn,Geom.Nperiods,Vmecrh.Nu,Vmecrh.Nw,'forceSize');

Vmecrh.B_ulist.m=double(wout.xm_nyq);
Vmecrh.B_ulist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
Vmecrh.B_ulist.cosparity=ones(size(Vmecrh.B_ulist.m));
Vmecrh.B_ulist.data=wout.bsubumnc(:,skrindh);
Vmecrh.B_umn=mnmat(Vmecrh.B_ulist,Vmecrh.Nu,Vmecrh.Nw,'forceSize');
Vmecrh.B_umntilde=remove00(Vmecrh.B_umn);

Vmecrh.B_wlist.m=double(wout.xm_nyq);
Vmecrh.B_wlist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
Vmecrh.B_wlist.cosparity=ones(size(Vmecrh.B_wlist.m)); 
Vmecrh.B_wlist.data=wout.bsubvmnc(:,skrindh) * signchange;
Vmecrh.B_wmn=mnmat(Vmecrh.B_wlist,Vmecrh.Nu,Vmecrh.Nw,'forceSize');
Vmecrh.B_wmntilde=remove00(Vmecrh.B_wmn);

Vmecrh.llist.m=double(wout.xm);
Vmecrh.llist.n=signchange*double(wout.xn)/Geom.Nperiods;
Vmecrh.llist.cosparity = 0*ones(size(Vmecrh.llist.m));
Vmecrh.llist.data=wout.lmns(:,skrindh);
Vmecrh.lmn=mnmat(Vmecrh.llist,Vmecrh.Nu,Vmecrh.Nw,'forceSize');
Vmecrh.l=ifftmn(Vmecrh.lmn,Geom.Nperiods,Vmecrh.Nu,Vmecrh.Nw,'forceSize');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Pest coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vmecrh.Dw=2*pi/Vmecrh.Nw/Geom.Nperiods;
Vmecrh.Du=2*pi/Vmecrh.Nu;
Vmecrh.uvec=(0:Vmecrh.Nu-1)*Vmecrh.Du;
Vmecrh.wvec=(0:Vmecrh.Nw-1)'*Vmecrh.Dw;
[Vmecrh.u,Vmecrh.w] = ndgrid(Vmecrh.uvec,Vmecrh.wvec);

Vmecrh.Dpthetau=ifftmn(Vmecrh.lmn,Geom.Nperiods,Vmecrh.Nu,Vmecrh.Nw,'forceSize');
Vmecrh.ptheta=Vmecrh.Dpthetau+Vmecrh.u; %This is the "transformation"
Vmecrh.pzeta=Vmecrh.w;

Pest.Nptheta=Vmecrh.Nu;
Pest.Npzeta=Vmecrh.Nw;
Pest.Dptheta=Vmecrh.Du;
Pest.Dpzeta=Vmecrh.Dw;
Pest.pzeta=Vmecrh.w;
Pest.ptheta=Vmecrh.u;

Pest.Dpthetau=griddatacyclic(Vmecrh.ptheta,Vmecrh.pzeta,Vmecrh.Dpthetau,Geom.Nperiods);
Pest.vmecu=Pest.ptheta-Pest.Dpthetau; %This is the "transformation"
Pest.vmecw=Pest.pzeta;

Pest.B=interp2_cyclic(Vmecrh.u,Vmecrh.w,Vmecrh.B,Pest.vmecu,Pest.vmecw,Geom.Nperiods);
Pest.R=interp2_cyclic(Vmecrh.u,Vmecrh.w,Vmecrh.R,Pest.vmecu,Pest.vmecw,Geom.Nperiods);
Pest.Z=interp2_cyclic(Vmecrh.u,Vmecrh.w,Vmecrh.Z,Pest.vmecu,Pest.vmecw,Geom.Nperiods);

Pest.mnmat.B=fftmn(Pest.B);
Pest.mnmat.R=fftmn(Pest.R);
Pest.mnmat.Z=fftmn(Pest.Z);

Pest.R00=mean(mean(Pest.R));
Pest.B00=mean(mean(Pest.B));



  
if 0
  fig(1)
  %surf(Vmecrh.u,Vmecrh.w,Vmecrh.Dzetaw);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  %surf(Vmecrh.u,Vmecrh.w,Vmecrh.zeta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  surf(Vmecrh.u,Vmecrh.w,Vmecrh.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])

  fig(2)
  %surf(Vmecrh.u,Vmecrh.w,Vmecrh.Dthetau);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  surf(Vmecrh.u,Vmecrh.w,Vmecrh.ptheta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  %surf(Vmecrh.u,Vmecrh.w,Vmecrh.u);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])

  %fig(5)
  %surf(Booz_u,Booz_w,Pest.Dzetaw);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  %surf(Booz_u,Booz_w,Pest.zeta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  %surf(Booz_u,Booz_w,Pest.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])


  fig(3)
  surf(Vmecrh.ptheta,Vmecrh.pzeta,Vmecrh.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  fig(7)
  surf(Pest.ptheta,Pest.pzeta,Pest.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])

end