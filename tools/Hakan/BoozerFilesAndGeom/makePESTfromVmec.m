function [Pest,Vmec]=makePESTfromVmec(woutin,s_wish,Nu,Nw,handedness)

% woutin is the wout file name or just the netcdf variables from the wout file
% (Too old matlab versions do not have the necessary netcdf routines.)
%
% s_wish is the wanted flux surface
%
% (Nu, Nw) is the (poloidal, toroidal) resolution
%
% handedness = 1  => right-handed output structs
% handedness =-1  => left-handed output structs


if not(isstruct(woutin))%if not already loaded, assume woutin is a string with the file name
  wout=struct();
  if isstr(woutin)
    tmp=ncinfo(woutin);
    for vi=1:length(tmp.Variables)
      wout=setfield(wout,tmp.Variables(vi).Name,...
                         ncread(woutin,tmp.Variables(vi).Name));
    end
  else
    error('Unknown input! woutin must be a file name or the netcdf variables from the wout file!')
  end
end


% if handedness=1 (right handed) let w = -v!
% else if handedness=-1 (left handed) let w = v!


signchange=handedness*double(wout.signgs); % is handedness*(-1), because vmec is left handed

Geom.headertext.input_extension = ...
    wout.input_extension; %!<  suffix of the vmec-input file: input.<input_extension>

Geom.m0b      = double(wout.mpol);       %!< number of poloidal fourier modes
Geom.n0b      = double(wout.ntor);       %!< upper bound toroidal fourier modes: -ntor <= n <= ntor
Geom.nsurf    = double(wout.ns);         %!< number of radial surfaces
Geom.Nperiods = double(wout.nfp);        %!< number of field periods
if not(isfield(wout,'iasym'))
  wout.iasym=isfield(wout,'bmns');
end
Geom.handedness=handedness;
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

Pest.handedness=handedness;
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

Vmec.handedness=handedness;
Vmec.Nvmecu=Nu;
Vmec.Nvmecw=Nw;

Vmec.Rmnlist.m=double(wout.xm);
Vmec.Rmnlist.n=signchange*double(wout.xn)/Geom.Nperiods;
Vmec.Rmnlist.cosparity=ones(size(Vmec.Rmnlist.m));
Vmec.Rmnlist.data=(wout.rmnc(:,rindf_minus)+wout.rmnc(:,rindf_plus))/2;
Vmec.Rmn=mnmat(Vmec.Rmnlist,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');
Vmec.R=ifftmn(Vmec.Rmn,Geom.Nperiods,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');

Vmec.Zmnlist.m=double(wout.xm);
Vmec.Zmnlist.n=signchange*double(wout.xn)/Geom.Nperiods;
Vmec.Zmnlist.cosparity = 0*ones(size(Vmec.Zmnlist.m));
Vmec.Zmnlist.data=(wout.zmns(:,rindf_minus)+wout.zmns(:,rindf_plus))/2;
Vmec.Zmn=mnmat(Vmec.Zmnlist,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');
Vmec.Z=ifftmn(Vmec.Zmn,Geom.Nperiods,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');

Vmec.Bmnlist.m=double(wout.xm_nyq);
Vmec.Bmnlist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
Vmec.Bmnlist.cosparity=ones(size(Vmec.Bmnlist.m));
Vmec.Bmnlist.data=wout.bmnc(:,skrindh);
Vmec.Bmn=mnmat(Vmec.Bmnlist,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');
Vmec.B=ifftmn(Vmec.Bmn,Geom.Nperiods,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');

Vmec.Jmnlist.m=double(wout.xm_nyq);
Vmec.Jmnlist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
Vmec.Jmnlist.cosparity=ones(size(Vmec.Jmnlist.m));
Vmec.Jmnlist.data=wout.gmnc(:,skrindh) * signchange;
Vmec.Jmn=mnmat(Vmec.Jmnlist,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');
Vmec.J=ifftmn(Vmec.Jmn,Geom.Nperiods,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');

Vmec.B_ulist.m=double(wout.xm_nyq);
Vmec.B_ulist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
Vmec.B_ulist.cosparity=ones(size(Vmec.B_ulist.m));
Vmec.B_ulist.data=wout.bsubumnc(:,skrindh);
Vmec.B_umn=mnmat(Vmec.B_ulist,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');
Vmec.B_umntilde=remove00(Vmec.B_umn);

Vmec.B_wlist.m=double(wout.xm_nyq);
Vmec.B_wlist.n=signchange*double(wout.xn_nyq)/Geom.Nperiods;
Vmec.B_wlist.cosparity=ones(size(Vmec.B_wlist.m)); 
Vmec.B_wlist.data=wout.bsubvmnc(:,skrindh) * signchange;
Vmec.B_wmn=mnmat(Vmec.B_wlist,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');
Vmec.B_wmntilde=remove00(Vmec.B_wmn);

Vmec.llist.m=double(wout.xm);
Vmec.llist.n=signchange*double(wout.xn)/Geom.Nperiods;
Vmec.llist.cosparity = 0*ones(size(Vmec.llist.m));
Vmec.llist.data=wout.lmns(:,skrindh);
Vmec.lmn=mnmat(Vmec.llist,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');
Vmec.l=ifftmn(Vmec.lmn,Geom.Nperiods,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Pest coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vmec.Dvmecw=2*pi/Vmec.Nvmecw/Geom.Nperiods;
Vmec.Dvmecu=2*pi/Vmec.Nvmecu;
Vmec.vmecuvec=(0:Vmec.Nvmecu-1)*Vmec.Dvmecu;
Vmec.vmecwvec=(0:Vmec.Nvmecw-1)'*Vmec.Dvmecw;
[Vmec.vmecu,Vmec.vmecw] = ndgrid(Vmec.vmecuvec,Vmec.vmecwvec);

Vmec.Dpthetavmecu=ifftmn(Vmec.lmn,Geom.Nperiods,Vmec.Nvmecu,Vmec.Nvmecw,'forceSize');
Vmec.ptheta=Vmec.Dpthetavmecu+Vmec.vmecu; %This is the "transformation"
Vmec.pzeta=Vmec.vmecw;

Pest.Nptheta=Vmec.Nvmecu;
Pest.Npzeta=Vmec.Nvmecw;
Pest.Dptheta=Vmec.Dvmecu;
Pest.Dpzeta=Vmec.Dvmecw;
Pest.pzeta=Vmec.vmecw;
Pest.ptheta=Vmec.vmecu;

Pest.Dpthetavmecu=griddatacyclic(Vmec.ptheta,Vmec.pzeta,Vmec.Dpthetavmecu,Geom.Nperiods);
Pest.vmecu=Pest.ptheta-Pest.Dpthetavmecu; %This is the "transformation"
Pest.vmecw=Pest.pzeta;

Pest.B=interp2_cyclic(Vmec.vmecu,Vmec.vmecw,Vmec.B,Pest.vmecu,Pest.vmecw,Geom.Nperiods);
Pest.R=interp2_cyclic(Vmec.vmecu,Vmec.vmecw,Vmec.R,Pest.vmecu,Pest.vmecw,Geom.Nperiods);
Pest.Z=interp2_cyclic(Vmec.vmecu,Vmec.vmecw,Vmec.Z,Pest.vmecu,Pest.vmecw,Geom.Nperiods);
%The routine interp2_cyclic can also be used to interpolate other quantities from
%vmec to pest or vice versa.

Pest.mnmat.B=fftmn(Pest.B);
Pest.mnmat.R=fftmn(Pest.R);
Pest.mnmat.Z=fftmn(Pest.Z);

Pest.R00=mean(mean(Pest.R));
Pest.B00=mean(mean(Pest.B));



  
if 0
  fig(1)
  %surf(Vmec.u,Vmec.w,Vmec.Dzetaw);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  %surf(Vmec.u,Vmec.w,Vmec.zeta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  surf(Vmec.u,Vmec.w,Vmec.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])

  fig(2)
  %surf(Vmec.u,Vmec.w,Vmec.Dthetau);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  surf(Vmec.u,Vmec.w,Vmec.ptheta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  %surf(Vmec.u,Vmec.w,Vmec.u);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])

  %fig(5)
  %surf(Booz_u,Booz_w,Pest.Dzetaw);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  %surf(Booz_u,Booz_w,Pest.zeta);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  %surf(Booz_u,Booz_w,Pest.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])


  fig(3)
  surf(Vmec.ptheta,Vmec.pzeta,Vmec.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])
  fig(7)
  surf(Pest.ptheta,Pest.pzeta,Pest.B);shading flat;view(0,90);colorbar;axis([-0.5,6.5,-0.2,1.4])

end