function [Geom,xf]=readBoozerXformfile(ncfile,min_Bmn,max_m,maxabs_n,symmetry)
switch nargin
 case 1
  min_Bmn=0;
  max_m=inf;
  maxabs_n=inf;
  symmetry='unknown';
 case 2
  max_m=inf;
  maxabs_n=inf;
  symmetry='unknown';
 case 3
  maxabs_n=inf;
  symmetry='unknown';
 case 4
  symmetry='unknown';
end

ver=version;
ind20=findstr(ver,'20');
if str2num(ver(ind20+2:ind20+3))<11
  error(['You must use a matlab version >2011 for netcdf to work!',...
         '2012b seems to be the only one that works on the hgw cluster!'])
end

%%% Read the netcdf file:
info=ncinfo(ncfile);
Nvars=length(info.Variables);
xf=struct();
xf.variableNames={info.Variables.Name};
xf.Nvariables=Nvars;

for vari=1:Nvars
  thisname=info.Variables(vari).Name;
  thisval=ncread(ncfile,['/',thisname]);
  xf=setfield(xf,thisname,thisval);
end

Geom.headertext.maincomment=sprintf('CC Converted from the Xform file \nCC %s',ncfile);
Geom.m0b      = double(xf.mboz_b)/2;
Geom.n0b      = double(xf.nboz_b);
Geom.nsurf    = length(xf.jlist);
Geom.Nperiods = double(xf.nfp_b);
Geom.psi_a    = xf.phi_b(end)/2/pi; %Note: the sign change is done below
Geom.torfluxtot=xf.phi_b(end); %Note: the sign change is done below
Geom.minorradiusW7AS=NaN;
Geom.majorradiusLastbcR00=NaN;
Geom.minorradiusVMEC=NaN;
Geom.majorradiusVMEC=NaN;
Geom.VolumeVMEC=NaN;
Geom.Bfilter.min_Bmn=min_Bmn;
Geom.Bfilter.max_m=max_m;
Geom.Bfilter.maxabs_n=maxabs_n;

Geom.m0b=min(Geom.m0b,max_m);
Geom.n0b=min(Geom.n0b,maxabs_n);



if not(strcmp(symmetry,'StelSym') || strcmp(symmetry,'unknown'))
  error('So far, this routine only treats stellarator symmetric cases!')
end
Geom.StelSym  =1;
sfull=xf.phi_b'/xf.phi_b(end);
shalf=(sfull(2:end)+sfull(1:end-1))/2;
%s = xf.phi_b(xf.jlist')/xf.phi_b(end);
s=shalf(xf.jlist'-1);
Geom.rnorm=sqrt(s);                     %on half mesh
Geom.s=s;                               %on half mesh
Geom.Bphi     = xf.bvco_b(xf.jlist)';   %on half mesh %Note: the sign change is done below
Geom.Btheta   = xf.buco_b(xf.jlist)';   %on half mesh
Geom.iota     = xf.iota_b(xf.jlist)';   %on half mesh
p=xf.pres_b(xf.jlist)';                 %on half mesh

mu0= 4*pi*1e-7;
dpds=[(p(2)-p(1))/(s(2)-s(1)),...
      (p(3:end)-p(1:end-2))./(s(3:end)-s(1:end-2)),...
      (p(end)-p(end-1))/(s(end)-s(end-1))];
Geom.dpds=dpds;                     %on half mesh  
%Geom.p=p;                          %on half mesh

Geom.dVdsoverNper=NaN*zeros(size(Geom.s));
%Geom.dVdsoverNper2=NaN*zeros(size(Geom.s));

Geom.B00=zeros(size(Geom.s));
Geom.R00=zeros(size(Geom.s));

%First, read all the data:

m0n0ind=find(xf.ixm_b==0 & xf.ixn_b==0);
for rind=1:Geom.nsurf
  Geom.nmodes(rind) = double(xf.mnboz_b);
  Geom.m{rind}      = double(xf.ixm_b)';
  Geom.n{rind}      = double(xf.ixn_b)'/Geom.Nperiods; %Note: the sign change is done below
  Geom.Bmn{rind}    = xf.bmnc_b(:,rind)';
  B00(rind)         = xf.bmnc_b(m0n0ind,rind);
  Geom.Bnorm{rind}  = Geom.Bmn{rind}/B00(rind);
  Geom.B00(rind)    = B00(rind);
  Geom.R00(rind)    = xf.rmnc_b(m0n0ind,rind);
  Geom.R{rind}      = xf.rmnc_b(:,rind)'; 
  Geom.Z{rind}      = xf.zmns_b(:,rind)'; 
  Geom.Dphi{rind}   = xf.pmns_b(:,rind)'...
                        *Geom.Nperiods/2/pi; %Note: the sign change is done below
  Geom_Jac{rind}    = xf.gmn_b(:,rind)';  %cos terms for g=Jacobian 
                                          %NB! g is not Jacobian^2 here
  Geom.dVdsoverNper(rind)=...
      -abs(4*pi^2*Geom_Jac{rind}(m0n0ind)*...
           Geom.psi_a/Geom.Nperiods.*Geom.iota(rind)); %This should be negative in a LH system
end

lastbcsurfind=find(abs(Geom.s-0.995)<5e-4);
if not(isempty(lastbcsurfind))
  Geom.majorradiusLastbcR00=Geom.R00(lastbcsurfind);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change the signs from left to right handed (r,poloidal,toroidal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signchange=-1; %sign changer
Geom.psi_a=Geom.psi_a*signchange;
Geom.torfluxtot=Geom.torfluxtot*signchange;
Geom.Bphi=Geom.Bphi*signchange;
Geom.Btheta=Geom.Btheta;
Geom.iota=Geom.iota*signchange;
Geom.dVdsoverNper=Geom.dVdsoverNper*signchange;
for tmpind=1:length(Geom.n)
  Geom.n{tmpind}=-Geom.n{tmpind};
  Geom.Dphi{tmpind}=-Geom.Dphi{tmpind}; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turn the Stellarator around to make the toroidal flux positive
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Geom.psi_a     =-Geom.psi_a;
Geom.torfluxtot=-Geom.torfluxtot;
Geom.Bphi      =-Geom.Bphi;
Geom.Btheta    =-Geom.Btheta;
for tmpind=1:length(Geom.n)
  Geom.Dphi{tmpind}=-Geom.Dphi{tmpind}; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now do the filtering:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for rind=1:Geom.nsurf
  m0inds=find(Geom.m{rind}==0);
  Geom.n{rind}(m0inds)=-abs(Geom.n{rind}(m0inds)); %This is just a convention
  
  good=find(Geom.m{rind}<=max_m & abs(Geom.n{rind})<=maxabs_n & ...
            abs(Geom.Bmn{rind})>=min_Bmn);
  
  good2=[];
  for m=0:max(Geom.m{rind}(good))
    minds=find(Geom.m{rind}(good)==m);
    if not(isempty(minds))
      [dummy,ii]=sort(Geom.n{rind}(good(minds)),2,'ascend');
      good2=[good2,good(minds(ii))];
    end
  end 
  
  Geom.nmodes(rind)= length(good2);
  Geom.m{rind}     = Geom.m{rind}(good2);
  Geom.n{rind}     = Geom.n{rind}(good2);
  Geom.Bmn{rind}   = Geom.Bmn{rind}(good2);
  Geom.Bnorm{rind} = Geom.Bnorm{rind}(good2);
  Geom.R{rind}     = Geom.R{rind}(good2);
  Geom.Z{rind}     = Geom.Z{rind}(good2);
  Geom.Dphi{rind}  = Geom.Dphi{rind}(good2);  
  %error('stop')
end
