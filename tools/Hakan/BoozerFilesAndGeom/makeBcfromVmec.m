function Geom=makeBcfromVmec(woutin,Nu,Nw,min_Bmn)
if nargin==1
  Nu=inf;
  Nw=inf;
  min_Bmn=0;
end

%The result can be saved with writeBoozerfile.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the .nc file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%woutin is the wout file name or just the netcdf variables from the wout file
%(Too old matlab versions do not have the necessary netcdf routines.)

if isstruct(woutin)
  wout=woutin;
else %if not already loaded, assume woutin is a string with the file name
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% let w = -v!

signchange=double(wout.signgs); % is -1, because vmec is left handed

Geom.headertext.input_extension = ...
    wout.input_extension; %!<  suffix of the vmec-input file: input.<input_extension>

StelSym=not(isfield(wout,'bmns'));

if StelSym
  % Use Joachim Geiger's JMC convention as default for stellarator symmetric files:
  % Toroidal flux is counted positive in the direction opposite to the toroidal coordinate.
  Geom.headertext.maincomment=sprintf(...
    ['CC Converted from VMEC using HS''s matlab routines\n',...
     'CC The phase convention (m theta - n phi) is used in this file.\n',...
     'CC The coordinate system (r,theta,phi) is left-handed.\n',...
     'CC Toroidal flux is counted positive in the direction opposite to the toroidal ', ...
     'coordinate (JG Style).']);
else  
  % Use Erika Strumberger's convention as default for non-stellarator symmetric files:
  % Toroidal flux is counted positive in the direction of the toroidal coordinate.
  Geom.headertext.maincomment=sprintf(...
    ['CC Converted from VMEC using HS''s matlab routines\n',...
     'CC The phase convention (m theta - n phi) is used in this file.\n',...
     'CC The coordinate system (r,theta,phi) is left-handed.\n',...
     'CC Toroidal flux is counted positive in the direction of the toroidal coordinate (ES Style).']);
  % Note: the comment that (r,theta,phi) is left-handed regards the 
  % file that is created when saving with writeBoozerfile.m. The struct
  % Geom produced by this function (makeBcfromVMEC) is right handed.
end

Geom.m0b      = NaN;%double(wout.mpol);       %!< number of poloidal fourier modes (Not used)
Geom.n0b      = NaN;%double(wout.ntor);       %!< upper bound toroidal fourier modes:
                                         %-ntor <= n <= ntor (Not used)
Geom.nsurf    = NaN;%double(wout.ns);    %!< number of radial surfaces, Set this below

Geom.Nperiods = double(wout.nfp);        %!< number of field periods
%if not(isfield(wout,'iasym'))
%  wout.iasym=isfield(wout,'bmns');
%end
%Geom.StelSym  = not(wout.iasym);%!<  defines stellarator symmetry for iasym=0,otherwise =

Geom.StelSym  = StelSym;
if StelSym
  Geom.newsigncorr=0; %This is the convention in SFINCS
end
Geom.torfluxtot = ...
   wout.phi(wout.ns)*signchange;%!<  total toroidal flux within the boundary (s=1)
Geom.psi_a=Geom.torfluxtot/2/pi;

Geom.minorradiusVMEC= wout.Aminor_p;  %!<  minor plasma radius
Geom.majorradiusLastbcR00 = NaN; %Calculate this below (not necessary)
Geom.minorradiusW7AS= NaN; %Calculate this below (not necessary)
Geom.majorradiusVMEC= wout.Rmajor_p;  %!<  major plasma radius

Geom.fullgrid.s     = wout.phi'/wout.phi(wout.ns);     %full grid
Geom.fullgrid.rnorm = sqrt(Geom.fullgrid.s);     %full grid

skip=1; %this is how many elements are skipped at low radii when going to half grid

Geom.s=(Geom.fullgrid.s(skip:end-1)+Geom.fullgrid.s(skip+1:end))/2; %half grid
Geom.rnorm=sqrt(Geom.s); %half grid
Geom.nsurf=length(Geom.s);

Geom.dpds=diff(wout.presf(skip:end))'./diff(Geom.fullgrid.s(skip:end));

%Geom.dVdsoverNper=dVdsoverNper*signchange;
Geom.Bphi  = wout.bvco*signchange;%direction switch sign
Geom.Btheta= wout.buco;%*signchange;

Geom.Bphi=Geom.Bphi(skip+1:end)';
Geom.Btheta=Geom.Btheta(skip+1:end)';

Geom.fullgrid.iota = wout.iotaf'*signchange; 
Geom.iota = wout.iotas'*signchange; %half mesh
Geom.iota = Geom.iota(skip+1:end);
Geom.Bfilter.min_Bmn=min_Bmn;

%%% Choose a suitable spatial resoulition
if Nu==inf
  %Use the Vmec mpol and mtor
  Nu = 1+2*max(abs(double(wout.xm)));
else
  Nu=2*floor(Nu/2)+1; %force to be odd
end
if Nw==inf
  %Use the Vmec mpol and mtor
  Nw=1+2*max(abs(double(wout.xn))); 
else
  Nw=2*floor(Nw/2)+1; %force to be odd
end
axisymm=0;
if Nw==1
  axisymm=1;
  Nw=3; %The subroutines cannot handle Nw=1
end

Geom.m0b=(Nu-1)/2;
Geom.n0b=(Nw-1)/2;
Geom.Bfilter.max_m=(Nu-1)/2;
Geom.Bfilter.maxabs_n=(Nw-1)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specific for each radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if Geom.StelSym
%  error('Stellarator symmetric case not implemented yet!')
%end

tic
for sind=1:length(Geom.s)
  fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b radius %i / %i',...
          sind,length(Geom.s));
  Booz=makeBoozfromVmec(woutin,Geom.s(sind),Nu,Nw);
  fprintf(1,' ,  %3.1f s',round(10*toc)/10);
  Bmnlist=mnlist(Booz.mnmat.B);
  Rmnlist=mnlist(Booz.mnmat.R);
  Zmnlist=mnlist(Booz.mnmat.Z);
  Dzetawmnlist=mnlist(Booz.mnmat.Dzetaw);
  
  if Geom.StelSym
    Geom.nmodes(1,sind)=length(Bmnlist.m);
    Geom.m{sind}  =Bmnlist.m;
    Geom.n{sind}  =Bmnlist.n;
    Geom.parity{sind}=Bmnlist.cosparity;
    Geom.Bmn{sind}=Bmnlist.data;
    Geom.R{sind}  =Rmnlist.data;
  else
    Geom.nmodes(1,sind)=length(Bmnlist.m)+1;
    Geom.m{sind}  =[Bmnlist.m,0];
    Geom.n{sind}  =[Bmnlist.n,0];
    Geom.parity{sind}=[Bmnlist.cosparity,0]; %Add a m=n=0 parity=0 (sin) component for Z00! 
    Geom.Bmn{sind}=[Bmnlist.data,0]; %Set to 0 rather than NaN, to avoid NaNs in bc file.
    Geom.R{sind}  =[Rmnlist.data,0]; %Set to 0 rather than NaN, to avoid NaNs in bc file.
  end
  Geom.B00(1,sind)=Booz.B00;
  Geom.Bnorm{sind}=Geom.Bmn{sind}/Geom.B00(sind);
  Geom.R00(1,sind)=Booz.R00;
  if not(Geom.StelSym)
    Geom.Z00(1,sind)=Booz.Z00;
  end
  Geom.Z{sind}=0*ones(size(Geom.R{sind}));
  for mnind=1:length(Zmnlist.m)
    ind=find(Geom.m{sind}==Zmnlist.m(mnind) & ...
             Geom.n{sind}==Zmnlist.n(mnind) & ...
             Geom.parity{sind}==not(Zmnlist.cosparity(mnind)));
    Geom.Z{sind}(ind)=Zmnlist.data(mnind);
  end
  
  Geom.Dphi{sind}=0*ones(size(Geom.R{sind}));
  for mnind=1:length(Dzetawmnlist.m)
    ind=find(Geom.m{sind}==Dzetawmnlist.m(mnind) & ...
             Geom.n{sind}==Dzetawmnlist.n(mnind) & ...
             Geom.parity{sind}==not(Dzetawmnlist.cosparity(mnind)));
    Geom.Dphi{sind}(ind)=Geom.Nperiods/2/pi*Dzetawmnlist.data(mnind);
  end
  
  %if Geom.StelSym
    ind=find(Geom.m{sind}==0 & Geom.n{sind}==0 & Geom.parity{sind}==1);
    Geom.Dphi{sind}(ind)=0; 
    Geom.Z{sind}(ind)=0; %These are NaN, but I set them to zero to avoid making .bc
                         %files with NaNs in them. 
  %end
  
  if axisymm
    good=find(Geom.n{sind}==0);
    Geom.nmodes(1,sind)=length(good);
    Geom.m{sind}=Geom.m{sind}(good);
    Geom.n{sind}=Geom.n{sind}(good);
    Geom.parity{sind}=Geom.parity{sind}(good);
    Geom.Bmn{sind}=Geom.Bmn{sind}(good);
    Geom.Bnorm{sind}=Geom.Bnorm{sind}(good);
    Geom.R{sind}=Geom.R{sind}(good);
    Geom.Z{sind}=Geom.Z{sind}(good);
    Geom.Dphi{sind}=Geom.Dphi{sind}(good);   
  end

  %Filter Bmns:
  good=find(abs(Geom.Bmn{sind})>min_Bmn);
  Geom.nmodes(1,sind)=length(good);
  Geom.m{sind}=Geom.m{sind}(good);
  Geom.n{sind}=Geom.n{sind}(good);
  Geom.parity{sind}=Geom.parity{sind}(good);
  Geom.Bmn{sind}=Geom.Bmn{sind}(good);
  Geom.Bnorm{sind}=Geom.Bnorm{sind}(good);
  Geom.R{sind}=Geom.R{sind}(good);
  Geom.Z{sind}=Geom.Z{sind}(good);
  Geom.Dphi{sind}=Geom.Dphi{sind}(good);   
  
  
  Geom.FSAB2(1,sind)=Booz.FSAB2; %=Nu*Nw/sum(sum(1./Booz.B.^2));
  
  Geom.mnmat{sind}=Booz.mnmat; %note that for axisymmetry I still keep Nw=3 here!
  
  Geom.dVdsoverNper(1,sind)=...
      abs(Geom.psi_a*4*pi^2/Geom.Nperiods*(Booz.G+Booz.iota.*Booz.I)/Booz.FSAB2);
end

if length(Geom.R00)==length(Geom.s) %just check that all were made.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate minorradiusW7AS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  accum=0;
  for m=1:wout.xm(end)
    mind=m+1;
    ii=find(wout.xm==m);
    start_modes=ii(1);
  end_modes=ii(end);
  ns=wout.xn(start_modes:end_modes)/double(wout.nfp);
  nnmat=(1+(-1).^(ns*ones(size(ns'))-ones(size(ns))*ns'))/2;
  accum=accum+...
        m*sum(sum((wout.rmnc(start_modes:end_modes,end)*...
                   wout.zmns(start_modes:end_modes,end)').*nnmat));
  end
  Geom.minorradiusW7AS=sqrt(abs(accum));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate majorradiusLastbcR00
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % We need to obtain the Boozer coordinates 
  % for the last half mesh flux surface
  % because Joachim defined his major radius to be 
  % R00 in Boozer coordinates at s=0.995
  %lastbcsurfind=find(abs(Geom.s-0.995)<5e-4);
  %maybe it is better to really take the last half mesh radius
  lastbcsurfind=length(Geom.s);

  if not(isempty(lastbcsurfind))
    Geom.majorradiusLastbcR00=Geom.R00(lastbcsurfind);
  else
    %We must construct the Boozer coordinates ourselfves!
    s_wish=0.995;
    Nu=121;
    Nv=121;
    Booz=makeBoozfromVmec(wout,s_wish,Nu,Nv);
    Geom.majorradiusLastbcR00=Booz.R00;
  end
end
fprintf(1,'\n')
%Geom.dVdsoverNper=