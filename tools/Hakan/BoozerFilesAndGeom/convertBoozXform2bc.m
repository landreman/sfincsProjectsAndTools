function convertBoozXform2bc(boozXfilename,bcfilename,aHat_or_woutfilename,varargin)

if nargin>3
  [Geomxf,xf]=readBoozerXformfile(boozXfilename,varargin{:});
else
  [Geomxf,xf]=readBoozerXformfile(boozXfilename);
end

%Throw away the first radius
if Geomxf.s(1)<0.009
  Geomxf=interpBoozer(Geomxf,'s',Geomxf.s(2:end),'nearest');
end

%Geomxf
%Geomxf_Bphi=Geomxf.Bphi(1:5)
%Geomxf_Btheta=Geomxf.Btheta(1:5)
%Geomxf_iota=Geomxf.iota(1:5)


if not(ischar(aHat_or_woutfilename))
  Geomxf.minorradiusW7AS=aHat_or_woutfilename;
  
  lastbcsurfind=find(abs(Geomxf.s-0.995)<5e-4);
  if isempty(lastbcsurfind)
    Geomxf.majorradiusLastbcR00=NaN;%Unknown but unimportant!
  else
    Geomxf.majorradiusLastbcR00=Geomxf.R00(lastbcsurfind);
  end
else
  woutfilename=aHat_or_woutfilename;
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %read the VMEC file
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  info=ncinfo(woutfilename);
  Nvars=length(info.Variables);
  wout=struct();
  wout.variableNames={info.Variables.Name};
  wout.Nvariables=Nvars;
  for vari=1:Nvars
    thisname=info.Variables(vari).Name;
    thisval=ncread(woutfilename,['/',thisname]);
    wout=setfield(wout,thisname,thisval);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calvulate minorradiusW7AS
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
  Geomxf.minorradiusW7AS=sqrt(abs(accum));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calvulate majorradiusLastbcR00
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % We need to obtain the Boozer coordinates 
  % for the last half mesh flux surface
  % because Joachim defined his major radius to be 
  % R00 in Boozer coordinates at s=0.995
  lastbcsurfind=find(abs(Geomxf.s-0.995)<5e-4);
  if not(isempty(lastbcsurfind))
    Geomxf.majorradiusLastbcR00=Geomxf.R00(lastbcsurfind);
  else
    %We must construct the Boozer coordinates ourselfves!
    s_wish=0.995;
    Nu=128;
    Nv=128;
    Booz=makeBoozfromVmec(wout,s_wish,Nu,Nv);
    Geomxf.majorradiusLastbcR00=Booz.R00;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Stor VMEC defined radii
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Geomxf.minorradiusVMEC=wout.Aminor_p;
  Geomxf.majorradiusVMEC=wout.Rmajor_p;
  Geomxf.VolumeVMEC=wout.volume_p;
end

disp(['Writing the file ',bcfilename]);
writeBoozerfile(Geomxf,bcfilename);