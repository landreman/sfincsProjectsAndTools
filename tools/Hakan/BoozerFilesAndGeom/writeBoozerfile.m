function writeBoozerfile(Geom,filename,varargin)

if filename(end-2:end)=='.bc'
  filetype='bc';
else
  if filename(end-3:end)=='.dat'
    filetype='dat';
    if nargin>2
      Nradii=varargin{1}; %HM's W7X calculations usually has 7
    else
      Nradii=Geom.nsurf;
    end
  else
    filename=[filename,'.bc'];
    filetype='bc';
  end
end


if not(isfield(Geom,'newsigncorr'))
   error('Geom did not contain any field ''newsigncorr''!')
end
  
if strcmp(filetype,'dat') && not(Geom.StelSym)
  error(['Non-stellarator symmetric case not implemented yet. JG and ES might have ', ...
         'different standards there. I don''t know'])
end


fid = fopen(filename,'wt');
if fid<0
  error(['Unable to open file: ',filename])
end

% Change the signs from right to left handed (r,poloidal,toroidal)
signchange=-1; %sign changer
if isfield(Geom,'psi_a')
  Geom.psi_a=Geom.psi_a*signchange;
else
  disp('Geom did not contain any field ''psi_a''!')
end
if isfield(Geom,'torfluxtot')
  Geom.torfluxtot=Geom.torfluxtot*signchange;
else
  disp('Geom did not contain any field ''torfluxtot''!')
end
Geom.Bphi=Geom.Bphi*signchange^2;
Geom.Btheta=Geom.Btheta*signchange;

if not(Geom.newsigncorr)
  Geom.Bphi=-Geom.Bphi;
  Geom.Btheta=-Geom.Btheta;
else
  Geom.torfluxtot=-Geom.torfluxtot;
end

Geom.iota=Geom.iota*signchange;
if isfield(Geom,'dVdsoverNper')
  Geom.dVdsoverNper=Geom.dVdsoverNper*signchange;
else
  disp('Geom did not contain any field ''dVdsoverNper''!')
end

for tmpind=1:length(Geom.n)
  Geom.n{tmpind}=-Geom.n{tmpind};
end

if isfield(Geom,'Dphi')
for tmpind=1:length(Geom.n)
  Geom.Dphi{tmpind}=-Geom.Dphi{tmpind};
end
else
  disp('Geom did not contain any field ''Dphi''!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(filetype,'bc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %Begin writing to the file

  if isfield(Geom,'headertext')  
    fprintf(fid,'%s\n',Geom.headertext.maincomment);

    if isfield(Geom.headertext,'globalvars')
      fprintf(fid,'%s\n',Geom.headertext.globalvars);
    else
      if Geom.StelSym %JG includes a little more info in this case
        fprintf(fid,[' m0b  n0b nsurf nper  flux/[Tm^2]     a/[m]     R/[m]   avol/[m]  ' ...
                     'Rvol/[m] Vol/[m^3]\n']);
      else %and ES a little less
        fprintf(fid,' m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]\n');
      end
    end
  end

  minorradiusVMEC_exists=0;
  if isfield(Geom,'minorradiusVMEC')
    if not(isnan(Geom.minorradiusVMEC))
      minorradiusVMEC_exists=1;
    end
  end
  majorradiusVMEC_exists=0;
  if isfield(Geom,'majorradiusVMEC')
    if not(isnan(Geom.majorradiusVMEC))
      majorradiusVMEC_exists=1;
    end
  end
  VolumeVMEC_exists=0;
  if isfield(Geom,'VolumeVMEC')
    if not(isnan(Geom.VolumeVMEC))
      VolumeVMEC_exists=1;
    end
  end


  if isnan(Geom.majorradiusLastbcR00)
    Geom.majorradiusLastbcR00=0;
  end


  if Geom.StelSym
    if minorradiusVMEC_exists && not(majorradiusVMEC_exists) && not(VolumeVMEC_exists)
      fprintf(fid,'%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f\n',...
              Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,Geom.minorradiusVMEC);
    elseif minorradiusVMEC_exists && majorradiusVMEC_exists && not(VolumeVMEC_exists)
      fprintf(fid,'%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n',...
              Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,...
              Geom.minorradiusVMEC,Geom.majorradiusVMEC,...
              pi*Geom.minorradiusVMEC^2*2*pi*Geom.majorradiusVMEC);    
    elseif minorradiusVMEC_exists && majorradiusVMEC_exists && VolumeVMEC_exists
      fprintf(fid,'%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n',...
              Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,...
              Geom.minorradiusVMEC,Geom.majorradiusVMEC,Geom.VolumeVMEC);        
    else
      fprintf(fid,'%3d%5d%6d%4d%16.6E%10.5f%10.5f\n',...
              Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00);
    end
  else
    %Difficult! I Think Erika uses VMEC quantities here but I am not sure!
    if minorradiusVMEC_exists && majorradiusVMEC_exists
      fprintf(fid,'%4d%6d%7d%6d%15.8f%15.8f%15.8f\n',...
              Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
              Geom.minorradiusVMEC,Geom.majorradiusVMEC);
    else
      error('Did not find the minorradiusVMEC and majorradiusVMEC!')
      %fprintf(fid,'%4d%6d%7d%6d%15.8f%15.8f%15.8f\n',...
      %        Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
      %        Geom.minorradiusW7AS,Geom.majorradiusLastbcR00);  
    end
  end
  for rind=1:Geom.nsurf
    if isfield(Geom.headertext,'surfvars')
      fprintf(fid,'%s\n',Geom.headertext.surfvars);
    else
      if Geom.StelSym
        fprintf(fid,['       s         iota  curr_pol/nper    curr_tor    pprime   sqrt ' ...
                     'g(0,0)\n']);
      else
        fprintf(fid,['        s               iota           Jpol/nper          Itor',...
                     '            pprime         sqrt g(0,0)\n']);
      end
    end
    if isfield(Geom.headertext,'surfvarunits')
      fprintf(fid,'%s\n',Geom.headertext.surfvarunits);
    else
      if Geom.StelSym
        fprintf(fid,['                            [A]            [A]   dp/ds,[Pa] ' ...
                     '(dV/ds)/nper\n']);
      else
        fprintf(fid,['                                          [A]           [A] ',...
                     '            [Pa]         (dV/ds)/nper\n']);
      end
    end
    if Geom.StelSym
      fprintf(fid,'%12.4E%12.4E%12.4E%12.4E%12.4E%12.4E\n',...
              Geom.s(rind),Geom.iota(rind),...
              Geom.Bphi(rind)/(Geom.Nperiods/2/pi*(4*pi*1e-7)),...
              Geom.Btheta(rind)/(1/2/pi*(4*pi*1e-7)),...
              Geom.dpds(rind),Geom.dVdsoverNper(rind));
    else
      fprintf(fid,'%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E\n',...
              Geom.s(rind),Geom.iota(rind),...
              Geom.Bphi(rind)/(Geom.Nperiods/2/pi*(4*pi*1e-7)),...
              Geom.Btheta(rind)/(1/2/pi*(4*pi*1e-7)),...
              Geom.dpds(rind),Geom.dVdsoverNper(rind));
    end
    if isfield(Geom.headertext,'datavars')
      fprintf(fid,'%s\n',Geom.headertext.datavars);
    else
      if Geom.StelSym
        fprintf(fid,['    m    n        r/[m]           z/[m] (phib-phi)*nper/twopi     ',...
                     'bmn/[T]\n']);
      else
        fprintf(fid,['    m    n      rmnc [m]         rmns [m]         zmnc [m]         ',...
                     'zmns [m]         vmnc [ ]         vmns [ ]         ',...
                     'bmnc [T]         bmns [T]\n']);
      end
    end
    if Geom.StelSym
      for ind=1:Geom.nmodes(rind)
        fprintf(fid,'%5d%5d%16.8E%16.8E%16.8E%16.8E\n',...
                Geom.m{rind}(ind),Geom.n{rind}(ind),...
                Geom.R{rind}(ind),Geom.Z{rind}(ind),...
                Geom.Dphi{rind}(ind),Geom.Bmn{rind}(ind));
      end
    else     
      allms=sort(Geom.m{rind});
      ms=allms([1,find(diff(allms))+1]);
      Nm=length(ms);
      for mind=1:Nm
        m=ms(mind);
        allns=sort(Geom.n{rind}(find(Geom.m{rind}==m)));
        ns=allns([1,find(diff(allns))+1]);
        Nn=length(ns);
        for nind=1:Nn
          n=ns(nind);
          mninds=find((Geom.n{rind}==n)&(Geom.m{rind}==m));
          if length(mninds)>2 || isempty(mninds)
            error('This is impossible!')
          end
          rmnc=0;rmns=0;zmnc=0;zmns=0;
          vmnc=0;vmns=0;bmnc=0;bmns=0;
          for iind=1:length(mninds)
            if Geom.parity{rind}(mninds(iind))==1 %1<=>cosinus Bmn (stelsym part)
              rmnc=Geom.R{rind}(mninds(iind));
              zmns=Geom.Z{rind}(mninds(iind));
              vmns=Geom.Dphi{rind}(mninds(iind));
              bmnc=Geom.Bmn{rind}(mninds(iind));
            else %0<=>sinus Bmn (non-stelsym part)
              rmns=Geom.R{rind}(mninds(iind));
              zmnc=Geom.Z{rind}(mninds(iind));
              vmnc=Geom.Dphi{rind}(mninds(iind));
              bmns=Geom.Bmn{rind}(mninds(iind));            
            end
          end
          fprintf(fid,'%5d%5d%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E\n',...
                  m,n,rmnc,rmns,zmnc,zmns,vmnc,vmns,bmnc,bmns);        
        end
      end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else %filetype='dat'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if Nradii~=Geom.nsurf
    %dr=(Geom.rnorm(end)-Geom.rnorm(1))/Nradii;
    rnorms=linspace(Geom.rnorm(1),Geom.rnorm(end),Nradii+1);
    rnorms=rnorms(1:end-1);
    as=rnorms.^2;
    %dr=1/Nradii;
    %as=(([1:Nradii]-0.5)*dr).^2;
    %new_s(1)=max(Geom.s(1),new_s(1));
    %if new_s(2)<Geom.s(1)
    %  error('unforeseen')
    %end

    Geom=interpBoozer(Geom,'s',as,'linear');
    for rind=1:Nradii
      Geom.dBmndr{rind}=Geom.dBmnds{rind}*2*rnorms(rind)/Geom.minorradiusW7AS;
    end
  else
    Geom=interpBoozer(Geom,'s',Geom.s,'linear')
    for rind=1:Nradii
      Geom.dBmndr{rind}=Geom.dBmnds{rind}*2*sqrt(Geom.s(rind))/Geom.minorradiusW7AS;
    end
  end

  
  %Begin writing to the file
  
  globalvars_exist=0;
  if isfield(Geom,'headertext')
    fprintf(fid,'%s\n',Geom.headertext.maincomment);
    globalvars_exist= isfield(Geom.headertext,'globalvars');
  end
  if globalvars_exist
    fprintf(fid,'cc %s\n',Geom.headertext.globalvars);
  else
    if isfield(Geom,'m0b') && isfield(Geom,'n0b') && isfield(Geom,'torfluxtot')
      %if Geom.StelSym %JG includes a little more info in this case
      fprintf(fid,['cc  m0b  n0b nsurf nper  flux/[Tm^2]     a/[m]     R/[m]   avol/[m]  ' ...
                   'Rvol/[m] Vol/[m^3]\n']);
      %else %and ES a little less
      %  fprintf(fid,' m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]\n');
      %end
    else
      fprintf(fid,['cc  nsurf nper    a/[m]     R/[m]   avol/[m]  ' ...
                   'Rvol/[m] Vol/[m^3]\n']);
    end
  end


  minorradiusVMEC_exists=0;
  if isfield(Geom,'minorradiusVMEC')
    if not(isnan(Geom.minorradiusVMEC))
      minorradiusVMEC_exists=1;
    end
  end
  majorradiusVMEC_exists=0;
  if isfield(Geom,'majorradiusVMEC')
    if not(isnan(Geom.majorradiusVMEC))
      majorradiusVMEC_exists=1;
    end
  end
  VolumeVMEC_exists=0;
  if isfield(Geom,'VolumeVMEC')
    if not(isnan(Geom.VolumeVMEC))
      VolumeVMEC_exists=1;
    end
  end

  if isnan(Geom.majorradiusLastbcR00)
    Geom.majorradiusLastbcR00=0;
  end

  if isfield(Geom,'m0b') && isfield(Geom,'n0b') && isfield(Geom,'torfluxtot')
    if minorradiusVMEC_exists && not(majorradiusVMEC_exists) && not(VolumeVMEC_exists)
      fprintf(fid,'cc %3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f\n',...
              Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,Geom.minorradiusVMEC);
    elseif minorradiusVMEC_exists && majorradiusVMEC_exists && not(VolumeVMEC_exists)
      fprintf(fid,'cc %3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n',...
              Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,...
              Geom.minorradiusVMEC,Geom.majorradiusVMEC,...
              pi*Geom.minorradiusVMEC^2*2*pi*Geom.majorradiusVMEC);    
    elseif minorradiusVMEC_exists && majorradiusVMEC_exists && VolumeVMEC_exists
      fprintf(fid,'cc %3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n',...
              Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,...
              Geom.minorradiusVMEC,Geom.majorradiusVMEC,Geom.VolumeVMEC);        
    else
      fprintf(fid,'cc %3d%5d%6d%4d%16.6E%10.5f%10.5f\n',...
              Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00);
    end
  else
    if minorradiusVMEC_exists && not(majorradiusVMEC_exists) && not(VolumeVMEC_exists)
      fprintf(fid,'cc %6d%4d%10.5f%10.5f%11.5f\n',...
              Geom.nsurf,Geom.Nperiods,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,Geom.minorradiusVMEC);
    elseif minorradiusVMEC_exists && majorradiusVMEC_exists && not(VolumeVMEC_exists)
      fprintf(fid,'cc %6d%4d%10.5f%10.5f%11.5f%11.5f%11.5f\n',...
              Geom.nsurf,Geom.Nperiods,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,...
              Geom.minorradiusVMEC,Geom.majorradiusVMEC,...
              pi*Geom.minorradiusVMEC^2*2*pi*Geom.majorradiusVMEC);    
    elseif minorradiusVMEC_exists && majorradiusVMEC_exists && VolumeVMEC_exists
      fprintf(fid,'cc %6d%4d%10.5f%10.5f%11.5f%11.5f%11.5f\n',...
              Geom.nsurf,Geom.Nperiods,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,...
              Geom.minorradiusVMEC,Geom.majorradiusVMEC,Geom.VolumeVMEC);        
    else
      fprintf(fid,'cc %6d%4d%10.5f%10.5f\n',...
              Geom.nsurf,Geom.Nperiods,...
              Geom.minorradiusW7AS,Geom.majorradiusLastbcR00);
    end    
  end
  fprintf(fid,' c     m, n, B_mn/B_00, B_mn, (dB_mn/dr)/B_00 [1/cm], R_mn, z_mn\n');
  for rind=1:Geom.nsurf

    fprintf(fid,'%6.2f%8.4f%3i%8.2f%8.2f%7.4f%8.2f%8.2f  r,io,N,ra,Rm,s,r_o,R00\n',...
            Geom.minorradiusW7AS*Geom.rnorm(rind)*100,...
            Geom.iota(rind),Geom.Nperiods,Geom.minorradiusW7AS*100,...
            Geom.majorradiusLastbcR00*100,Geom.s(rind),...
            Geom.minorradiusW7AS*Geom.rnorm(rind)*100,Geom.R00(rind)*100);
    
    %Bphi(rind)  =surfheader2(1)*1e6*Nperiods/2/pi*(4*pi*1e-7); %Tesla*meter
    %Btheta(rind)=surfheader2(2)*1e6/2/pi*(4*pi*1e-7);          %Tesla*meter
    %B00(rind)   =surfheader2(3); %Tesla
    fprintf(fid,' >  %13.4E%13.4E%13.4E     cur_pol=g, cur_tor=I, B_00\n',...
            Geom.Bphi(rind)/(1e6*Geom.Nperiods/2/pi*(4*pi*1e-7)),...
            Geom.Btheta(rind)/(1e6/2/pi*(4*pi*1e-7)),...
            Geom.B00(rind));
    %        tmp=sscanf(tmp_str,'%d %d %f %f %f %f %f',7);
    %    if (abs(tmp(3))>min_Bmn)&&(tmp(1)<=max_m)&&(abs(tmp(2))<=maxabs_n)
    %      modeind=modeind+1;
    %      modesm{rind}(modeind)=tmp(1);
    %      modesn{rind}(modeind)=tmp(2);
    %      modesbnorm{rind}(modeind)=tmp(3);
    %      modesb{rind}(modeind)=tmp(4);
    %      modesdbdrnorm{rind}(modeind)=tmp(5)*100; %convert cm^-1 -> m^-1
    %      modesr{rind}(modeind)=tmp(6);
    %      modesz{rind}(modeind)=tmp(7);
    %    end

    m0inds=find(Geom.m{rind}==0);
    [dum,sorti]=sort(Geom.n{rind}(m0inds));
    for ind=1:length(m0inds)
      fprintf(fid,'%4i%4i%11.6f%11.6f%12.3E%11.5f%11.5f\n', ...
              Geom.m{rind}(m0inds(sorti(ind))),Geom.n{rind}(m0inds(sorti(ind))),...
              Geom.Bnorm{rind}(m0inds(sorti(ind))),...
              Geom.Bmn{rind}(m0inds(sorti(ind))),...
              Geom.dBmndr{rind}(m0inds(sorti(ind)))/Geom.B00(rind)/100,...
              Geom.R{rind}(m0inds(sorti(ind)))*100,...
              Geom.Z{rind}(m0inds(sorti(ind)))*100);
    end
    for m=1:max(Geom.m{rind})
      thisminds=find(Geom.m{rind}==m);
      [dum,sorti]=sort(Geom.n{rind}(thisminds),'descend');
      for ind=1:length(thisminds)
        fprintf(fid,'%4i%4i%11.6f%11.6f%12.3E%11.5f%11.5f\n',...
                Geom.m{rind}(thisminds(sorti(ind))),Geom.n{rind}(thisminds(sorti(ind))),...
                Geom.Bnorm{rind}(thisminds(sorti(ind))),...
                Geom.Bmn{rind}(thisminds(sorti(ind))),...
                Geom.dBmndr{rind}(thisminds(sorti(ind)))/Geom.B00(rind)/100,...
                Geom.R{rind}(thisminds(sorti(ind)))*100,...
                Geom.Z{rind}(thisminds(sorti(ind)))*100)    ;
      end
    end
    
    fprintf(fid,'  -1   0 0.0 0.0 0.0 0.0 0.0     end of input\n');
  end
end
fclose(fid);