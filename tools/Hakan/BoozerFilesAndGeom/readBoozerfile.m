function Geom=readBoozerfile(filename,min_Bmn,max_m,maxabs_n,symmetry,varargin)

newsigncorrectionmethod=NaN; %This is the default from 20160610

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
 case 6
  if cmpstr(varargin{1},'signcorr1')
    newsigncorrectionmethod=0;
  elseif cmpstr(varargin{1},'signcorr2')
    newsigncorrectionmethod=1;
  else
    error('Sign correction method for JG files not recognised!')
  end
end

fid = fopen(filename);
if fid<0
 error(['Unable to open file: ',filename])
end

if filename(end-3:end)=='.dat'
  %Henning Maassberg type file
  filetype='HM';
elseif filename(end-2:end)=='.bc'
  %Joachim Geiger type file
  filetype='JG';
else %default to JG
  %Joachim Geiger type file
  filetype='JG';
end

if strcmp(filetype,'JG')

  tmp_str=fgetl(fid);
  concat_str=tmp_str;
  while tmp_str(1:2)=='CC'
    tmp_str=fgetl(fid);
    if tmp_str(1:2)=='CC'
      concat_str = sprintf([concat_str,'\n',tmp_str]); %comment line
    end
  end
  Geom.headertext.maincomment=concat_str;
  Geom.headertext.globalvars=tmp_str;
  header_d=fscanf(fid,'%d',4);
  %header_f=fscanf(fid,'%f',3);
  tmpstr=fgetl(fid);
  header_f=sscanf(tmpstr,'%f');
  %header_f=fscanf(fid,'%f',3);
  %fgetl(fid); %just skip the return character
  Geom.headertext.surfvars=fgetl(fid);  %Variable name line

  Geom.m0b        = header_d(1);
  Geom.n0b        = header_d(2);
  Geom.nsurf      = header_d(3);
  Geom.Nperiods   = header_d(4);
  Geom.psi_a=NaN; %Insert psi_a at this place in the list, but set it later.
  Geom.torfluxtot = header_f(1); %Note that this is not per pol. angle
  Geom.minorradiusW7AS        = header_f(2); 
  Geom.majorradiusLastbcR00 = header_f(3);
  if length(header_f)>3
    Geom.minorradiusVMEC=header_f(4);
    if length(header_f)>4
      Geom.majorradiusVMEC=header_f(5);
      if length(header_f)>5
        Geom.VolumeVMEC=header_f(6);
      end
    end
  end
  
  eof=0;
  rind=0;
  
  while not(eof)
    rind=rind+1;
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b%5i/%5i',rind,Geom.nsurf)
    Geom.headertext.surfvarunits=fgetl(fid);
    surfheader=fscanf(fid,'%f',6);
    torfluxnorm(rind) = surfheader(1);
    iota(rind)  = surfheader(2);
    Bphi(rind)  = surfheader(3)*Geom.Nperiods/2/pi*(4*pi*1e-7); %Tesla*meter
    Btheta(rind)= surfheader(4)/2/pi*(4*pi*1e-7);               %Tesla*meter
    dpds(rind)  = surfheader(5);
    dVdsoverNper(rind)  = surfheader(6);
    fgetl(fid); %just skip the return character
    tmpstr=fgetl(fid); %units line
    if rind==1
      Geom.headertext.datavars=tmpstr;
      Geom.StelSym=isempty(strfind(Geom.headertext.datavars, 'rmnc'));
      if not(strcmp(symmetry,'unknown'))
        if Geom.StelSym && not(strcmp(symmetry,'StelSym'))
          error(['Boozer file is stellarator symmetric, but input to readBoozerfile ', ...
                 'said it should be non-stellarator symmetric!'])
        elseif not(Geom.StelSym) && strcmp(symmetry,'StelSym')
          error(['Boozer file is non-stellarator symmetric, but input to readBoozerfile ', ...
                 'said it should be stellarator symmetric!'])          
        end
      end
    end
    position=ftell(fid);
    tmp_str=fgetl(fid);
    if Geom.StelSym
      tmp=sscanf(tmp_str,'%d %d %f %f %f %f',6);
      while not(tmp(1)==0 && tmp(2)==0)
        tmp_str=fgetl(fid);
        tmp=sscanf(tmp_str,'%d %d %f %f %f %f',6);    
      end
      B00(rind)=tmp(6);
      R00(rind)=tmp(3);
    else
      tmp=sscanf(tmp_str,'%d %d %f %f %f %f %f %f %f %f',10);
      while not(tmp(1)==0 && tmp(2)==0)
        tmp_str=fgetl(fid);
        tmp=sscanf(tmp_str,'%d %d %f %f %f %f %f %f %f %f',10);
      end
      B00(rind)=tmp(9);
      R00(rind)=tmp(3);     
    end
    fseek(fid,position,'bof'); %Rewind to beginning of surface data
    
    proceed=1;
    modeind=0;
    while proceed
      tmp_str=fgetl(fid);
      if length(tmp_str)==1
        if tmp_str==-1 %End of file has been reached
          proceed=0;
          eof=1;
        end
      elseif not(isempty(find(tmp_str=='s'))) %Next flux surface has been reached
        proceed=0;
      else
        if Geom.StelSym
          tmp=sscanf(tmp_str,'%d %d %f %f %f %f',6);
          if (abs(tmp(6))/B00(rind)>min_Bmn)&&(tmp(1)<=max_m)&&(abs(tmp(2))<=maxabs_n)
            modeind=modeind+1;
            modesm{rind}(modeind)=tmp(1);
            modesn{rind}(modeind)=tmp(2);
            modesr{rind}(modeind)=tmp(3);
            modesz{rind}(modeind)=tmp(4);
            modesp{rind}(modeind)=tmp(5);
            modesb{rind}(modeind)=tmp(6);
            modesbnorm{rind}(modeind)=tmp(6)/B00(rind);
          end
        else
          tmp=sscanf(tmp_str,'%d %d %f %f %f %f %f %f %f %f',10);
          if (tmp(1)<=max_m)&&(abs(tmp(2))<=maxabs_n)
            if (abs(tmp(9))/B00(rind)>min_Bmn)
              %Cosinus component
              modeind=modeind+1;
              modesm{rind}(modeind)=tmp(1);
              modesn{rind}(modeind)=tmp(2);
              modesr{rind}(modeind)=tmp(3);
              modesz{rind}(modeind)=tmp(6);
              modesp{rind}(modeind)=tmp(8);
              modesb{rind}(modeind)=tmp(9);
              modesbnorm{rind}(modeind)=tmp(9)/B00(rind);
              modespar{rind}(modeind)=1; %parity 1 <=> Cosinus component
            end
            if  (abs(tmp(10))/B00(rind)>min_Bmn)
              %Sinus component
              modeind=modeind+1;
              modesm{rind}(modeind)=tmp(1);
              modesn{rind}(modeind)=tmp(2);
              modesr{rind}(modeind)=tmp(4);
              modesz{rind}(modeind)=tmp(5);
              modesp{rind}(modeind)=tmp(7);
              modesb{rind}(modeind)=tmp(10);
              modesbnorm{rind}(modeind)=tmp(10)/B00(rind);
              modespar{rind}(modeind)=0; %parity 0 <=> Sinus component
            end
          end
        end
      end
    end
    no_of_modes(rind)=modeind;
    %modesbnorm{rind}=modesb{rind}/B00(rind);
  end
  fprintf(1,'\n')

  if any(dVdsoverNper>0)
    error(['The coordinate system in the Boozer file should be left handed,'...
          ' but it has a positive Jacobian. Something is wrong!'])
  end
  

  if Geom.torfluxtot*Bphi(1)>0
    if isnan(newsigncorrectionmethod)
      error('You must specify the signcorrection method signcorr1 or signcorr2')
    end
    if not(newsigncorrectionmethod)
      Geom.newsigncorr=0;
      disp(['This is a stellarator symmetric file from Joachim Geiger.'...
            ' It will now be turned 180 degrees around a ' ...
            'horizontal axis <=> flip the sign of G and I, so that it matches the sign ' ...
            'of its total toroidal flux.'])
      Geom.headertext.maincomment=sprintf([Geom.headertext.maincomment,'\n',...
           'CC This stellarator symmetric file has been turned 180 degrees around an\n', ...
           'CC horizontal axis by flipping the signs of Itor and Ipol compared with the\n',...
           'CC original file. This was done in order for these signs to be consistent with\n',...
           'CC the total toroidal flux, which has the same sign as in the original file.']);
      
      % Such a rotation does not make it necessary to change the Fourier coefficients, 
      % in the case of Stellarator symmetry
      if not(Geom.StelSym)
        error('Rotating 180 degrees only allowed for stellarator symmetric cases!')
      end
      Bphi=-Bphi;
      Btheta=-Btheta;
    else % Use newsigncorrectionmethod
      Geom.newsigncorr=1;
      disp(['This is a stellarator symmetric file from Joachim Geiger.'...
          ' The sign of the toroidal flux will be flipped so that it it is counted positive ' ...
          ' in the direction of the toroidal coordinate. '])
      Geom.headertext.maincomment=sprintf([Geom.headertext.maincomment,'\n',...
           'CC In this stellarator symmetric file the sign of the toroidal flux\n', ...
           'CC has been flipped so that it it is counted positive in the direction\n',...
           'CC of the toroidal coordinate.']);
      Geom.torfluxtot=-Geom.torfluxtot;      
    end
  end
    
  rthetazeta_righthanded=sign(Geom.torfluxtot*Bphi(1)); %This is -1, because
                                                        %the Boozer file is supposed to be left handed
  
  if rthetazeta_righthanded==1
    error('The coordinate system in the Boozer file was right handed')
  end
  
  Geom.torfluxtot=Geom.torfluxtot*rthetazeta_righthanded;
  
  Geom.rnorm=sqrt(torfluxnorm);
  Geom.s=torfluxnorm;
  Geom.Bphi=Bphi*rthetazeta_righthanded^2;%amperes law minus sign and direction switch sign
  Geom.Btheta=Btheta*rthetazeta_righthanded;%amperes law minus sign
  Geom.iota=iota*rthetazeta_righthanded;
  Geom.dpds=dpds;
  if Geom.dpds(end)==0 %Joachim sets the last point to zero for some reason
    Geom.dpds(end)=2*Geom.dpds(end-1)-Geom.dpds(end-2); %Extrapolate!
  end
  Geom.dVdsoverNper=dVdsoverNper*rthetazeta_righthanded;
  Geom.nmodes=no_of_modes;
  Geom.m=modesm;
  Geom.n=modesn; %sign is switched below
  Geom.Bmn=modesb;
  Geom.Bnorm=modesbnorm;
  Geom.B00=B00;
  Geom.Bfilter.min_Bmn  = min_Bmn;
  Geom.Bfilter.max_m    = max_m;
  Geom.Bfilter.maxabs_n = maxabs_n;
  Geom.R00=R00;
  %Geom.Z00=zeros(size(R00));
  Geom.R=modesr;
  Geom.Z=modesz;
  Geom.Dphi=modesp;
  %for rind=1:Geom.nsurf
  %  m0n0ind=find(Geom.m{rind}==0 & Geom.n{rind}==0);
  %  Geom.Z00(rind)=Geom.Z{rind}(m0n0ind); %This is always zero (also in Erika's files)
  %end
  if rthetazeta_righthanded==-1
    for tmpind=1:length(Geom.n)
      Geom.n{tmpind}=-Geom.n{tmpind};        %This assumes the argument  
                                             %(m theta - n N phi) in (r,theta,phi) left-handed 
                                             %Correct the Boozerfile first if another
                                             %convention was used there
      Geom.Dphi{tmpind}=-Geom.Dphi{tmpind}; 
    end
  end
  if not(Geom.StelSym)
    Geom.parity=modespar;
  end
  Geom.psi_a=Geom.torfluxtot/2/pi;
  
  if not(Geom.StelSym)
    %If non-stellarator symmetric, then assume that it is from Erika.
    %Then the minor and major radii are VMEC and not Joachims definitions
    Geom.minorradiusVMEC=Geom.minorradiusW7AS;
    Geom.majorradiusVMEC=Geom.majorradiusLastbcR00;
    Geom.VolumeVMEC=pi*Geom.minorradiusVMEC^2*2*pi*Geom.majorradiusVMEC;
    Geom.minorradiusW7AS=NaN;
    Geom.majorradiusLastbcR00=NaN;
  end
  
  %fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b')
  fprintf(1,'\n')
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Henning Maassberg type Boozer file  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif strcmp(filetype,'HM')
  Geom.StelSym=1;
  
  tmp_str=fgetl(fid);
  while tmp_str(2)=='c' || tmp_str(2)=='C'; %Skip comment line
    tmp_str=fgetl(fid);
  end

  eof=0;
  rind=0;
  
  while not(eof)
    rind=rind+1;
    surfheader1=sscanf(tmp_str,'%f %f %d %f %f %f %f %f',8);
    radii(rind) = surfheader1(1)/100; %(cm)
    iota(rind)  = surfheader1(2);
    Nperiods=surfheader1(3);
    minorradiusW7AS = surfheader1(4)/100; %(cm)
    majorradiusLastbcR00  = surfheader1(5)/100; %(cm)
    if length(surfheader1)>5
      torfluxnorm(rind) = surfheader1(6);
      R00(rind)   = surfheader1(8)/100;
    else
      torfluxnorm(rind) = NaN;
      R00(rind) =NaN;
    end
    
    tmp_str=fgetl(fid);  %next line
    if tmp_str(2)=='>'
      surfheader2=sscanf(tmp_str(3:end),'%f %f %f',3);
      
      Bphi(rind)  =surfheader2(1)*1e6*Nperiods/2/pi*(4*pi*1e-7); %Tesla*meter
      Btheta(rind)=surfheader2(2)*1e6/2/pi*(4*pi*1e-7);          %Tesla*meter
      B00(rind)   =surfheader2(3); %Tesla
      tmp_str=fgetl(fid);
    else
      Bphi(rind)  =NaN;
      Btheta(rind)=NaN;
      B00(rind)   =NaN;
    end
      
    proceed=1;
    modeind=0;
    while proceed
      if sscanf(tmp_str,'%d',1)==-1 %Next flux surface has been reached
        proceed=0;
      else
        tmp=sscanf(tmp_str,'%d %d %f %f %f %f %f',7);
        if (abs(tmp(3))>min_Bmn)&&(tmp(1)<=max_m)&&(abs(tmp(2))<=maxabs_n)
          modeind=modeind+1;
          modesm{rind}(modeind)=tmp(1);
          modesn{rind}(modeind)=tmp(2);
          modesbnorm{rind}(modeind)=tmp(3);
          modesb{rind}(modeind)=tmp(4);
          modesdbdrnorm{rind}(modeind)=tmp(5)*100; %convert cm^-1 -> m^-1
          modesr{rind}(modeind)=tmp(6);
          modesz{rind}(modeind)=tmp(7);
        end
      end
      tmp_str=fgetl(fid); %get the next line or surface header line
    end
    no_of_modes(rind)=modeind;
    
    if length(tmp_str)==1
      if tmp_str==-1 %End of file has been reached
          eof=1;
      end
    end      
  end

  %Geom.torfluxtot is not stored in Henning Maassberg's files. Because they are 
  %based on Joachim Geiger's .bc files, however, they are left-handed (r,pol,tor) and
  %inconsistent with the sign of torfluxtot in the .bc file. The same correction is
  %therefore made here as for the .bc files above.
  rthetazeta_righthanded=-1;
  
  Geom.nsurf=length(radii);
  Geom.Nperiods=Nperiods;
  Geom.minorradiusW7AS=minorradiusW7AS;
  Geom.majorradiusLastbcR00=majorradiusLastbcR00;
  Geom.rnorm=radii/Geom.minorradiusW7AS;
  Geom.s=torfluxnorm;
  Geom.iota=iota*rthetazeta_righthanded;
  Geom.Bphi=Bphi;
  Geom.Btheta=Btheta*rthetazeta_righthanded;
  Geom.nmodes=no_of_modes;
  Geom.m=modesm;
  Geom.n=modesn; %sign is switched below
  Geom.Bmn=modesb;
  Geom.Bnorm=modesbnorm;
  Geom.B00=B00;
  Geom.Bfilter.min_Bmn  = min_Bmn;
  Geom.Bfilter.max_m    = max_m;
  Geom.Bfilter.maxabs_n = maxabs_n;
  Geom.R00=R00;
  Geom.R=modesr;
  Geom.Z=modesz;
  if rthetazeta_righthanded==-1
    for tmpind=1:length(Geom.n)
      Geom.n{tmpind}=-Geom.n{tmpind};
    end
  end
   
end

fclose(fid);
