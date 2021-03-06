function Geom=readBoozerfile(filename,min_Bmn,max_m,maxabs_n,symmetry,varargin)
%
% Geom=readBoozerfile(filename,min_Bmn,max_m,maxabs_n,symmetry,varargin)
%
% This reads the Boozer file of the .bc type
% Note that the .bc file is left-handed and the output struct Geom
% is right-handed.
%
% Only Fourier components with abs(Bmn)>min_Bmn, m<=max_m, |n|<=maxabs_n
% are read.
%
% The input symmetry can be 'StelSym' to double-check that the 
% input is stellarator symmetric
%
% If the .bc file was produced by JMC, one has to supply an extra argument
% choosing which type of correction of the sign of the total flux to make.
% varargin{1}='signcorr1': Bphi and Btheta will get sign changes
% varargin{1}='signcorr2': Total toroidal flux will get a sign change

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
  if strcmp(varargin{1},'signcorr1')
    newsigncorrectionmethod=0;
  elseif strcmp(varargin{1},'signcorr2')
    newsigncorrectionmethod=1;
  else
    error('Sign correction method for JG files not recognised!')
  end
end
% some checks
if nargin>=2 && isstr(min_Bmn)
  error('Something is wrong. The min_Bmn input was given as a string. It should be a number!')
end
if nargin>=3 && isstr(max_m)
  error('Something is wrong. The max_m input was given as a string. It should be a number!')
end
if nargin>=4 && isstr(maxabs_n)
  error('Something is wrong. The maxabs_n input was given as a string. It should be a number!')
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

  YTsign=1; %1 means no sign change, 
            %-1 means sign change because YT has resaved JG's file
  tmp_str=fgetl(fid);
  concat_str=tmp_str;
  if strcmp(tmp_str(1:2),'CC')
    if not(isempty(strfind(tmp_str,'CStconfig')))
      YTsign=-1;
    end
    while strcmp(tmp_str(1:2),'CC')
      tmp_str=fgetl(fid);
      if strcmp(tmp_str(1:2),'CC')
        if not(isempty(strfind(tmp_str,'CStconfig')))
          YTsign=-1;
        end
        concat_str = sprintf([concat_str,'\n',tmp_str]); %comment line
      end
    end
    Geom.headertext.maincomment=concat_str;
  else
    Geom.headertext.maincomment='CC ----------------------------------------------';
  end
  Geom.headertext.globalvars=tmp_str;
  header_d=fscanf(fid,'%d',4);
  %header_f=fscanf(fid,'%f',3);
  tmpstr=fgetl(fid);
  header_f=sscanf(tmpstr,'%f');
  %header_f=fscanf(fid,'%f',3);
  %fgetl(fid); %just skip the return character
  Geom.headertext.surfvars=fgetl(fid);  %Variable name line  
  YTstyle=0;
  if not(isempty(strfind(Geom.headertext.surfvars,'[A]')))
    %This is a file from Yuriy Turkin. Correct this line to the JG standard
    YTstyle=1; %indicates Yuriy Turkin style
    Geom.headertext.surfvars=...
        '       s         iota  curr_pol/nper    curr_tor    pprime   sqrt g(0,0)';
    Geom.headertext.surfvarunits=...
        '                            [A]            [A]   dp/ds,[Pa] (dV/ds)/nper';   
  end
  Geom.m0b        = header_d(1);
  Geom.n0b        = header_d(2);
  Geom.nsurf      = header_d(3);
  Geom.Nperiods   = header_d(4);
  Geom.psi_a=NaN; %Insert psi_a at this place in the list, but set it later.
  Geom.torfluxtot = header_f(1)*YTsign; %Note that this is not per pol. angle,
                                        %note: possible YTsign
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
    if not(YTstyle) %isempty(strfind(Geom.headertext.surfvars,'[A]'))
      Geom.headertext.surfvarunits=fgetl(fid); %unit line only in JG files
    end
    surfheader=fscanf(fid,'%f',6);
    torfluxnorm(rind) = surfheader(1);
    iota(rind)  = surfheader(2);
    Bphi(rind)  = surfheader(3)*Geom.Nperiods/2/pi*(4*pi*1e-7); %Tesla*meter
    Btheta(rind)= surfheader(4)/2/pi*(4*pi*1e-7);               %Tesla*meter
    dpds(rind)  = surfheader(5);
    dVdsoverNper(rind)  = surfheader(6);
    fgetl(fid); %just skip the return character
    tmpstrunits=fgetl(fid); %units line
    if rind==1
      Geom.headertext.datavars=tmpstrunits;
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
    tmp_str1=fgetl(fid);
    if Geom.StelSym
      tmp=sscanf(tmp_str1,'%d %d %f %f %f %f',6);
      while not(tmp(1)==0 && tmp(2)==0)
        tmp_str1=fgetl(fid);
        tmp=sscanf(tmp_str1,'%d %d %f %f %f %f',6);    
      end
      B00(rind)=tmp(6);
      R00(rind)=tmp(3);
    else
      tmp=sscanf(tmp_str1,'%d %d %f %f %f %f %f %f %f %f',10);
      while not(tmp(1)==0 && tmp(2)==0)
        tmp_str1=fgetl(fid);
        tmp=sscanf(tmp_str1,'%d %d %f %f %f %f %f %f %f %f',10);
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
    if no_of_modes(rind)==0
      error(['no modes found for rind=',num2str(rind)])
    end
  end
  fprintf(1,'\n')
  
  if any(dVdsoverNper>0)
    error(['The coordinate system in the Boozer file should be left handed,'...
           ' but it has a positive Jacobian. Something is wrong!'])
  end
  
  if Geom.torfluxtot*Bphi(1)>0
    if not(Geom.StelSym)
      error('Unknown sign convention in non-stellarator-symmetric file!') 
    end
    if isnan(newsigncorrectionmethod)
      error('You must specify the signcorrection method signcorr1 or signcorr2')
    end
    if not(newsigncorrectionmethod)
      Geom.newsigncorr=0;
      if 1 %Optional output to screen
        disp(['This was a stellarator symmetric file from Joachim Geiger.'...
              ' It has been turned 180 degrees around a ' ...
              'horizontal axis <=> flip the sign of G and I, so that it matches the sign ' ...
              'of its total toroidal flux.'])
      end
      if 0 %Not necessary to write, because writeBoozerfile.m put the signs back again
        Geom.headertext.maincomment=sprintf([Geom.headertext.maincomment,'\n',...
                            'CC This stellarator symmetric file has been turned 180 degrees around an\n', ...
                            'CC horizontal axis by flipping the signs of Itor and Ipol compared with the\n',...
                            'CC original file. This was done in order for these signs to be consistent with\n',...
                            'CC the total toroidal flux, which has the same sign as in the original file.']);
      end
      % Such a rotation does not make it necessary to change the Fourier coefficients, 
      % in the case of Stellarator symmetry
      if not(Geom.StelSym)
        error('Rotating 180 degrees only allowed for stellarator symmetric cases!')
      end
      Bphi=-Bphi;
      Btheta=-Btheta;
    else % Use newsigncorrectionmethod
      Geom.newsigncorr=1;
      if 0 %Optional output to screen
        disp(['This was a stellarator symmetric file from Joachim Geiger.'...
              ' The sign of the toroidal flux has been flipped so that it it is counted positive' ...
              ' in the direction of the toroidal coordinate irrespectively of handedness.'])
      end
      if 0 %Not necessary to write, because writeBoozerfile.m put the signs back again
        Geom.headertext.maincomment=sprintf([Geom.headertext.maincomment,'\n',...
                            'CC In this stellarator symmetric file the sign of the toroidal flux\n', ...
                            'CC has been flipped so that it it is counted positive in the direction\n',...
                            'CC of the toroidal coordinate.']);
      end
      Geom.torfluxtot=-Geom.torfluxtot;      
    end
  elseif Geom.StelSym %and Geom.torfluxtot*Bphi(1)<0
    error('Unknown sign convention in tellarator symmetric file! Neither YT nor JG standard.')  
  end
  
  rthetazeta_righthanded=sign(Geom.torfluxtot*Bphi(1)); %This is -1, because
                                                        %the Boozer file is supposed to be left handed
  
  if rthetazeta_righthanded==1
    error('The coordinate system in the Boozer file was right handed')
  end
  
  Geom.torfluxtot=Geom.torfluxtot*rthetazeta_righthanded;  
  Geom.psi_a=Geom.torfluxtot/2/pi;
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
  Geom.FSAB2=abs(4*pi^2*Geom.psi_a/Geom.Nperiods*...
                 (Geom.Bphi+Geom.iota.*Geom.Btheta)./Geom.dVdsoverNper);
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
  %Dphi =N/(2*pi)*(zeta-(-geomang))
  
  if not(Geom.StelSym)
    Geom.parity=modespar;
  end
  
  if not(Geom.StelSym)
    %If non-stellarator symmetric, then assume that it is from Erika.
    %Then the minor and major radii are VMEC and not Joachims definitions
    Geom.minorradiusVMEC=Geom.minorradiusW7AS;
    Geom.majorradiusVMEC=Geom.majorradiusLastbcR00;
    Geom.VolumeVMEC=pi*Geom.minorradiusVMEC^2*2*pi*Geom.majorradiusVMEC;
    Geom.minorradiusW7AS=NaN;
    Geom.majorradiusLastbcR00=NaN;
  end
  
  %%%% Some security checks 
  Ntheta=301;Nzeta=303;rind=1;
  Bmn=mnmat(Geom,rind,'B',Ntheta,Nzeta,'forceSize');
  B=ifftmn(Bmn);
  FSAB2_test=Ntheta*Nzeta/sum(sum(B.^(-2)));
  if abs(FSAB2_test./Geom.FSAB2(rind)-1)>0.05
    warning('dVdsoverNper not correct in bc file')
    Geom.dVdsoverNper=NaN*Geom.dVdsoverNper;
    Geom.FSAB2=NaN*Geom.FSAB2;
    abs(4*pi^2*Geom.psi_a/Geom.Nperiods*(Geom.Bphi+Geom.iota.*Geom.Btheta)./Geom.dVdsoverNper);
    %Calculate the correct values
    for rind=1:Geom.nsurf
      Bmn=mnmat(Geom,rind,'B',Ntheta,Nzeta,'forceSize');
      B=ifftmn(Bmn);
      Geom.FSAB2(rind)=Ntheta*Nzeta/sum(sum(B.^(-2)));
    end
    Geom.dVdsoverNper=abs(4*pi^2*Geom.psi_a/Geom.Nperiods*(Geom.Bphi+Geom.iota.*Geom.Btheta)./Geom.FSAB2);
  end
  
  
  %fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b')
  fprintf(1,'\n')
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Henning Maassberg type Boozer file  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif strcmp(filetype,'HM')
  Geom.StelSym=1;
  
  if isnan(newsigncorrectionmethod)
    newsigncorrectionmethod=0; 
    Geom.newsigncorr=0;
  else
    Geom.newsigncorr=newsigncorrectionmethod;
  end
  
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
   
  if not(Geom.newsigncorr)
    Geom.Bphi=-Geom.Bphi;
    Geom.Btheta=-Geom.Btheta;
  end
end

fclose(fid);
