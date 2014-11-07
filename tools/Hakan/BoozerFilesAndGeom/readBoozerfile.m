function Geom=readBoozerfile(filename,min_Bmn,max_m,maxabs_n,symmetry)

if nargin<2
  min_Bmn=0;
end
if nargin<4
  max_m=inf;
  maxabs_n=inf;
end
if nargin<5
  symmetry='unknown';
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
  header_f=fscanf(fid,'%f',3);
  fgetl(fid); %just skip the return character
  Geom.headertext.surfvars=fgetl(fid);  %Variable name line

  Geom.m0b        = header_d(1);
  Geom.n0b        = header_d(2);
  Geom.nsurf      = header_d(3);
  Geom.Nperiods   = header_d(4);
  Geom.torfluxtot = header_f(1); %No this is not what it is
  Geom.minorradius= header_f(2);
  Geom.majorradius= header_f(3);

  eof=0;
  rind=0;
  
  while not(eof)
    rind=rind+1;
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b%i/%i',rind,Geom.nsurf)
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
      Geom.StelSym=isempty(strfind(Geom.headertext.datavars, 'rmnc'))
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
  
  Geom.rnorm=sqrt(torfluxnorm);
  Geom.s=torfluxnorm;
  Geom.iota=iota;
  Geom.Bphi=Bphi;
  Geom.Btheta=Btheta;
  Geom.dpds=dpds;
  Geom.dVdsoverNper=dVdsoverNper;
  Geom.nmodes=no_of_modes;
  Geom.nmodes=no_of_modes;
  Geom.m=modesm;
  Geom.n=modesn;
  Geom.Bmn=modesb;
  Geom.Bnorm=modesbnorm;
  Geom.B00=B00;
  Geom.Bfilter.min_Bmn  = min_Bmn;
  Geom.Bfilter.max_m    = max_m;
  Geom.Bfilter.maxabs_n = maxabs_n;
  Geom.R00=R00;
  Geom.R=modesr;
  Geom.Z=modesz;
  Geom.Dphi=modesp;
  if not(Geom.StelSym)
    Geom.parity=modespar;
  end
  %fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b')
  fprintf(1,'\n')
  
  
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
    minorradius = surfheader1(4)/100; %(cm)
    majorradius = surfheader1(5)/100; %(cm)
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

  Geom.nsurf=length(radii);
  Geom.Nperiods=Nperiods;
  %Geom.torfluxtot=?
  Geom.minorradius=minorradius;
  Geom.majorradius=majorradius;
  Geom.rnorm=radii/Geom.minorradius;
  Geom.s=torfluxnorm;
  Geom.iota=iota;
  Geom.Bphi=Bphi;
  Geom.Btheta=Btheta;
  Geom.nmodes=no_of_modes;
  Geom.m=modesm;
  Geom.n=modesn;
  Geom.Bmn=modesb;
  Geom.Bnorm=modesbnorm;
  Geom.B00=B00;
  Geom.Bfilter.min_Bmn  = min_Bmn;
  Geom.Bfilter.max_m    = max_m;
  Geom.Bfilter.maxabs_n = maxabs_n;
  Geom.R00=R00;
  Geom.R=modesr;
  Geom.Z=modesz;
   
end

fclose(fid);
