function make_monoenergeticrunspecfile_w7x(pathtoall)


%THIS FILE GENERATES A runspec.dat FILE
% Either give the .bc file as input
% or the path to where the gidkes.cnf file resides.

outputonscreenonly=0;

%
% This produces the file runspec.dat which is used by the ruby script 
% and programMode = 21 in sfincs to launch mono-energetic runs
%
bcfile_was_given=0;
if length(pathtoall)>=3
  if pathtoall(end-2:end)=='.bc'
    bcfile_was_given=1;
    B_filename=pathtoall;
    pathtoall='';
  end
end
if not(strcmp(pathtoall,''))
  if pathtoall(end)~='/'
    pathtoall=[pathtoall,'/'];
  end
end


%Before running this file, add the paths to
addpath([getenv('SFINCS_PROJECTSANDTOOLS_HOME'),'/tools/Hakan/BoozerFilesAndGeom']) %This is where readBoozerfile.m resides.
addpath([getenv('SFINCS_PROJECTSANDTOOLS_HOME'),'/tools/Hakan/version3Scans']) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hennings experience with the configuration
% w7x-sc1-ecb2.bc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HM.cmul=[10,3,1,0.3,0.1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5,...
         3e-6,1e-6,3e-7];
HM.efield_norm=[1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3,2e-3,5e-3,...
                0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.8,1.0,...
                1.2,1.5,2.0,3.0,5.0];
%Below this EFIELD the results are very close to the E=0 results:
HM.efield_min=[1e-5,3e-5,1e-4,3e-4,1e-3,1e-3,1e-3,3e-4,3e-4,1e-4,...
               3e-5,1e-5,3e-6,1e-6,3e-7,1e-7];


if bcfile_was_given
  useHMnormEr=1;
  useHMminEr=1;
  useHMcmul=1;
else
  useHMnormEr_str=...
      input(['Use Henning''s Er values',...
             ' (return=yes, n=no) ?'],'s');
  if isempty(useHMnormEr_str)
    useHMnormEr=1;
  elseif strcmp(useHMnormEr_str,'y')
    useHMnormEr=1;
  else
    useHMnormEr=0;
  end
  useHMnormEr=isempty(useHMnormEr_str);

  useHMminEr_str=...
      input(['Use Henning''s minimum Er depending on cmul',...
             ' (return=yes, n=no) ?'],'s');
  if isempty(useHMminEr_str)
    useHMminEr=1;
  elseif strcmp(useHMminEr_str,'y')
    useHMminEr=1;
  else
    useHMminEr=0;
  end
  useHMminEr=isempty(useHMminEr_str);

  useHMcmul_str=...
      input(['Use Henning''s cmul values',...
             ' (return=yes, n=no=take from gidkes.cnf) ?'],'s');
  if isempty(useHMcmul_str)
    useHMcmul=1;
  elseif strcmp(useHMcmul_str,'y')
    useHMcmul=1;
  else
    useHMcmul=0;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sfincs experience with the configuration
% w7x-sc1-ecb2.bc
% note that transportCoeffs(2,1) has converged but not 
% transportCoeffs(1,2) at high collisionality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%r/a=0.5, Er=0:
SF.NuPrimes     =  [3e-5,1e-4,3e-4,1e-3,1e-2,0.1,0.3,  1, 10,100];
SF.Nthetas      =  [  41,  41,  35,  41,  47, 47, 47, 45, 37, 37];
SF.Nzetas       =  [ 260, 200, 150,  90,  47, 37, 37, 37, 37, 37];
SF.Nxis         =  [ 220, 180, 150,  90,  37, 21, 15,  9,  7,  5];
SF.hydra.Nnodes =  [   16, 16,   8,   1,   1,  1,  1,  1,  1,  1];
SF.hydra.Nppernode=[    4,  4,   4,   4,   2,  2,  1,  1,  1,  1];
SF.hydra.ConsCpus =[    4,  4,   4,   4,   4,  4,  4,  4,  4,  4];
SF.hydra.walClock =[   90, 60,  60,  60,  30, 30, 30, 30, 30, 30];
%SF.hgw.Nproc    =  [ NaN, NaN, NaN, NaN, NaN,NaN,NaN,NaN,NaN,NaN];

%for nuPrime=1e-4 16 nodes goes under 30 min, and 8 nodes also work but take ~60 min

%r/a=0.9, Er=0:
%SF.NuPrimes     =  [3e-5,1e-4,3e-4,1e-3,1e-2,0.1,0.3,  1, 10,100];
%SF.Nthetas      =  [                      47, 47, 41, 41, 41, 41];
%SF.Nzetas       =  [                      47, 41, 41, 41, 41, 41];
%SF.Nxis         =  [                      37, 21, 15,  9,  5,  5];
%SF.hydra.Nnodes =  [   16, 16,   8,   1,   1,  1,  1,  1,  1,  1];
%SF.hydra.Nppernode=[    4,  4,   4,   4,   2,  2,  1,  1,  1,  1];
%SF.hydra.ConsCpus =[    4,  4,   4,   4,   4,  4,  4,  4,  4,  4];
%SF.hydra.walClock =[   90, 60,  60,  60,  30, 30, 30, 30, 30, 30];

%r/a=0.1, Er=0:
%SF.NuPrimes     =  [3e-5,1e-4,3e-4,1e-3,1e-2,0.1,0.3,  1, 10,100];
%SF.Nthetas      =  [                                            ];
%SF.Nzetas       =  [                                            ];
%SF.Nxis         =  [                                            ];
%SF.hydra.Nnodes =  [   16, 16,   8,   ?,   1,  1,  1,  1,  1,  1];
%SF.hydra.Nppernode=[    4,  4,   4,   4,   2,  2,  1,  1,  1,  1];
%SF.hydra.ConsCpus =[    4,  4,   4,   4,   4,  4,  4,  4,  4,  4];
%SF.hydra.walClock =[   90, 60,  60,  60,  30, 30, 30, 30, 30, 30];


%SF.solverTolerance = 1e-6;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read input from gidkes.cnf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bcfile_was_given
  min_Bmn = 1e-6;                        % minimum Bnorm to include
  max_m   = 10;                          % maximum m to include
  maxabs_n= 20;                          % maximum |n| to include
else %read gidkes.cnf in the dir pathtoall
     %% Some default values
     %b0     = 1.0;     % magnetic field strength in T (not read from gidkes.cnf)
     %eps_l  = 1.e-4;   % minimum B_mn to include  (read from gidkes.cnf)
     %mp_u   = 10;      % maximum m to include     (read from gidkes.cnf)
     %nt_u   = 20;      % maximum |n| to include   (read from gidkes.cnf)

  gidkes_fid=fopen([pathtoall,'gidkes.cnf'],'r');

  ncmul  = 0;       % no. of CMUL values (collision frequencies to loop over)
  nefield= 0;       % no. of EFIELD values (E fields to loop over)

  str=fgetl(gidkes_fid); 
  strs= textscan(str,'%s');          % filename with Fourier coefficients for |B|
  B_filename=strs{1}{1};

  disp(['Input taken from the file ',B_filename]);

  b0dummy    = fscanf(gidkes_fid,'%f',1); % magnetic field strength in T
                                          %->not used, always use 1.0
  fgetl(gidkes_fid);                      % skip to next line

  val        = fscanf(gidkes_fid,'%f, %d, %d'); % read three values
  fgetl(gidkes_fid);                      % skip to next line
  min_Bmn = val(1);                          % minimum Bnorm to include
  max_m   = val(2);                          % maximum m to include
  maxabs_n= val(3);                          % maximum |n| to include

  ncmul_fromcnf = fscanf(gidkes_fid,'%d',1); % no. of CMUL values
  fgetl(gidkes_fid);                      % skip to next line
  str=fgetl(gidkes_fid);
  format='%g';
  for i=1:ncmul_fromcnf-1
    format=[format,',%g'];
  end
  cmul_fromcnf = sscanf(str,format,ncmul); % read the cmul values (collision freqs to loop over)
  

  nefield_fromcnf = fscanf(gidkes_fid,'%d',1);  % no. of CMUL values
  fgetl(gidkes_fid);                       % skip to next line
  str=fgetl(gidkes_fid);
  format='%g';
  for i=1:nefield_fromcnf-1
    format=[format,',%g'];
  end
  efield_fromcnf = sscanf(str,format,ncmul);  % read the efield values (E fields to loop over)

  fclose(gidkes_fid);                     % end of the input from the
                                          % gidkes file      
end

if useHMcmul
  cmul = HM.cmul;
  ncmul=length(cmul);
else
  cmul  = cmul_fromcnf;
  ncmul = ncmul_fromcnf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read the Boozer files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bcfilename=B_filename;
if B_filename(end-1:end)=='bc'
  B_filetype='JG';
  strconf=B_filename(1:end-3);
else
  B_filetype='HM';
  if B_filename(end-3:end)~='.dat'
    bcfilename=[B_filename,'.dat'];
    strconf=B_filename;
  else
    strconf=B_filename(1:end-4);
  end  
end
disp([' Magnetic field Fourier coefficients taken from ',bcfilename])  
fprintf(1,'\n')

Geom=readBoozerfile([pathtoall,bcfilename],min_Bmn,max_m,maxabs_n,'StelSym');
Nr=length(Geom.rnorm);
disp(['Number of radii = ',num2str(Nr)])
if Nr>20
  error('Are you sure that you want that many radii?')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Open the runspec file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
runspecfile='runspec.dat';
if outputonscreenonly
  runspec_fid=1;
else
  runspec_fid=fopen(runspecfile,'w');
end

%HEADER
fprintf(runspec_fid,'! This file was generated by make_runspecfile_w7x.m\n');
fprintf(runspec_fid,'!\n');
fprintf(runspec_fid,['!        rN_wish',...
                    '     nuPrime',...
                    '         EStar',...
                    ' Ntheta',...
                    ' Nzeta',...
                    '   Nxi',...
                    ' Nnodes',...
                    ' Nppernode',...
                    ' ConsCpus',...
                    ' WalClock\n']);
%                    ' Nx',...
%                    '         Delta',...
%                    '         omega',...
%                    ' collisionOperator',...
%                    ' programMode',...
%                    '     RHSMode\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Set some constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delta = 4.5694e-3;
omega = 2.2847e-3;
collisionOperator=1;
programMode=1;
RHSMode=3;
Nx=1;
Chandra1=(erf(sqrt(1))-sqrt(1)*2/sqrt(pi).*exp(-1))./2./1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% LOOP over radii %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for rind=1:Nr 
  normradius_wish=Geom.rnorm(rind);
  B00 = Geom.B00(rind);
  iota= Geom.iota(rind);
  G   = Geom.Bphi(rind);
  I   = Geom.Btheta(rind);
  %dPsidr=2*Geom.torfluxtot/2/pi/Geom.minorradius*Geom.rnorm(rind);
  dPsidr_GIDKES=Geom.minorradius*Geom.rnorm(rind)*B00;
  nuPrime_over_cmul=4/3/sqrt(pi)*(G+iota*I)/B00/(erf(1)-Chandra1);
  
  disp(['--------- Radius ',num2str(rind),' of ',...
        num2str(Nr),',  nuPrime/cmul=', ...
        num2str(nuPrime_over_cmul),'-----------'])
  for cmulind=1:ncmul
    
    nuPrime=nuPrime_over_cmul*cmul(cmulind);
    
    if nuPrime < min(abs(SF.NuPrimes))
      disp(['nuPrime=',num2str(nuPrime),...
            ' is less than SFINCS can run and is ignored.'])
    else
      
      Ntheta=ceil(interp1(log(SF.NuPrimes),SF.Nthetas,log(abs(nuPrime))));
      Nzeta =ceil(interp1(log(SF.NuPrimes),SF.Nzetas, log(abs(nuPrime))));
      Nxi   =ceil(interp1(log(SF.NuPrimes),SF.Nxis,   log(abs(nuPrime))));
      nuPvind=floor(interp1(log(SF.NuPrimes),1:length(SF.NuPrimes),log(abs(nuPrime))));
      Nnodes=SF.hydra.Nnodes(nuPvind);
      Nppernode=SF.hydra.Nppernode(nuPvind);
      ConsCpus=SF.hydra.ConsCpus(nuPvind);
      WalClock=SF.hydra.WalClock(nuPvind);

      if useHMnormEr
        efield=HM.efield_norm'*B00*iota*...
            Geom.rnorm(rind)*Geom.minorradius/Geom.majorradius;
      else
        efield=efield_fromcnf;
      end

      
      efield_min=exp(interp1(log(HM.cmul),log(HM.efield_min),...
                             log(cmul(cmulind)),'pchip'));
      if useHMminEr
        efield_forthis_cmul=[0;efield(find(abs(efield)>=efield_min))];
        if length(efield_forthis_cmul)==1
          efield_forthis_cmul=...
              [0;efield(find(abs(efield)==max(abs(efield))))];
        end
      else
        efield_min=NaN;
        efield_forthis_cmul=[0;efield];
      end
      nefld=length(efield_forthis_cmul);
      
      for efldind=1:nefld
        
        EStar = -G/iota/dPsidr_GIDKES*efield_forthis_cmul(efldind);
        %EStar = -G/iota/dPsidr*efield_forthis_cmul(efldind)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Now write to the runspecfile
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %RUNSPEC DATA

        fprintf(runspec_fid,...
                ['%14.6e%14.6e%14.6e',...
                 '%7i%6i%6i',...
                 '%7i%10i%9i%9i\n'],...
                normradius_wish,nuPrime,EStar,...
                Ntheta,Nzeta,Nxi,...
                Nnodes,Nppernode,ConsCpus,WalClock);
      end
    end
  end



end
if not(outputonscreenonly)
  fclose(runspec_fid);
end



