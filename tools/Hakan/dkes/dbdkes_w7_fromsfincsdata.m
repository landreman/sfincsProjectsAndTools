function dbdkes_w7_fromsfincsdata(runsdir)
addpath([getenv('SFINCS_PROJECTSANDTOOLS_HOME'),'/tools/Hakan/BoozerFilesAndGeom'])
addpath([getenv('SFINCS_PROJECTSANDTOOLS_HOME'),'/tools/Hakan/dkes'])
addpath([getenv('SFINCS_PROJECTSANDTOOLS_HOME'),'/tools/Hakan/version3Scans'])

[runs,missing]=getresults(runsdir);

% determine the path to the .bc file
input_namelist=runs.run.input_namelist;
indeqF=findstr(input_namelist,'equilibriumFile =');
strstart=indeqF+17-1+1+findstr(input_namelist(indeqF+17:indeqF+25),'"');
strend=strstart-1+findstr(input_namelist(strstart:strstart+50),'"');
bcfile=input_namelist(strstart:strend-1);
indslash=findstr(bcfile,'/');
if isempty(indslash)
  confnamestart=1;
else
  confnamestart=indslash(end)+1;
end
indpbc=findstr(bcfile,'.bc');
if isempty(indpbc)
  confnameend=length(bcfile);
else
  confnameend=indpbc-1;
end
strconf=bcfile(confnamestart:confnameend);

disp(['Loading the boozer file ',bcfile])
Geom=readBoozerfile(bcfile);

%bcfile=[CONFBC,strconf,'.bc'];
%datfile=[CONF,strconf,'.dat'];

%%%%%%% Check if the runs finished %%%%%%%%%%%%%%
if not(isempty(missing))
  answer=input('All DKES runs are not completed! Continue? (return or n): ')
  if not(isempty(answer))
    error('Stopped by user!')
  end
end

%a mathematical constant:
Chandra1=(erf(sqrt(1))-sqrt(1)*2/sqrt(pi).*exp(-1))./2./1;

%%%%%%%%%%%%%%%% sort the runs and calculate cmul and efield %%%%%%%%%%%%%%%%%%%%%%%%%

[rN1,isort1]=sort(runs.rN);
nuPrime1=runs.nuPrime(isort1);
EStar1=runs.EStar(isort1);
transportCoeffs1=runs.transportCoeffs(isort1,:,:);
rindstarts=[1,find(diff(rN1))+1];
dkdata.Nruns=diff([rindstarts,length(rN1)+1]);
dkdata.Nradii=length(rindstarts);
dkdata.rN=rN1(rindstarts);

for rind=1:length(dkdata.Nruns)
  rN=rN1(rindstarts(rind));
  [~,Grind(rind)]=min((Geom.rnorm-rN).^2);
  G=Geom.Bphi(Grind(rind));
  I=Geom.Btheta(Grind(rind));
  iota=Geom.iota(Grind(rind));
  B00=Geom.B00(Grind(rind));
  dPsidr_GIDKES=rN*Geom.minorradius*B00;
  
  nuPrimes=nuPrime1(rindstarts(rind):rindstarts(rind)+dkdata.Nruns(rind)-1);
  EStars=EStar1(rindstarts(rind):rindstarts(rind)+dkdata.Nruns(rind)-1);
  transportCoeffss=transportCoeffs1(rindstarts(rind):rindstarts(rind)+dkdata.Nruns(rind)-1,:,:);
  
  dkdata.cmul{rind}=nuPrimes*3*sqrt(pi)/4*(erf(1)-Chandra1)*B00/(G+iota*I);
  dkdata.efld{rind}=-EStars*iota/G*dPsidr_GIDKES;
  dkdata.d11{rind}=squeeze(transportCoeffss(:,1,1))' ...
        /dPsidr_GIDKES^2*sqrt(pi)/8*G^2/(G+iota*I)/B00;
  dkdata.d31{rind}=squeeze(transportCoeffss(:,2,1)+transportCoeffss(:,1,2))'/2 ...
        /dPsidr_GIDKES*sqrt(pi)/4*G;
  dkdata.d33{rind}=squeeze(transportCoeffss(:,2,2))' ...
      *(-sqrt(pi)/2)*(G+iota*I)*B00;
  dkdata.d11e{rind}=zeros(size(dkdata.d11{rind}));
  dkdata.d31e{rind}=zeros(size(dkdata.d31{rind}));
  dkdata.d33e{rind}=zeros(size(dkdata.d33{rind}));
  dkdata.fmn{rind}=NaN*zeros(size(dkdata.d11{rind}));
  dkdata.lalpha{rind}=NaN*zeros(size(dkdata.d11{rind}));
  
  for runind=1:dkdata.Nruns(rind)
    %example:    
    %  1.000E+01   0.000E+00  -1.49453E+01  -4.52591E-06  -6.65745E-02
    % >3            49    60   3.66407E-10   9.92864E-16   7.74470E-17
    dkdata.line1{rind}{runind}=sprintf('%11.3E %11.3E %13.5E %13.5E %13.5E',...
                                       dkdata.cmul{rind}(runind),...
                                       dkdata.efld{rind}(runind),...
                                       dkdata.d11{rind}(runind),...
                                       dkdata.d31{rind}(runind),...
                                       dkdata.d33{rind}(runind));   
    dkdata.line2{rind}{runind}=sprintf(' >3       %6d %6d %13.5E %13.5E %13.5E',...
                                       dkdata.cmul{rind}(runind),...
                                       dkdata.efld{rind}(runind),...
                                       dkdata.d11{rind}(runind),...
                                       dkdata.d31{rind}(runind),...
                                       dkdata.d33{rind}(runind));
  end
  %[cmul2,isort2]=sort(dkdata.cmul{rind});
  %efields2=dkdata.efield{rind}(isort2);
  %transportCoeffs2=transportCoeffss(isort2,:,:);
  % 
  %cmulstarts=[1,find(diff(cmul2))+1];
  %Ncmuls=diff([cmulstarts,length(cmul2)+1]);
  %
  %for cmulind=1:length(Ncmuls)
  %  cmul=cmul2(cmulstarts(cmulind));
  %  data.cmul{rind}=cmul;
  %  
  %  efields_unsrt=efields2(cmulstarts(cmulind):cmulstarts(cmulind)+Ncmuls(cmulind)-1);
  %  transportCoeffs_unsrt=transportCoeffs2(cmulstarts(cmulind):cmulstarts(cmulind)+ ...
  %                                         Ncmuls(cmulind)-1,:,:);
  %  [efield3,isort3]=sort(efields_unsrt);
  %  transportCoeffs3=transportCoeffs_unsrt(isort3);
  %  
  %  dkdata.cmulfirst.efld{rind}{cmulind}=efield3;
  %  %data.transportCoeffs{rind}{cmulind}=transportCoeffs3;
  %  dkdata.cmulfirst.d11{rind}{cmulind}=transportCoeffs3(:,1,1)...
  %      /dPsidr_GIDKES^2*sqrt(pi)/8*G^2/(G+iota*I)/B00;
  %  dkdata.cmulfirst.d31{rind}{cmulind}=(transportCoeffs3(:,2,1)+transportCoeffs3(:,1,2))/2...
  %      /dPsidr_GIDKES*sqrt(pi)/4*G;
  %  dkdata.cmulfirst.d33{rind}{cmulind}=transportCoeffs3(:,2,2)...
  %      *(-sqrt(pi)/2)*(G+iota*I)*B00;
  %end
end



% %%%%%%%%%%%%%%%% Load the runs (DKES version) %%%%%%%%%%%%%%%%%%%%%%%%%
% oldradiusind=1;
% runindforthisr=0;
% for runind=1:length(runs.radius_str)
%   rundk_fid=fopen([dkesoutputpath,runs.dkfilename{runind}]);
%   radiusind=str2num(runs.radius_str{runind});
%   if oldradiusind~=radiusind
%     dkdata.Nruns(oldradiusind)=runindforthisr;
%     runindforthisr=0;
%   end
%   oldradiusind=radiusind;
%   for entryind=1:runs.Nentries(runind)
%     runindforthisr=runindforthisr+1;
%     dkdata.line1{radiusind}{runindforthisr}=fgetl(rundk_fid);
%     dkdata.line2{radiusind}{runindforthisr}=fgetl(rundk_fid);
%     use_this_entry=1;
%     if (dkdata.line1{radiusind}{runindforthisr}(2)=='c' ||...
%         dkdata.line1{radiusind}{runindforthisr}(2)=='C')
%       if deactivate_data_with_error_warning
%         runindforthisr=runindforthisr-1; %disregard commented data
%         use_this_entry=0;
%       else
%         dkdata.line1{radiusind}{runindforthisr}=[' ',...
%                             dkdata.line1{radiusind}{runindforthisr}(3:end)];
%         dkdata.line2{radiusind}{runindforthisr}=[' ',...
%                             dkdata.line2{radiusind}{runindforthisr}(3:end)];
%       end
%     end
%     if use_this_entry
%       line1vals=sscanf(dkdata.line1{radiusind}{runindforthisr},'%E',5);
%       line2vals=sscanf(dkdata.line2{radiusind}{runindforthisr}(5:end),...
%                        '%d%d%E%E%E',5);
%       isadouble=0;
%       if runindforthisr>1
%         doubleind=find(line1vals(1)==dkdata.cmul{radiusind}(1:runindforthisr-1) & ...
%                        line1vals(2)==dkdata.efld{radiusind}(1:runindforthisr-1));
%         if not(isempty(doubleind))
%           isadouble=1;
%           runindforthisr=runindforthisr-1; %disregard double data
%           fprintf(1,['Found double data (loading only one):\n',...
%                      '[%8.4E,%8.4E,%8.4E,%8.4E,%8.4E\n',...
%                      ' %8.4E,%8.4E,%8.4E,%8.4E,%8.4E]\n'],...
%                   line1vals(1),line1vals(2),line1vals(3),line1vals(4),line1vals(5),...
%                   dkdata.cmul{radiusind}(doubleind(1)),...
%                   dkdata.efld{radiusind}(doubleind(1)),...
%                   dkdata.d11{radiusind}(doubleind(1)),...
%                   dkdata.d13{radiusind}(doubleind(1)),...
%                   dkdata.d33{radiusind}(doubleind(1)))                  
%         end
%       end
%       if not(isadouble)
%         dkdata.cmul{radiusind}(runindforthisr)=line1vals(1);
%         dkdata.efld{radiusind}(runindforthisr)=line1vals(2);
%         dkdata.d11{radiusind}(runindforthisr) =line1vals(3);
%         dkdata.d13{radiusind}(runindforthisr) =line1vals(4);
%         dkdata.d33{radiusind}(runindforthisr) =line1vals(5);
%         dkdata.fmn{radiusind}(runindforthisr) =line2vals(1);
%         dkdata.lalpha{radiusind}(runindforthisr) =line2vals(2);
%         dkdata.d11e{radiusind}(runindforthisr) =line2vals(3);
%         dkdata.d13e{radiusind}(runindforthisr) =line2vals(4);
%         dkdata.d33e{radiusind}(runindforthisr) =line2vals(5);
%       end
%     end
%   end
%   fclose(rundk_fid);
% end
% dkdata.Nruns(radiusind)=runindforthisr; %Nruns for the last radius
% dkdata.Nradii=str2num(runs.radius_str{end});

for rind=1:dkdata.Nradii
  %structure the data according to {cmulind}(efldind)
  %begin by resorting (S)
  [cmulS,indS]=sort(dkdata.cmul{rind},2,'descend');
  startinds=[0,find(diff(cmulS))]+1;
  dkdata.cmulfirst{rind}.Neflds=diff([startinds,length(dkdata.cmul{rind})+1]);
  dkdata.cmulfirst{rind}.Ncmuls=length(dkdata.cmulfirst{rind}.Neflds);
  dkdata.cmulfirst{rind}.cmul=dkdata.cmul{rind}(indS(startinds));
  for cmulind=1:dkdata.cmulfirst{rind}.Ncmuls
    unsortinds=find(dkdata.cmul{rind}==dkdata.cmulfirst{rind}.cmul(cmulind));
    unsort_efld=dkdata.efld{rind}(unsortinds);
    [~,sortsubinds]=sort(abs(unsort_efld));
    sortinds=unsortinds(sortsubinds);
    dkdata.cmulfirst{rind}.cmulinds{cmulind}=sortinds;
    dkdata.cmulfirst{rind}.efld{cmulind}=dkdata.efld{rind}(sortinds);
    dkdata.cmulfirst{rind}.d11{cmulind}=dkdata.d11{rind}(sortinds);
    dkdata.cmulfirst{rind}.d11e{cmulind}=dkdata.d11e{rind}(sortinds);
  end
  
  
  %structure the data according to {efldind}(cmulind)
  %begin by resorting (S)
  [efldS,indS]=sort(dkdata.efld{rind});
  startinds=[0,find(diff(efldS))]+1;
  dkdata.efldfirst{rind}.Ncmuls=diff([startinds,length(efldS)+1]);
  dkdata.efldfirst{rind}.Neflds=length(dkdata.efldfirst{rind}.Ncmuls);
  dkdata.efldfirst{rind}.efld=dkdata.efld{rind}(indS(startinds));
  for efldind=1:dkdata.efldfirst{rind}.Neflds
    unsortinds=find(dkdata.efld{rind}==dkdata.efldfirst{rind}.efld(efldind));
    unsort_cmul=dkdata.cmul{rind}(unsortinds);
    [~,sortsubinds]=sort(abs(unsort_cmul),2,'descend');
    sortinds=unsortinds(sortsubinds);
    dkdata.efldfirst{rind}.efldinds{efldind}=sortinds;
    dkdata.efldfirst{rind}.cmul{efldind}=...
        dkdata.cmul{rind}(sortinds);
    dkdata.efldfirst{rind}.d11{efldind}=...
        dkdata.d11{rind}(sortinds);
    dkdata.efldfirst{rind}.d11e{efldind}=...
        dkdata.d11e{rind}(sortinds);
    if not(isempty(find(diff(dkdata.efldfirst{rind}.cmul{efldind})>0)))
      error('cmuls are in wrong order!')
    end
  end
end


%%%%%%%%%%%% Interactive fit to the D_11 data %%%%%%%%%%%%%%%%%%%

%%%%Geomdat=readBoozerfile(datfile);
%%%%% Make reduced data struct with only the chosen radii:
Geomdat.rnorm=Geom.rnorm(Grind);
Geomdat.B00=Geom.B00(Grind);
Geomdat.R00=Geom.R00(Grind);
Geomdat.iota=Geom.iota(Grind);
Geomdat.Bphi=Geom.Bphi(Grind);
Geomdat.Btheta=Geom.Btheta(Grind);
Geomdat.majorradius=Geom.majorradius;
Geomdat.minorradius=Geom.minorradius;
Geomdat.m={Geom.m{Grind}};
Geomdat.n={Geom.n{Grind}};
Geomdat.Bnorm={Geom.Bnorm{Grind}};
Geomdat.Bmn={Geom.Bmn{Grind}};

dkdata

d11fits=dkes_ftd11(Geomdat,dkdata);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Begin writing to the .dk file %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' THE .DK FILE IS BEING WRITTEN')

dk_fid=fopen([strconf,'.dk'],'w');
fprintf(dk_fid,'cc   %s:  configuration nickname for DKES database file\n',strconf);
fprintf(dk_fid,'cc   Generated by program dbdkes_w7_fromsfincsdata.m at %s\n',datestr(now));

% Copy some comments from the .bc file to the .dk file
bcfid=fopen(bcfile);
if bcfid~=-1
  bcfile_exists=1;
  bc_line=fgetl(bcfid);
  Geiger_date='';
  while lower(bc_line(1:2))=='cc'
    if lower(bc_line(4:12))=='created: '
      Geiger_date=bc_line(4:end);
    elseif lower(bc_line(4:12)=='based on ')
      fprintf(dk_fid,'cc   VMEC %s\n',bc_line(13:end));
      fprintf(dk_fid,'cc   %s by Joachim Geiger\n',Geiger_date);
    elseif lower(bc_line(4:9))=='<beta>'
      fprintf(dk_fid,'cc   %s\n',bc_line(4:end));
    end
    bc_line=fgetl(bcfid);
  end 
  fclose(bcfid);
else
  bcfile_exists=0;
end


%%%%%%%%%%%% Loop over radii %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for rind=1:dkdata.Nradii
  %ind11=find(Geomdat.m{rind}==1 & Geomdat.n{rind}==1);
  %ind01=find(Geomdat.m{rind}==0 & Geomdat.n{rind}==1);
  ind10=find(Geomdat.m{rind}==1 & Geomdat.n{rind}==0);
  %b_11=Geomdat.Bnorm{rind}(ind11);
  %b_01=Geomdat.Bnorm{rind}(ind01);
  b_10=Geomdat.Bnorm{rind}(ind10);
  
  epsilon=Geomdat.rnorm(rind)*Geomdat.minorradius/Geomdat.R00(rind);
  xkn=abs(b_10/epsilon);

  [trappedfrac,FSAB,FSAB2]=trappedPartFrac(Geomdat,rind,100,100);

  B00norm=1.0;
  fprintf(dk_fid,'%7.4f%8.4f%5.1f%8.4f%8.5f%8.5f%8.5f  r,R,B,io,xkn,ft,<b^2>\n',...
          Geomdat.rnorm(rind)*Geomdat.minorradius,...
          Geomdat.R00(rind),...
          B00norm,Geomdat.iota(rind),xkn,trappedfrac,FSAB2/Geomdat.B00(rind)^2);
  if rind==1
    fprintf(dk_fid,' c         eps_eff     g11_ft      efield_u    g11_er    ex_er\n');
  end
  fprintf(dk_fid,' cfit  %12.3E%12.3E%12.3E%12.3E%7.3f\n',...
          d11fits.eps_efh(rind),d11fits.g11_ft(rind),d11fits.er_u(rind),...
          d11fits.g11_er(rind),d11fits.ex_er(rind));
  %for ind=1:dkdata.Nruns(rind)
  %  fprintf(dk_fid,'%s\n',dkdata.line1{rind}{ind});
  %  fprintf(dk_fid,'%s\n',dkdata.line2{rind}{ind});
  %end
  for cmulind=1:dkdata.cmulfirst{rind}.Ncmuls
    for efldlind=1:dkdata.cmulfirst{rind}.Neflds(cmulind)
      runind=dkdata.cmulfirst{rind}.cmulinds{cmulind}(efldlind);
      fprintf(dk_fid,'%s\n',dkdata.line1{rind}{runind});
      fprintf(dk_fid,'%s\n',dkdata.line2{rind}{runind});
    end
  end
  fprintf(dk_fid,' e\n');
end
%%%%%%%%%%%%%%%% End of Loop over radii %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% close open files
fclose(dk_fid);
disp(' FINISHED')