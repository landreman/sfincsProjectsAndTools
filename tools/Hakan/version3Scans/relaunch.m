function relaunch(thefullpath)

%This is for draco. Use relaunch_hydra.m on hydra

overheadLossFactor=0.93; %Estimate how much is left of node memory when overhead is removed

if thefullpath(end)=='/'
  thefullpath=thefullpath(1:end-1);
end

if not(exist(thefullpath,'dir'))
  thefullpath=[getenv('SFINCS_HOME'),'/fortran/version3/',thefullpath];
end
listing = dir(thefullpath)

if isempty(listing)
  error('No subdirectories exist!')
end

disp('Subdirectories:')
disp('---------------------')
gind=0;
dirs=listing([]);
isok=[];
for ind=1:length(listing)
  if listing(ind).isdir

    if not(strcmp(listing(ind).name,'.') || strcmp(listing(ind).name,'..'))

      if not(isempty(str2num(listing(ind).name))) || ...
            strcmp(listing(ind).name,'baseCase') || ...
            listing(ind).name(1)=='N'|| ...
            strcmp(listing(ind).name(1:4),'dPhi')|| ...
            not(isempty(strfind(listing(ind).name,'solverTolerance')))|| ...
            not(isempty(strfind(listing(ind).name,'nhats')))
        if inQueueOrRunning_slurm([thefullpath,'/',listing(ind).name])
          disp([listing(ind).name,' : is in queue or running.'])
        else
          gind=gind+1;
          dirs(gind)=listing(ind);
          %[thefullpath,'/',listing(ind).name]
          %readTime([thefullpath,'/',listing(ind).name])
          %[minutes,runnumbers]=readTime([thefullpath,'/',listing(ind).name]);
          %if minutes(1)==-1
          %  disp([dirs(gind).name,' : run #',num2str(runnumbers(1)),', had error,',...
          %        ' was cancelled',...
          %        ' or reached wall clock limit'])
          %else
          %  disp([dirs(gind).name,' : run #',num2str(runnumbers(1)),', time=', ...
          %        num2str(minutes(1)),' minutes'])
	  fid = fopen([thefullpath,'/',listing(ind).name,'/Sfincs.out']);
          if fid<0
            error(['Could not open ',thefullpath,'/',listing(ind).name,'/Sfincs.out'])
          end
          %lind=0;
          tline = fgetl(fid);
          fini=0;
          errorsfound=0;
          while ischar(tline)
            %lind=lind+1;
            if not(isempty(strfind(tline,'ERROR RETURN ** FROM DMUMPS INFO(1)=')))
              errorsfound=1;
            end
            if not(isempty(strfind(tline,'Goodbye')))
              fini=1;
            end
            tline = fgetl(fid);
          end
          fclose(fid);
          
          if fini==0 || errorsfound==1
            disp([dirs(gind).name,' had error,',...
                  ' was cancelled',...
                  ' or reached wall clock limit'])
            isok(gind)=0;
          else
            disp([dirs(gind).name,' is finished'])
            isok(gind)=1;
          end
        end
      end
    end
  end
end

if isempty(dirs) || all(isok)
  error('No subdirectories that need relaunch present!')
end

for ind=1:length(dirs)
  dirs(ind).inputnamelistFile=[dirs(ind).name,'/input.namelist'];
  dirs(ind).jobFile=[dirs(ind).name,'/job.sfincsScan'];
  dirs(ind).sfincsOutFile=[dirs(ind).name,'/Sfincs.out'];
  if not(exist([thefullpath,'/',dirs(ind).inputnamelistFile],'file'))
    error(['The file input.namelist is missing in ',dirs(ind).name,])
  end
  if not(exist([thefullpath,'/',dirs(ind).jobFile],'file'))
    error(['The file job.sfincsScan is missing in ',dirs(ind).name,])
  end
  if not(exist([thefullpath,'/',dirs(ind).sfincsOutFile],'file'))
    error(['The file Sfincs.out is missing in ',dirs(ind).name,])
  end
end

%Do not continue unless everything is there.
wcl_h=input('New wall clock limit, hours: ');
wcl_m=input('New wall clock limit, minutes: ');
%newConsumableMemoryGB=input('Which nodes do you want (ConsumableMemory in GB): ');
%if not(newConsumableMemoryGB==120 || newConsumableMemoryGB==56)
%  error('Nonexistent nodes!')
%end
ConsumableMemoryGB=120;
availableMemoryMB=112000;


for ind=1:length(dirs)
  disp('----------------------------------------')
  disp(['Directory name                 : ',dirs(ind).name]);
  %First, read the file Sfincs.out:
  fid = fopen([thefullpath,'/',dirs(ind).sfincsOutFile]);
  dirs(ind).finished=0;
  dirs(ind).errorsfound=0;
  dirs(ind).mbytes=-1;
  lind=0;
  tline = fgetl(fid);
  while ischar(tline)
    lind=lind+1;
    %SfinsfOutLines{lind}=tline;
    if not(isempty(strfind(tline,'Goodbye')))
      dirs(ind).finished=1;
      disp('The calculation is finished')
    end
    if not(isempty(strfind(tline,'ERROR RETURN ** FROM DMUMPS INFO(1)=')))
      dirs(ind).errorsfound=1;
    end
    if not(isempty(strfind(tline,...
              '** TOTAL     space in MBYTES for IC factorization')))
      colonind=strfind(tline,':');
      dirs(ind).mbytes=str2num(tline(colonind+1:end));
      if isempty(dirs(ind).mbytes)
        error(['Error reading MBYTES in ',dirs(ind).sfincsOutFile])
      end
    end
    tline = fgetl(fid);
  end
  fclose(fid);
  
  if dirs(ind).mbytes == -1
    dirs(ind).mbytes=NaN;
  end
  dirs(ind).goodToGo=(not(dirs(ind).finished) || dirs(ind).errorsfound) && dirs(ind).mbytes~=-1;
  
  if dirs(ind).goodToGo
    mumpsGB = dirs(ind).mbytes/1e3;
    
    fid = fopen([thefullpath,'/',dirs(ind).jobFile]);
    node_line_ind=-1;
    taskspernode_line_ind=-1;
    wallclock_line_ind=-1;
    ConsMem_line_ind=-1;
    launch_line_ind=-1;
    prev_wall_clock_limit_str='';
    prevNumNodes=-1;
    lind=1;
    tline = fgetl(fid);
    joblines={}; %clear old data.
    joblines{lind}=tline;
    cpupertask_line_ind=NaN;
    ompnumthreads_line_ind=NaN;
    
    while ischar(tline)
      %if not(isempty(strfind(tline,'ConsumableMemory')))
      %  ConsMem_line_ind=lind;
      %  leftparind =strfind(tline,'(');
      %  gbind=strfind(tline,'gb');
      %  if isempty(leftparind) || isempty(gbind)
      %    error(['leftparind or gbind is empty! In ',dirs(ind).jobFile])
      %  end
      %  ConsumableMemoryGB=str2num(tline(leftparind+1:gbind-1));
      %  if isempty(ConsumableMemoryGB)
      %    error(['Error readin ConsumableMemory in ',dirs(ind).jobFile])
      %  end
      %end
      if not(isempty(strfind(tline,'#SBATCH --nodes=')))
        node_line_ind=lind;
        prevNumNodes=str2num(tline(strfind(tline,'=')+1:end));
      end
      if not(isempty(strfind(tline,'#SBATCH --ntasks-per-node=')))
        taskspernode_line_ind=lind;
        tasksPerNode=str2num(tline(strfind(tline,'=')+1:end));
      end
      if not(isempty(strfind(tline,'#SBATCH --time=')))
        wallclock_line_ind=lind;
        prev_wall_clock_limit_str=tline(strfind(tline,'=')+1:end);
      end
      if not(isempty(strfind(tline,'#SBATCH --cpus-per-task=')))
        cpupertask_line_ind=lind;
        prev_cpupertask=str2num(tline(strfind(tline,'=')+1:end));
      end
      if not(isempty(strfind(tline,'export OMP_NUM_THREADS=')))
        ompnumthreads_line_ind=lind;
        prev_ompnumthreads=str2num(tline(strfind(tline,'=')+1:end));
      end
      if not(isempty(strfind(tline,'#SBATCH --partition=')))
        partition_line_ind=lind;
        prev_partition=tline(strfind(tline,'=')+1:end);
      end

      if length(tline)>1
        if tline(1)~='#'
          tmpind=strfind(tline,'-mat_mumps_icntl_23');
          if not(isempty(tmpind))
            launch_line_ind=lind;
            tmp = textscan(tline, '%s', 'delimiter', ' ');  
            lineparts=tmp{1};
            flagind=find(strcmp(lineparts,'-mat_mumps_icntl_23'));
            lineparts{flagind+1}=num2str(ceil(availableMemoryMB/tasksPerNode));
            newLaunchLine=lineparts{1};
            for indx=2:length(lineparts)
              newLaunchLine=[newLaunchLine,' ',lineparts{indx}];
            end
          end
        end
      end
      lind=lind+1;
      tline = fgetl(fid);
      joblines{lind}=tline;
    end    
    fclose(fid);
    if prevNumNodes==-1
      error(['Could not read number of nodes in ',dirs(ind).jobFile])
    end
    if wallclock_line_ind==-1
      error(['Could not read wall_clock_limit in ',dirs(ind).jobFile])
    end
    %if ConsMem_line_ind==-1
    %  error(['Could not read ConsumableMemory in ',dirs(ind).jobFile])
    %end
    
    newNumNodes=ceil(mumpsGB/(ConsumableMemoryGB*overheadLossFactor));

    disp(['Previously used wall_clock_limit: ',prev_wall_clock_limit_str]);
    disp(['New chosen used wall_clock_limit: ',...
          num2str(wcl_h,'%02d'),':',num2str(wcl_m,'%02d'),':00']);
    %disp(['Previous ConsumableMemory       : ',num2str(ConsumableMemoryGB),' GB']);
    %disp(['New ConsumableMemory            : ',num2str(newConsumableMemoryGB),' GB']);
    disp(['MUMPS memory requirement        : ',num2str(mumpsGB),' GB']);
    disp(['Previously used number of nodes : ',num2str(prevNumNodes)]);
    disp(['Recommended number of nodes     : ',num2str(newNumNodes)]);

    if newNumNodes*ConsumableMemoryGB<=prevNumNodes*ConsumableMemoryGB
      disp(['Recommended #nodes*ConsumableMemory <= used #nodes*ConsumableMemory ! ',...
            ' Recommend manual launch with higher wall_clock_limit!'])
      %disp(['Recommended #nodes <= used #nodes !  Not launching ',dirs(ind).jobFile,...
      %      '. Recommend manual launch with higher wall_clock_limit!'])
    end
    %  dirs(ind).goodToGo=0;
    %else
      numnodeanswer=input(['Do you want ',num2str(newNumNodes),...
          ' nodes? (ret=yes, n=do not launch, or give the wanted number) :'],'s');
      if strcmp(numnodeanswer,'n')
        dirs(ind).goodToGo=0;
      else
        if not(isempty(str2num(numnodeanswer)))
          wantedNewNumNodes=str2num(numnodeanswer);
        else
          wantedNewNumNodes=newNumNodes;
        end
        % Write cell joblines into job file
        disp(['CHANGING THE FILE ',dirs(ind).jobFile])
        joblines{node_line_ind}     =['#SBATCH --nodes=',num2str(wantedNewNumNodes)];
        %joblines{ConsMem_line_ind}  =...
        %    ['# @ node_resources = ConsumableMemory(',num2str(newConsumableMemoryGB),'gb)'];
        joblines{wallclock_line_ind}=...
            ['#SBATCH --time=',num2str(wcl_h,'%02d'),':',num2str(wcl_m,'%02d'),':00'];
        if not(isnan(cpupertask_line_ind))
          joblines{cpupertask_line_ind}=...
              ['#SBATCH --cpus-per-task=',num2str(floor(32/tasksPerNode))];
        end
        
        if not(isnan(ompnumthreads_line_ind))
          joblines{ompnumthreads_line_ind}=...
              ['#export OMP_NUM_THREADS=',num2str(floor(32/tasksPerNode))];
        end
        
        if wcl_h==0 && wcl_m<=30
	  partitionstr='express';
	elseif wcl_h<=3 || (wcl_h==4 && wcl_m==0)
	  partitionstr='short';
	else
          partitionstr='general';
	end
        joblines{partition_line_ind}=...
	  ['#SBATCH --partition=',partitionstr];

        joblines{launch_line_ind}= newLaunchLine;
        
        fid = fopen([thefullpath,'/',dirs(ind).jobFile], 'w');
        for i = 1:numel(joblines)
          if joblines{i+1} == -1
            fprintf(fid,'%s', joblines{i});
            break
          else
            fprintf(fid,'%s\n', joblines{i});
          end
        end 
        fclose(fid);
      end
    %end
  end
end
disp('----------------------------------------')

numbertolaunch=0;
for ind=1:length(dirs)
  numbertolaunch=numbertolaunch + dirs(ind).goodToGo;
end

if numbertolaunch==0
  disp('There are no jobs to launch.')
else
  answer=input(['Do you want to launch the ',...
                num2str(numbertolaunch),' jobs in the changed files (y/n): '],'s');
  while not(strcmp(answer,'y') || strcmp(answer,'n'))
    answer=input('Please answer y or n: ','s');
  end
  
  if strcmp(answer,'y') %The user wants to launch the jobs
    currdir=pwd;
    cd(thefullpath);
    for ind=1:length(dirs)
      if dirs(ind).goodToGo
        cd(dirs(ind).name)
        status=system('sbatch job.sfincsScan');
        if status~=0
          disp(['Could NOT submit the job: ',thefullpath,'/',dirs(ind).jobFile]);
        else
          disp(['Submitted the job       : ',thefullpath,'/',dirs(ind).jobFile]);
          %pwd
        end
        cd('..')
      end
    end
    cd(currdir);
    [status,sdirresult]=system('sdir');
    disp(sdirresult)
  end
end


