function runs_out=loadallh5(dirpath,sortafter)
%
% If "dirpath" is a directory name, then the routine 
% loads all the sfincsOutput.h5 files stored in the subdirectories to
% "dirpath" which have numerical names (typically called 00, 01, 02,...)
% or which have names beginning with "baseCase", "Ntheta", ...
%
% If "dirpath" is a cell array this routine loads all the
% sfincsOutput.h5 files stored in
% the directories listed in the array.
%

presentdir=pwd;

if nargin==0
  dirpath=presentdir;
end

if not(iscell(dirpath)) 
  if not(isempty(strfind(presentdir,dirpath(2:end-1))))
    if dirpath(end)=='/'
      tmp=dirpath(1:end-1);
    else
      tmp=dirpath;
    end
    if tmp(end-3:end)==presentdir(end-3:end)
      dirpath=pwd %We are already standing in that directory
    end
  end
  if not(exist(dirpath,'dir'))
    if dirpath(1)=='/'
    dirpath=[getenv('SFINCS_HOME'),'/fortran/version3',dirpath];
    else    
      dirpath=[getenv('SFINCS_HOME'),'/fortran/version3/',dirpath];
    end
  end
  list=dir(dirpath);
  if not(dirpath(end)=='/')
    dirpath=[dirpath,'/'];
  end

  if isempty(list)
    error(['Nothing in the directory ',dirpath])
  end
  
  runlist={};
  runind=0;
  for ind=1:length(list)
    if list(ind).name(1)~='.'
      if not(isempty(str2num(list(ind).name))) || ...
            strcmp(list(ind).name,'baseCase') || ...
            list(ind).name(1)=='N' || ...
            strcmp(list(ind).name(1:2),'Er') || ...
            strcmp(list(ind).name(1:4),'dPhi') || ...
            not(isempty(strfind(list(ind).name,'solverTolerance'))) || ...
            not(isempty(strfind(list(ind).name,'nhats')))
        if not(exist([dirpath,list(ind).name,'/sfincsOutput.h5'],'file'))
          disp(['Skipping the directory ',dirpath,list(ind).name])
          disp('The file sfincsOutput.h5 did not exist!')
        else         
          runind=runind+1;
          runlist{runind}=list(ind).name;
        end
      end
    end
  end
  
  
  
  if runind==0
    runs{1}='no subdirs';
    warning(['No subdirectories in ',dirpath])
  else
    for ind=1:runind
      try
        runs{ind}=h5load([dirpath,runlist{ind},'/sfincsOutput.h5']);
        runs{ind}.dirpath=dirpath;
        runs{ind}.rundir=runlist{ind};
      catch err
        if strcmp(err.identifier,'MATLAB:imagesci:hdf5lib:libraryError')
          runs{ind}=[runlist{ind},' error'];
        elseif strcmp(err.identifier,'MATLAB:imagesci:validate:fileOpen')
          runs{ind}=[runlist{ind},' empty'];
        else
          disp(['Problem opening the file: ',dirpath,runlist{ind},'/sfincsOutput.h5'])
          disp('rethrown error:')
          rethrow(err);
        end
      end
    end
  end
  
else %Input was a cell array
  runlist=dirpath;
  for ind=1:length(runlist)
    if not(isempty(strfind(pwd,runlist{ind}(2:end-1))))
      runlist{ind}=pwd; %We are already standing in that directory
    end
    if not(exist(runlist{ind},'dir'))
      if runlist{ind}(1)=='/'
        runlist{ind}=[getenv('SFINCS_HOME'),'/fortran/version3',runlist{ind}];
      else    
        runlist{ind}=[getenv('SFINCS_HOME'),'/fortran/version3/',runlist{ind}];
      end
    end
    try
      runs{ind}=h5load([runlist{ind},'/sfincsOutput.h5']);
      runs{ind}.dirpath='';
      runs{ind}.rundir=runlist{ind};
    catch err
      if strcmp(err.identifier,'MATLAB:imagesci:hdf5lib:libraryError')
        runs{ind}=[runlist{ind},' error'];
      elseif strcmp(err.identifier,'MATLAB:imagesci:validate:fileOpen')
        runs{ind}=[runlist{ind},' empty'];
      else
        disp('rethrown error:')
        rethrow(err);
      end
    end
  end
end

if nargin<2 %do not sort
  runs_out=runs;
else %sort the loaded results
  N=length(runs);
  good=[];
  vals=[];
  notgood=[];
  for ind=1:N
    if isfield(runs{ind},sortafter)
      good=[good,ind];
      vals=[vals,getfield(runs{ind},sortafter)];
    else
      notgood=[notgood,ind];
    end
  end
  [dummy,perm]=sort(vals);
  good=good(perm);
  runs_out={runs{[good,notgood]}};
end