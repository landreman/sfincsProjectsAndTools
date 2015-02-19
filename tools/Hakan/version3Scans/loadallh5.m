function runs=loadallh5(dirpath)
% Loads all the sfincsOutput.h5 files stored in the subdirectories to
% "dirpath" which have numerical names (typically called 00, 01, 02,...)
% or which have names beginning with "baseCase", "Ntheta", ...

if nargin==1
  list=dir(dirpath);
  if not(dirpath(end)=='/')
    dirpath=[dirpath,'/'];
  end
else
  list=dir;
  dirpath='';
end

runlist={};
runind=0;
for ind=1:length(list)
  if not(isempty(str2num(list(ind).name))) || ...
        strcmp(list(ind).name,'baseCase') || ...
        list(ind).name(1)=='N'
    runind=runind+1;
    runlist{runind}=list(ind).name;
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
        disp('rethrown error:')
        rethrow(err);
      end
    end
  end
end