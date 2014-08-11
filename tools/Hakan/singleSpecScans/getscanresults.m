function runs=getscanresults(dirpath,sortafter)
% For a combined convergence and parameter scan in Sfincs, 
% this loads the results stored in the subdirectories to "dirpath" 
% typically named baseCase/, Ntheta13/, Ntheta19/, Nx6/,...

%addpath('~/Forskning/Stellarator/sfincs/matlab')

list=dir(dirpath);
if not(dirpath(end)=='/')
  dirpath=[dirpath,'/'];
end

paramlist={};
pind=0;
for ind=1:length(list)
  if list(ind).isdir
    if list(ind).name(1)~='.'
      pind=pind+1;
      paramlist{pind}=list(ind).name;
    end
  end
end
Nparam=pind;

runs={};
rind=0;
for pind=1:Nparam
  if nargin==2
    run=getresults([dirpath,paramlist{pind},'/'],sortafter);    
  else
    run=getresults([dirpath,paramlist{pind},'/']);
  end
  if run.NumElements>0
    rind=rind+1;
    runs{rind}=run;
    runs{rind}.path=[dirpath,paramlist{pind},'/'];
    runs{rind}.dirname=paramlist{pind};
    foundnum=0;
    for nind=1:length(runs{rind}.dirname)
      if not(foundnum)
        charnum=str2num(runs{rind}.dirname(nind));
        foundnum=not(isempty(charnum) || runs{rind}.dirname(nind)=='i'...
                     || runs{rind}.dirname(nind)=='j');
        if not(foundnum)
          len=nind;
        end
      end
    end
    runs{rind}.scanParam=runs{rind}.dirname(1:len);
  end
end
  
