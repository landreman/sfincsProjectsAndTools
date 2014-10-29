function [runs,P,Nval,missing]=getcomboscanresults(dirpath,sortafter)
% For a combined convergence and parameter scan in Sfincs, 
% this loads the results stored in the subdirectories to "dirpath" 
% typically named baseCase/, Ntheta13/, Ntheta19/, Nx6/,...

% P a struct array with Nval values for the 'sortafter' parameter stored in 
% P(1:Nval).val . For each value of the 'sortafter' parameter there is a list of the
% names of the scanned disretisation parameter in P(i).scanParam. For the j'th such
% disretisation parameter P(i).scanParam{j}, the calculated value of, e.g., tauhat_s
% can be retrieved as 
% runs{P(i).ind(j,1)}.tauhat_s(P(i).ind(j,2))
% where i lies in (1:Nval), and j lies in (1:length(P(i).scanParam))

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

missing=[];
runs={};
rind=0;
for pind=1:Nparam
  if nargin==2
    [run,miss]=getresults([dirpath,paramlist{pind},'/'],sortafter);    
  else
    [run,miss]=getresults([dirpath,paramlist{pind},'/']);
  end
  missing=[missing,miss];
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

if  length(runs)==0
  disp('No runs were found!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make the list P of distinct sortafter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(runs)==0
  Nval=0;
  P=[];
else
  for rind=1:length(runs)
    runpvals=getfield(runs{rind},sortafter);
    for eind=1:runs{rind}.NumElements
      if rind==1 && eind==1
        valind=1;
        pvals=getfield(runs{rind},sortafter);
        P(valind).val=pvals(1);
        P(valind).ind=[1,1]; %rind,eind
        P(valind).scanParam={};
        P(valind).scanParam{1}=getfield(runs{rind},'scanParam');
      else      
        
        found=0;
        for vind=1:valind
          if (abs((runpvals(eind)-P(vind).val)/P(vind).val)<1e-3) || ...
                (runpvals(eind)==0 && P(vind).val==0)
            %We consider the values to be the same
            %Add the index of that run to the list for this value
            P(vind).ind=[P(vind).ind; rind,eind]; 
            P(vind).scanParam{length(P(vind).scanParam)+1}=...
                getfield(runs{rind},'scanParam');
            found=1;
          end
        end
        if not(found)
          %New value
          valind=valind+1;
          P(valind).val=runpvals(eind);
          P(valind).ind=[rind,eind];
          P(valind).scanParam={};
          P(valind).scanParam{1}=getfield(runs{rind},'scanParam');
        end
      end
    end
  end
  Nval=valind;
end