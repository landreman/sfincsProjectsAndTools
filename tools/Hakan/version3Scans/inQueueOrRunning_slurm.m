function running=inQueueOrRunning(thepath)

%[status,longlist]=system(['llq -l -u ',getenv('USER')]);
%[status,longlist]=system(['squeue -u ',getenv('USER')]);
[status,longlist]=system(['sdir']);

longlist=longlist(find(longlist~=char(10)));
longlist=longlist(find(longlist~=char(13)));

if iscell(thepath)
  for pind=1:length(thepath)
    start=findstr(thepath{pind},'runs');
    if not(isempty(start))
      thispath=thepath{pind}(start(1):end);
    else
      thispath=thepath{pind};
    end
    running(pind)=not(isempty(findstr(longlist,thispath)));
  end
else
   start=findstr(thepath,'runs');
    if not(isempty(start))
      thispath=thepath(start(1):end);
    else
    thispath=thepath;
  end
  running=not(isempty(findstr(longlist,thispath)));
end
