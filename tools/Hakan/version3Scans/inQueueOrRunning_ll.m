function running=inQueueOrRunning(thepath)

[status,longlist]=system(['llq -l -u ',getenv('USER')]);

longlist=longlist(find(longlist~=char(10)));
longlist=longlist(find(longlist~=char(13)));

if iscell(thepath)
  for pind=1:length(thepath)
    if thepath{pind}(1)=='~';
      thispath=thepath{pind}(2:end);
    else
      thispath=thepath{pind};
    end
    running(pind)=not(isempty(findstr(longlist,thispath)));
  end
else
  if thepath(1)=='~';
    thispath=thepath(2:end);
  else
    thispath=thepath;
  end
  running=not(isempty(findstr(longlist,thispath)));
end
