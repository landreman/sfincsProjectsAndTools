function [t,p,files]=runtime(rundirpath)

if not(rundirpath(end)=='/')
  rundirpath=[rundirpath,'/'];
end
list=dir(rundirpath);
files=[];

found=0;
for ind=1:length(list)
  if strcmp(list(ind).name,'sfincsOutput.h5')
    found=1;
    output_datevec=datevec(list(ind).datenum);
  end
end
if not(found)
  t=NaN;
  return
  %error('No sfincsOutput.h5 in that directory')
end
%datestr(output_datevec)


filind=0;
for ind=1:length(list)
  %list(ind).name
  %strfilind(list(ind).name,'.o')
  if not(isempty(strfind(list(ind).name,'.o')))
    %datestr(list(ind).datenum)
    %abs(etime(datevec(list(ind).datenum),output_datevec))
    if abs(etime(datevec(list(ind).datenum),output_datevec))<5*60 %< 5 min between storage
      filind=filind+1;
      files(filind).filename=list(ind).name;
      files(filind).starttime=list(ind).datenum;
      vec=datevec(files(filind).starttime);
      files(filind).starthh=vec(4);
      files(filind).startmm=vec(5);
      files(filind).startss=vec(6);
    end
  end
end

Nfiles=filind;
t=NaN*zeros(Nfiles,1);
p=NaN*zeros(Nfiles,1);

for ind=1:Nfiles
  fid=fopen([rundirpath,files(ind).filename],'r');
  if fid<0
    error(['Could not open ',rundirpath,files(ind).filename])
  else
    %disp(['Opened ',rundirpath,files(ind).filename])
  end
  
  foundtime=0;
  foundproc=0;
  while not(feof(fid) || (foundtime && foundproc))
    str=fgetl(fid);
    if not(isempty(strfind(str,'Start time: ')))
      foundtime=1;
      timestr=str(14:end);
    end
    if not(isempty(strfind(str,'Parallel job (')))
      foundproc=1;
      procstr=str(16:20);
      %disp(['xxx',procstr,'xxx'])
    end
  end
  if foundtime
    hh=str2num(timestr(1:2));
    mm=str2num(timestr(3:4));
    ss=str2num(timestr(5:6));
    s=ss+60*mm+3600*hh;
    files(ind).runtime=files(ind).starthh*3600+files(ind).startmm*60+files(ind).startss-s;
    t(ind)=files(ind).runtime;
    if t(ind)<0
      t(ind)=t(ind)+24*3600;
    end
    tmp=t(ind);
    runss=rem(tmp,60);
    tmp=(tmp-runss)/60;
    runmm=rem(tmp,60);
    tmp=(tmp-runmm)/60;
    files(ind).runhh=rem(tmp,60);
    files(ind).runmm=runmm;
    files(ind).runss=runss;
  end
  if foundproc
    p(ind)=str2num(procstr);
  end  
  fclose(fid);
end

good=find(not(isnan(t)));
t=t(good);
files=files(good);
