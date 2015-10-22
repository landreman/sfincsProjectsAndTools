function [minutes,runnumbers]=readTime(thepath)

if thepath(end)=='/'
  thepath=thepath(1:end-1);
end

lista = dir(thepath);

errfiles=lista([]);
minutes=[];
runnumbers=[];
eind=0;
for ind=1:length(lista)
  if length(lista(ind).name) > 5
    if lista(ind).name(1:4)=='err.'
      eind=eind+1;
      runnumbers(eind,1)=str2num(lista(ind).name(5:end));
      minutes(eind,1)=-1; %default value
      fid = fopen([thepath,'/',lista(ind).name]);
      tline = fgetl(fid);
      while ischar(tline)
        if length(tline) > 5
          if tline(1:4)=='real'
            mind=strfind(tline,'m');
            minutes(eind,1)=str2num(tline(5:mind-1));
          end
        end
        tline = fgetl(fid);
      end
      fclose(fid);
    end
  end
end

%Put the latest run first:
%disp('-----------------------------------------')
%runnumbers
%minutes


[runnumbers,ind]=sort(runnumbers);
minutes=minutes(ind);

minutes=flipud(minutes);
runnumbers=flipud(runnumbers);


%runnumbers
%minutes
