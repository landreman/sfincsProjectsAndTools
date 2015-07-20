function [runs,miss]=plot_radscanresults(directory)

if directory(end-3:end)=='.dat'
  dirlist={};
  fid = fopen(directory);
  tline = fgetl(fid);
  while ischar(tline)
    if tline(1)~='%' && tline(1)~='!'
      dirlist={dirlist{:},tline};
    end
    tline = fgetl(fid);
  end
  fclose(fid);
  [runs,miss]=getresults(dirlist);
else
  [runs,miss]=getresults(directory);
end

nbar=1e20;
Tbar=1.6022e-19*1e3;
pbar=nbar*Tbar;

fig(1)

plot(runs.rN,pbar*runs.NTV(:,1),'r',runs.rN,pbar*runs.NTV(:,2),'b',...
     runs.rN,pbar*runs.NTVfromFlux(:,1),'r--',runs.rN,pbar*runs.NTVfromFlux(:,2),'b--')
legend('NTV e','NTV i','NTVfromFlux e','NTVfromFlux i')
title('Terms of NTV')
xlabel('r / a')

fig(2)

s=sqrt(runs.rN);
NTVtot=pbar*(runs.NTV(:,1)+runs.NTV(:,2));
%plot(runs.rN,NTVtot,'g',...
%     runs.rN,pbar*(runs.NTVfromFlux(:,1)+runs.NTVfromFlux(:,2)),'g--')
plot(s,NTVtot,'g',...
     s,pbar*(runs.NTVfromFlux(:,1)+runs.NTVfromFlux(:,2)),'g--')
title('NTV (both species together)')
xlabel('s')

integr=NTVtot'.*4*pi./runs.FSABHat2.*(runs.GHat+runs.iota.*runs.IHat);

fig(3)
plot(s,integr)
title('integrand for total NTV')
xlabel('s')


%NTVtot0=NTVtot(1)-s(1)*(NTVtot(2)-NTVtot(1))/(s(2)-s(1));
%NTVtot1=NTVtot(end)+(1-s(end))*(NTVtot(end)-NTVtot(end-1))/(s(end)-s(end-1));



if 0 %EXTRAPOLATE
  integr0=integr(1)-s(1)*(integr(2)-integr(1))/(s(2)-s(1));
  integr1=integr(end)+(1-s(end))*(integr(end)-integr(end-1))/(s(end)-s(end-1));
  s=[0,s,1];
  integr=[integr0,integr,integr1];
end

NTV_Nm=trapz(s,integr)


%%%%%%%%%% This is for step 1 %%%%%%%%%%
fig(5)
plot(runs.rN,runs.FSABFlow(:,1),runs.rN,runs.FSABFlow(:,2))
title('FSABFlow')
xlabel('r / a')
legend('e','i')