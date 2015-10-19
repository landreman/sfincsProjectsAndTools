function compareSingleMulti(singleDir,multiDir)

if singleDir(end-3:end)=='.dat'
  dirlist={};
  fid = fopen(singleDir);
  tline = fgetl(fid);
  while ischar(tline)
    if tline(1)~='%' && tline(1)~='!'
      dirlist={dirlist{:},tline};
    end
    tline = fgetl(fid);
  end
  fclose(fid);
  [runsS,missS]=getresults(dirlist);
else
  [runsS,missS]=getresults(singleDir);
end

if multiDir(end-3:end)=='.dat'
  dirlist={};
  fid = fopen(multiDir);
  tline = fgetl(fid);
  while ischar(tline)
    if length(tline)>0
      if tline(1)~='%' && tline(1)~='!'
        dirlist={dirlist{:},tline};
      end
    end
    tline = fgetl(fid);
  end
  fclose(fid);
  [runsM,missM]=getresults(dirlist);
else
  [runsM,missM]=getresults(multiDir);
end

Rbar=1;
nbar=1e20;
Tbar=1.6022e-19*1e3;
pbar=nbar*Tbar;
psiAHat=runsM.psiAHat(1);
if runsS.psiAHat(1) ~= psiAHat
  error(' psiAHat is different in the two inputs! ')
end

fig(1)

  plot(runsM.rN,pbar*runsM.NTV(:,1),'g-.',...
       runsM.rN,pbar*runsM.NTVfromFlux(:,1),'c-.',...
       runsM.rN,pbar*runsM.NTV(:,2),'r',...
       runsM.rN,pbar*runsM.NTVfromFlux(:,2),'b',...
       runsS.rN,pbar*runsS.NTV(:,1),'k--',...
       runsS.rN,pbar*runsS.NTVfromFlux(:,1),'m--')
  
  legend('NTV e','NTVfromFlux e',...
         'NTV i','NTVfromFlux i',...
         'NTV i, single','NTVfromFlux i, single')
%title('Terms of NTV')
xlabel('r / a')

fig(2)

sS=runsS.rN.^2;
sM=runsM.rN.^2;


NTVtotM=pbar*(runsM.NTV(:,1)+runsM.NTV(:,2));
NTVfromFluxtotM=pbar*(runsM.NTVfromFlux(:,1)+runsM.NTVfromFlux(:,2));

NTVtotS=pbar*runsS.NTV(:,1);
NTVfromFluxtotS=pbar*runsS.NTVfromFlux(:,1);

plot(sM,NTVtotM,'r',...
     sM,NTVfromFluxtotM,'b',...
     sS,NTVtotS,'k--',...
     sS,NTVfromFluxtotS,'m--')
title('NTV total (all species together)')
xlabel('s')
legend('multi direct','multi fromFlux','single direct','single fromFlux')

integrM=NTVtotM'.*abs(4*pi^2./runsM.FSABHat2.*...
                      (runsM.GHat+runsM.iota.*runsM.IHat)*...
                      psiAHat*Rbar^3);
integrS=NTVtotS'.*abs(4*pi^2./runsS.FSABHat2.*...
                      (runsS.GHat+runsS.iota.*runsS.IHat)*...
                      psiAHat*Rbar^3);
%abs is taken because previously psiAHat had the wrong sign!

fig(3)
plot(sM,integrM,'r-',sS,integrS,'k--')
title('integrand for total NTV')
xlabel('s')
legend('multi','single')
NTV_Nm_multi =trapz(sM,integrM)
NTV_Nm_single=trapz(sS,integrS)


%%%%%%%%%% This is for step 1 %%%%%%%%%%
fig(5)

plot(runsM.rN,runsM.FSABFlow(:,1),runsM.rN,runsM.FSABFlow(:,2),...
     runsS.rN,runsS.FSABFlow(:,1),'--')
legend('e','i','i single spec.',4)
title('FSABFlow')
xlabel('r / a')

