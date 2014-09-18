function runs=plotscanresults(dirpath,param)
%This funtion plots the results of a combined convergence and parameter scan
%"dirpath" is the path to the directory where "baseCase" and the different directories
%for convergence tested parameters.
%"param" is the name of the scanned parameter

%dirpath='~/Forskning/Stellarator/sfincs/gitsfincs/fortran/singleSpecies/runs/ntv/AUG28061/r0p6/Er23kVconvScan';
%dirpath='~/Forskning/Stellarator/sfincs/gitsfincs/fortran/singleSpecies/runs/ntv/AUG28061/r0p6/Er0convScan';
%dirpath='~/Forskning/Stellarator/sfincs/gitsfincs/fortran/singleSpecies/runs/ntv/AUG28061/r0p3';
%dirpath='~/Forskning/Stellarator/sfincs/gitsfincs/fortran/singleSpecies/runs/ntv/AUG28061r0p3_old';
%dirpath='~/Forskning/Stellarator/sfincs/gitsfincs/fortran/singleSpecies/runs/ntv/AUG28061/r0p3DKES';
%dirpath=['~/Forskning/Stellarator/sfincs/gitsfincs/fortran/singleSpecies/runs/' ...
%         'fig2ErScan_ecb2_conv/']; %I am running the base case with the new NTVMatrix storage
%param='EStar';
%param='dPhiHatdpsi';

[runs,P,Nval,missing]=getscanresults(dirpath,param);

if not(runs{1}.RHSMode==2)
  error('I assumed RHSMode==2 but it is not!')
end
if length(missing)>0
  disp('Missing in action:')
  for mind=1:length(missing)
    disp(missing(mind).dir)
  end
end

rat1=NaN*zeros(1,Nval);
rat2=NaN*zeros(1,Nval);
rat3=NaN*zeros(1,Nval);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
showparam={'Ntheta','Nzeta','Nx','Nxi'};%,'solverTolerance'};
Nrow=length(showparam);
Ncol=6+1;

makeplots=1
if makeplots
close all
%For each value of 'param', make different scan variable plots
for vind=1:Nval
  fig(50+vind)
  pos=get(gcf,'Position');
  pos(3:4)=[200*Ncol,200*Nrow];
  set(gcf,'Position',pos);
  for shind=1:Nrow
    good=[];
    Ngood=0;base=[];
    for gind=1:length(P(vind).scanParam)
      if strcmp(P(vind).scanParam(gind),showparam{shind})
        good=[good,gind];
        Ngood=Ngood+1;
      elseif strcmp(P(vind).scanParam(gind),'baseCase')
        good=[good,gind];
        Ngood=Ngood+1;
        base=Ngood;
      end
    end
    goodind=P(vind).ind(good,:);
    
    L=NaN*zeros(Ngood,3,3);
    LNTV=NaN*zeros(Ngood,3);
    %tauhat_s=zeros(Ngood,1);
    scanval=zeros(Ngood,1);
    for n=1:Ngood
      scanvalsrun=getfield(runs{goodind(n,1)},showparam{shind});
      scanval(n)=scanvalsrun(goodind(n,2));
      L(n,:,:)=runs{goodind(n,1)}.transportMatrix(goodind(n,2),:,:);
      if isfield(runs{goodind(n,1)},'NTVMatrix')
        LNTV(n,:)=runs{goodind(n,1)}.NTVMatrix(goodind(n,2),:);
      end
      %tauhat_s(n)=runs{goodind(n,1)}.tauhat_s(goodind(n,2));
    end
   
    other=[1:base-1,base+1:Ngood];
    
    Lind=1; %L11
    subplot(Nrow,Ncol,(shind-1)*Ncol+Lind)
    plot(scanval(base),L(base,1,1),'r+',scanval(other),L(other,1,1),'g+',...
         scanval(base),LNTV(base,1),'ro',scanval(other),LNTV(other,1),'go')
    xlabel(showparam{shind})
    ylabel('L11 +, LNTV1 o')
    title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L12
    subplot(Nrow,Ncol,(shind-1)*Ncol+Lind)
    plot(scanval(base),L(base,1,2),'r+',scanval(other),L(other,1,2),'g+',...
         scanval(base),LNTV(base,2),'ro',scanval(other),LNTV(other,2),'go')
    xlabel(showparam{shind})
    ylabel('L12 +, LNTV2 o')
    title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L22
    subplot(Nrow,Ncol,(shind-1)*Ncol+Lind)
    plot(scanval(base),L(base,2,2),'r+',scanval(other),L(other,2,2),'g+')%,...
         %scanval(base),LNTV(base,3),'ro',scanval(other),LNTV(other,3),'go')
    xlabel(showparam{shind})
    ylabel('L22')% +, LNTV3 o')
    title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L13
    subplot(Nrow,Ncol,(shind-1)*Ncol+Lind)
    plot(scanval(base),L(base,1,3),'r+',scanval(other),L(other,1,3),'g+')%,...
         %[min(scanval)*0.97,max(scanval)*1.03],-1/runs{1}.iota(1)*[1,1],'k:')
    xlabel(showparam{shind})
    ylabel('L13')
    title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L23
    subplot(Nrow,Ncol,(shind-1)*Ncol+Lind)
    plot(scanval(base),L(base,2,3),'r+',scanval(other),L(other,2,3),'g+')%,...
        % [min(scanval),max(scanval)],-4.57*/runs{1}.iota*[1,1],'k:')
        %Matt's value 4.57 could be wrong here
    xlabel(showparam{shind})
    ylabel('L23')
    title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L33
    subplot(Nrow,Ncol,(shind-1)*Ncol+Lind)
    plot(scanval(base),L(base,3,3),'r+',scanval(other),L(other,3,3),'g+')
    xlabel(showparam{shind})
    ylabel('L33')
    title([param,' = ',num2str(P(vind).val)])
    
    %Lind=Lind+1;
    %subplot(Nrow,Ncol,(shind-1)*Ncol+Lind)
    %plot(scanval(base),LNTV(base,1),'r+',scanval(other),LNTV(other,1),'g+')
    %xlabel(showparam{shind})
    %ylabel('LNTV_1')
    %title([param,' = ',num2str(P(vind).val)])
    
    %Lind=Lind+1;
    %subplot(Nrow,Ncol,(shind-1)*Ncol+Lind)
    %plot(scanval(base),LNTV(base,2),'r+',scanval(other),LNTV(other,2),'g+')
    %xlabel(showparam{shind})
    %ylabel('LNTV_2')
    %title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1;
    subplot(Nrow,Ncol,(shind-1)*Ncol+Lind)
    plot(scanval(base),LNTV(base,3),'r+',scanval(other),LNTV(other,3),'g+')
    xlabel(showparam{shind})
    ylabel('LNTV_3')
    title([param,' = ',num2str(P(vind).val)])
    
  end
  if not(isempty(base))
    rat1(vind)=LNTV(base,1)/L(base,1,1);
    rat2(vind)=LNTV(base,2)/L(base,1,2);
    rat3(vind)=LNTV(base,3)/L(base,1,3);
  end
end
    
rat1
rat2
rat3

end

%%%%%%%%%%%%%%%%% run time vs discretisation
curvefitworked=0;
Nthetabase=NaN;
Nzetabase=NaN;
Nxibase=NaN;
Nxbase=NaN;
for runind=1:length(runs)
  if strcmp(runs{runind}.scanParam,'baseCase')
    Nthetabase=runs{runind}.Ntheta(1)
    Nzetabase=runs{runind}.Nzeta(1)
    Nxibase=runs{runind}.Nxi(1)
    Nxbase=runs{runind}.Nx(1)
  end
end

time16=[];
Ntheta16=[];
Nzeta16 =[];
Nxi16 =[];
Nx16 =[];
time36 =[];
Ntheta36 =[];
Nzeta36 =[];
Nxi36 =[];
Nx36=[];
timeNN =[];
NthetaNN =[];
NzetaNN =[];
NxiNN =[];
NxNN=[];

ind16=0;
ind36=0;
indNN=0;
for rind=1:length(runs)
  for uind=1:runs{rind}.NumElements
    if runs{rind}.run(uind).proc==16
      ind16=ind16+1;
      time16(ind16)=runs{rind}.run(uind).time/60; %minutes
      Ntheta16(ind16)=runs{rind}.Ntheta(uind);
      Nzeta16(ind16)=runs{rind}.Nzeta(uind);
      Nxi16(ind16)=runs{rind}.Nxi(uind);
      Nx16(ind16)=runs{rind}.Nx(uind);
    elseif runs{rind}.run(uind).proc==36
      ind36=ind36+1;
      time36(ind36)=runs{rind}.run(uind).time/60; %minutes
      Ntheta36(ind36)=runs{rind}.Ntheta(uind);
      Nzeta36(ind36)=runs{rind}.Nzeta(uind);
      Nxi36(ind36)=runs{rind}.Nxi(uind);
      Nx36(ind36)=runs{rind}.Nx(uind);      
    else
      indNN=indNN+1;
      timeNN(indNN)=runs{rind}.run(uind).time/60; %minutes
      NthetaNN(indNN)=runs{rind}.Ntheta(uind);
      NzetaNN(indNN)=runs{rind}.Nzeta(uind);
      NxiNN(indNN)=runs{rind}.Nxi(uind);
      NxNN(indNN)=runs{rind}.Nx(uind);      
    %else
    %  disp(['Number of cores: ',num2str(runs{rind}.run(uind).proc)])
    %  error('Not implemented yet!')
    end
  end
end

Ntot16=Nx16.*Nxi16.*Nzeta16.*Ntheta16;
Ntot36=Nx36.*Nxi36.*Nzeta36.*Ntheta36;
NtotNN=NxNN.*NxiNN.*NzetaNN.*NthetaNN;
time16(find(time16<0))=NaN;
time36(find(time36<0))=NaN;
timeNN(find(timeNN<0))=NaN;

if not(isnan(Nzetabase))
  options.Display='off';lb=[0,-500];ub=[200,500];
  fun = @(k,xdata)k(1)*xdata+k(2);
  try
    if not(isempty(Ntot16))
      k = lsqcurvefit(fun,[200/3e6,0],Ntot16,time16,lb,ub,options);
    else
      k = lsqcurvefit(fun,[200/3e6,0],NtotNN,timeNN,lb,ub,options);  
    end
    curvefitworked=1;
    disp(['minutes per Ntot: ',num2str(k(1))])
    disp(['minutes per Ntheta (crude): ',num2str(k(1)*Nzetabase*Nxibase*Nxbase)])
    disp(['minutes per Nzeta  (crude): ',num2str(k(1)*Nthetabase*Nxibase*Nxbase)])
    disp(['minutes per Nxi    (crude): ',num2str(k(1)*Nthetabase*Nzetabase*Nxbase)])
    disp(['minutes per Nx     (crude): ',num2str(k(1)*Nthetabase*Nzetabase*Nxibase)])
  catch me
    curvefitworked=0;
    disp('Could not do the curve fit.')
  end
  
  if not(isempty(Ntot16))
    inds=find((Nzeta16==Nzetabase) & (Nxi16==Nxibase) & (Nx16==Nxbase));
    if length(inds)>=2
      ktheta=lsqcurvefit(fun,[k(1)*Nzetabase*Nxibase*Nxbase,0],Ntheta16(inds),time16(inds),lb,ub,options);
      disp(['16 cores: minutes per Ntheta: ',num2str(ktheta(1))])
    end
    inds=find((Ntheta16==Nthetabase) & (Nxi16==Nxibase) & (Nx16==Nxbase));
    if length(inds)>=2
      kzeta=lsqcurvefit(fun,[k(1)*Nthetabase*Nxibase*Nxbase,0],Nzeta16(inds),time16(inds),lb,ub,options);
      disp(['16 cores: minutes per Nzeta : ',num2str(kzeta(1))])
    end
    inds=find((Ntheta16==Nthetabase) & (Nzeta16==Nzetabase) & (Nx16==Nxbase));
    if length(inds)>=2
      kxi=lsqcurvefit(fun,[k(1)*Nthetabase*Nzetabase*Nxbase,0],Nxi16(inds),time16(inds),lb,ub,options);
      disp(['16 cores: minutes per Nxi   : ',num2str(kxi(1))])
    end
    inds=find((Ntheta16==Nthetabase) & (Nzeta16==Nzetabase) & (Nxi16==Nxibase));
    if length(inds)>=2
      kx=lsqcurvefit(fun,[k(1)*Nthetabase*Nzetabase*Nxbase,0],Nx16(inds),time16(inds),lb,ub,options);
      disp(['16 cores: minutes per Nx    : ',num2str(kx(1))])
    end
  else
    inds=find((NzetaNN==Nzetabase) & (NxiNN==Nxibase) & (NxNN==Nxbase));
    if length(inds)>=2
      ktheta=lsqcurvefit(fun,[k(1)*Nzetabase*Nxibase*Nxbase,0],NthetaNN(inds),timeNN(inds),lb,ub,options);
      disp(['NN cores: minutes per Ntheta: ',num2str(ktheta(1))])
    end
    inds=find((NthetaNN==Nthetabase) & (NxiNN==Nxibase) & (NxNN==Nxbase));
    if length(inds)>=2
      kzeta=lsqcurvefit(fun,[k(1)*Nthetabase*Nxibase*Nxbase,0],NzetaNN(inds),timeNN(inds),lb,ub,options);
      disp(['NN cores: minutes per Nzeta : ',num2str(kzeta(1))])
    end
    inds=find((NthetaNN==Nthetabase) & (NzetaNN==Nzetabase) & (NxNN==Nxbase));
    if length(inds)>=2
      kxi=lsqcurvefit(fun,[k(1)*Nthetabase*Nzetabase*Nxbase,0],NxiNN(inds),timeNN(inds),lb,ub,options);
      disp(['NN cores: minutes per Nxi   : ',num2str(kxi(1))])
    end
    inds=find((NthetaNN==Nthetabase) & (NzetaNN==Nzetabase) & (NxiNN==Nxibase));
    if length(inds)>=2
      kx=lsqcurvefit(fun,[k(1)*Nthetabase*Nzetabase*Nxbase,0],NxNN(inds),timeNN(inds),lb,ub,options);
      disp(['NN cores: minutes per Nx    : ',num2str(kx(1))])
    end
  end
end

fig(7)
if not(isempty(Ntot16)&isempty(Ntot36))
  plot(Ntot16,time16,'b+',Ntot36,time36,'r+')
  legend('16 proc data','36 proc',2)
else
  plot(NtotNN,timeNN,'b+')
  legend('data',2)
end
if curvefitworked
  hold on
  Ntotv=linspace(0,3e6,1000);
  plot(Ntotv,k(1)*Ntotv+k(2),'b-')
  hold off
end
xlabel('Nx * Nxi * Nzeta * Ntheta')
ylabel('time (min)')
title('simulation time')
