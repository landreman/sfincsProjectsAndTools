function [runs,missing]=getconvergencescanresults(dirpath,makeplots)
% For a simple convergence scan in Sfincs, 
% this loads the results stored in the subdirectories named 00/, 01/,... in "dirpath" 
% Plotting is made if makeplots=1.

[runs,missing]=getresults(dirpath);

convParams={'Ntheta','Nzeta','Nxi','Nx','solverTolerance','NL','NxPotentialsPerVth'};
isScanned=[];
for pind=1:length(convParams)
  vals(pind,:)=getfield(runs,convParams{pind});
  if all(vals(pind,:)==vals(pind,1))
    baseVals(pind)=vals(1,pind);
    scanVals{pind}=[];
    scanRuns{pind}=[];
  else
    isScanned=[isScanned,pind];
    sorted=sort(vals(pind,:));
    newind=[1,find(diff(sorted))+1];
    [dummy,maxind]=max(diff([newind,length(sorted)+1]));
    baseVals(pind)=sorted(newind(maxind));
    scanVals{pind}=sorted([newind(1:maxind-1),newind(maxind+1:end)]);
    for ind=1:length(scanVals{pind})
      rind=find(vals(pind,:)==scanVals{pind}(ind));
      if length(rind)>1
        warning('length(rind)>1')
      end
      scanRuns{pind}(ind)=rind(1);
    end
  end
end
NParam=length(isScanned);
Nruns=size(vals,2);
convParams={convParams{isScanned}};
baseVals=baseVals(isScanned);
scanVals={scanVals{isScanned}};
scanRuns={scanRuns{isScanned}};
vals=vals(isScanned,:);
baseRun=find(prod(vals==(baseVals'*ones(1,Nruns))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%makeplots=1;
if makeplots
  close all
  %For each value of 'param', make different scan variable plots

  Nrow=NParam;
  Ncol=6+1;
  fig(60)
  pos=get(gcf,'Position');
  pos(3:4)=[200*Ncol,200*Nrow];
  set(gcf,'Position',pos);

  L=runs.transportMatrix;
  LNTV=runs.NTVMatrix;

  for pind=1:NParam
    
    Lind=1; %L11
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    plot(baseVals(pind),L(baseRun,1,1),'r+',scanVals{pind},L(scanRuns{pind},1,1),'g+',...
         baseVals(pind),LNTV(baseRun,1),'ro',scanVals{pind},LNTV(scanRuns{pind},1),'go')
    xlabel(convParams{pind})
    ylabel('L11 +, LNTV1 o')
    %title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L12
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    plot(baseVals(pind),L(baseRun,1,2),'r+',scanVals{pind},L(scanRuns{pind},1,2),'g+',...
         baseVals(pind),LNTV(baseRun,2),'ro',scanVals{pind},LNTV(scanRuns{pind},2),'go')
    xlabel(convParams{pind})
    ylabel('L12 +, LNTV2 o')
    %title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L22
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    plot(baseVals(pind),L(baseRun,2,2),'r+',scanVals{pind},L(scanRuns{pind},2,2),'g+')%,...
                                                                                      %baseVals(pind),LNTV(baseRun,3),'ro',scanVals{pind},LNTV(scanRuns{pind},3),'go')
    xlabel(convParams{pind})
    ylabel('L22')% +, LNTV3 o')
                 %title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L13
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    plot(baseVals(pind),L(baseRun,1,3),'r+',scanVals{pind},L(scanRuns{pind},1,3),'g+')%,...
                                                                                      %[min(scanval)*0.97,max(scanval)*1.03],-1/runs{1}.iota(1)*[1,1],'k:')
    xlabel(convParams{pind})
    ylabel('L13')
    %title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L23
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    plot(baseVals(pind),L(baseRun,2,3),'r+',scanVals{pind},L(scanRuns{pind},2,3),'g+')%,...
                                                                                      % [min(scanval),max(scanval)],-4.57*/runs{1}.iota*[1,1],'k:')
                                                                                      %Matt's value 4.57 could be wrong here
    xlabel(convParams{pind})
    ylabel('L23')
    %title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L33
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    plot(baseVals(pind),L(baseRun,3,3),'r+',scanVals{pind},L(scanRuns{pind},3,3),'g+')
    xlabel(convParams{pind})
    ylabel('L33')
    %title([param,' = ',num2str(P(vind).val)])
    
    %Lind=Lind+1;
    %subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    %plot(baseVals(pind),LNTV(baseRun,1),'r+',scanVals{pind},LNTV(scanRuns{pind},1),'g+')
    %xlabel(convParams{pind})
    %ylabel('LNTV_1')
    %title([param,' = ',num2str(P(vind).val)])
    
    %Lind=Lind+1;
    %subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    %plot(baseVals(pind),LNTV(baseRun,2),'r+',scanVals{pind},LNTV(scanRuns{pind},2),'g+')
    %xlabel(convParams{pind})
    %ylabel('LNTV_2')
    %title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1;
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    plot(baseVals(pind),LNTV(baseRun,3),'r+',scanVals{pind},LNTV(scanRuns{pind},3),'g+')
    xlabel(convParams{pind})
    ylabel('LNTV_3')
    %title([param,' = ',num2str()])

  end
end