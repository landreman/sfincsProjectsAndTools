function [runs,missing]=getconvergencescanresults(dirpath,makeplots)
% For a simple convergence scan in Sfincs, 
% this loads the results stored in the subdirectories named 00/, 01/,... in "dirpath" 
% Plotting is made if makeplots=1.

[runs,missing]=getresults(dirpath);

if runs.NumElements==0
  warning('There were no successful runs !!')
else
  
convParams={'Ntheta','Nzeta','Nxi','Nx','solverTolerance','NL','NxPotentialsPerVth'};
isScanned=[];
vals=NaN*zeros(length(convParams),runs.NumElements);
for pind=1:length(convParams)
  vals(pind,:)=getfield(runs,convParams{pind});
  if all(vals(pind,:)==vals(pind,1))
    baseVals(pind)=vals(pind,1);
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
  didItConverge=runs.didItConverge;
  if didItConverge(baseRun)==1
    baseColor='r';
  else
    baseColor='y';
  end
  
  for pind=1:NParam
    
    Nscanruns=length(scanRuns{pind});
    colors=[];
    for j=1:Nscanruns
      colors=[colors,'b'];
    end
    colors(find(didItConverge(scanRuns{pind})~=1))='c';
    
    Lind=1; %L11
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    if not(isempty(baseRun))
      plot(baseVals(pind),L(baseRun,1,1),[baseColor,'+']);hold on;
      plot(baseVals(pind),LNTV(baseRun,1),[baseColor,'o'])
    end
    for srind=1:Nscanruns
      plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),1,1),...
           [colors(srind),'+'])    
      plot(scanVals{pind}(srind),LNTV(scanRuns{pind}(srind),1),...
           [colors(srind),'o'])
    end
    hold off
    xlabel(convParams{pind})
    ylabel('L11 +, LNTV1 o')
    %title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L12
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    if not(isempty(baseRun))
      plot(baseVals(pind),L(baseRun,1,2),[baseColor,'+']);hold on;
      plot(baseVals(pind),LNTV(baseRun,2),[baseColor,'o'])
    end
    for srind=1:Nscanruns
      plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),1,2),...
           [colors(srind),'+'])    
      plot(scanVals{pind}(srind),LNTV(scanRuns{pind}(srind),2),...
           [colors(srind),'o'])
    end
    hold off
    xlabel(convParams{pind})
    ylabel('L12 +, LNTV2 o')
    %title([param,' = ',num2str(P(vind).val)])
    
    Lind=Lind+1; %L22
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    if not(isempty(baseRun))
      plot(baseVals(pind),L(baseRun,2,2),[baseColor,'+']);hold on;
    end
    for srind=1:Nscanruns
      plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),2,2),...
           [colors(srind),'+'])    
    end
    hold off
    xlabel(convParams{pind})
    ylabel('L22')   
    
    Lind=Lind+1; %L13
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    if not(isempty(baseRun))
      plot(baseVals(pind),L(baseRun,1,3),[baseColor,'+']);hold on;
    end
    for srind=1:Nscanruns
      plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),1,3),...
           [colors(srind),'+'])    
    end
    hold off
    xlabel(convParams{pind})
    ylabel('L13')
    
    Lind=Lind+1; %L23
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    if not(isempty(baseRun))
      plot(baseVals(pind),L(baseRun,2,3),[baseColor,'+']);hold on;
    end
    for srind=1:Nscanruns
      plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),2,3),...
           [colors(srind),'+'])    
    end
    hold off
    % [min(scanval),max(scanval)],-4.57*/runs{1}.iota*[1,1],'k:')
    %Matt's value 4.57 could be wrong here
    xlabel(convParams{pind})
    ylabel('L23')
    
    Lind=Lind+1; %L33
    subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
    if not(isempty(baseRun))
      plot(baseVals(pind),L(baseRun,3,3),[baseColor,'+']);hold on;
    end
    for srind=1:Nscanruns
      plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),3,3),...
           [colors(srind),'+'])    
    end
    hold off
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
    if not(isempty(baseRun))
      plot(baseVals(pind),LNTV(baseRun,3),[baseColor,'+']);hold on;
    end
    for srind=1:Nscanruns
      plot(scanVals{pind}(srind),LNTV(scanRuns{pind}(srind),3),...
           [colors(srind),'+'])    
    end
    hold off
    %plot(baseVals(pind),LNTV(baseRun,3),'r+',scanVals{pind},LNTV(scanRuns{pind},3),'g+')
    xlabel(convParams{pind})
    ylabel('LNTV_3')
  end
end
end