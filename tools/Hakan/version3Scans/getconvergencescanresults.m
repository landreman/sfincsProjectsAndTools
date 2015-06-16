function [runs,missing]=getconvergencescanresults(dirpath,makeplots)
% For a simple convergence scan in Sfincs, 
% this loads the results stored in the subdirectories named 00/, 01/,... in "dirpath" 
% Normal plotting is made if makeplots=1.
% Plotting of values relative to baseCase is made if makeplots=2.
if nargin==1
  makeplots=1; %optional input
end


[runs,missing]=getresults(dirpath);

if not(isempty(missing))
  missing_runs={missing.message}
  if not(isempty(input('Runs are missing. Continue (return=yes, n=no) ?','s')))
    error('Stopped by user')
  end
end

if runs.NumElements==0
  warning('There were no successful runs !!')
else
  
if all(runs.RHSMode==3)
  convParams={'Ntheta','Nzeta','Nxi','solverTolerance','NL'};
else
  convParams={'Ntheta','Nzeta','Nxi','Nx','solverTolerance','NL','NxPotentialsPerVth'}; 
end

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
baseRun=find(prod(double(vals==(baseVals'*ones(1,Nruns)))));
if length(baseRun)>1
  baseRun=baseRun(1);
  warning('More than one base case runs were found!');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%makeplots=1;
if makeplots
  close all
  %For each value of 'param', make different scan variable plots

  if all(runs.RHSMode==3) %Mono-energetic
    D=runs.transportCoeffs;
    Nrow=3;
    Ncol=NParam;
    stringForTop=sprintf(...
      'r/a = %4.2f,  nuPrime = %7.2e,  EStar = %7.2e',...
      runs.rN(1),runs.nuPrime(1),...
      runs.EStar(1));
  elseif all(runs.RHSMode==2)
    L=runs.transportMatrix;
    NTVisstored=0;
    if isfield(runs,'NTVMatrix')
      NTVisstored=1;
      LNTV=runs.NTVMatrix;
    else
      LNTV=NaN*ones(size(L,1),3);
    end
    Nrow=NParam;
    Ncol=6+1;
    stringForTop=sprintf(...
      'r/a = %4.2f,  nuPrime = %7.2e,  EStar = %7.2e',...
      runs.rN(1),runs.nuPrime(1),...
      runs.EStar(1));
  else
    Nrow=NParam;
    Ncol=runs.Nspecies(1)*4+1;
    stringForTop=sprintf(...
      'r/a = %4.2f, dPhiHatdpsiN = %7.2e',...
      runs.rN(1),...
      runs.dPhiHatdpsiN(1));
  end
  finished=runs.finished;
  if finished(baseRun)==1
    baseColor='r';
  else
    baseColor='y';
  end
  fig(60)
  pos=get(gcf,'Position');
  pos(3:4)=[200*Ncol,200*Nrow];
  set(gcf,'Position',pos);
  
  if all(runs.RHSMode==3) %monoenergetic
    for pind=1:NParam

      Nscanruns=length(scanRuns{pind});
      colors=[];
      for j=1:Nscanruns
        colors=[colors,'b'];
      end
      colors(find(finished(scanRuns{pind})~=1))='c';
      
      Lind=1; %D11
      subplot(Nrow,Ncol,(Lind-1)*Ncol+pind)
      if makeplots==1
        if not(isempty(baseRun))
          plot(baseVals(pind),D(baseRun,1,1),[baseColor,'+']);hold on
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),D(scanRuns{pind}(srind),1,1),...
               [colors(srind),'+']);hold on    
        end
      elseif makeplots==2
        if not(isempty(baseRun))
          plot(baseVals(pind),1,[baseColor,'+']);hold on
        else
          error('No base run to compare with! (makeplots=2)')
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),D(scanRuns{pind}(srind),1,1)/D(baseRun,1,1),...
               [colors(srind),'+']);hold on    
        end        
      end
      hold off
      xlabel(convParams{pind})
      ylabel('D11 +')
      %title([param,' = ',num2str(P(vind).val)])
      
      Lind=Lind+1; %D13
      subplot(Nrow,Ncol,(Lind-1)*Ncol+pind)
      if makeplots==1
        if not(isempty(baseRun))
          plot(baseVals(pind),D(baseRun,1,2),[baseColor,'+']);hold on;
          plot(baseVals(pind),D(baseRun,2,1),[baseColor,'o'])
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),D(scanRuns{pind}(srind),1,2),...
               [colors(srind),'+']);hold on    
          plot(scanVals{pind}(srind),D(scanRuns{pind}(srind),2,1),...
               [colors(srind),'o'])    
        end
      elseif makeplots==2
        if not(isempty(baseRun))
          plot(baseVals(pind),1,[baseColor,'+']);hold on;
          plot(baseVals(pind),D(baseRun,2,1)/D(baseRun,1,2),[baseColor,'o'])
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),D(scanRuns{pind}(srind),1,2)/D(baseRun,1,2),...
               [colors(srind),'+']);hold on    
          plot(scanVals{pind}(srind),D(scanRuns{pind}(srind),2,1)/D(baseRun,1,2),...
               [colors(srind),'o'])    
        end
      end
      hold off
      xlabel(convParams{pind})
      ylabel('D13 +, D31 o')
      %title([param,' = ',num2str(P(vind).val)])
      
      Lind=Lind+1; %D33
      subplot(Nrow,Ncol,(Lind-1)*Ncol+pind)
      if makeplots==1
        if not(isempty(baseRun))
          plot(baseVals(pind),D(baseRun,2,2),[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),D(scanRuns{pind}(srind),2,2),...
               [colors(srind),'+']);hold on    
        end
      elseif makeplots==2
        if not(isempty(baseRun))
          plot(baseVals(pind),1,[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),D(scanRuns{pind}(srind),2,2)/D(baseRun,2,2),...
               [colors(srind),'+']);hold on    
        end
      end
      hold off
      xlabel(convParams{pind})
      ylabel('D33')   
      
    end
    
    
  elseif all(runs.RHSMode==2)
    for pind=1:NParam
      
      Nscanruns=length(scanRuns{pind});
      colors=[];
      for j=1:Nscanruns
        colors=[colors,'b'];
      end
      colors(find(finished(scanRuns{pind})~=1))='c';
      
      Lind=1; %L11
      subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
      if makeplots==1
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
      elseif makeplots==2
        if not(isempty(baseRun))
          plot(baseVals(pind),1,[baseColor,'+']);hold on;
          plot(baseVals(pind),LNTV(baseRun,1)/L(baseRun,1,1),[baseColor,'o'])
        else
          error('No base run to compare with! (makeplots=2)')
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),1,1)/L(baseRun,1,1),...
               [colors(srind),'+'])    
          plot(scanVals{pind}(srind),LNTV(scanRuns{pind}(srind),1)/L(baseRun,1,1),...
               [colors(srind),'o'])
        end
      end
      hold off
      xlabel(convParams{pind})
      ylabel('L11 +, LNTV1 o')
      %title([param,' = ',num2str(P(vind).val)])
      
      Lind=Lind+1; %L12
      subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
      if makeplots==1
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
      elseif makeplots==2
        if not(isempty(baseRun))
          plot(baseVals(pind),1,[baseColor,'+']);hold on;
          plot(baseVals(pind),LNTV(baseRun,2)/L(baseRun,1,2),[baseColor,'o'])
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),1,2)/L(baseRun,1,2),...
               [colors(srind),'+'])    
          plot(scanVals{pind}(srind),LNTV(scanRuns{pind}(srind),2)/L(baseRun,1,2),...
               [colors(srind),'o'])
        end        
      end
      hold off
      xlabel(convParams{pind})
      ylabel('L12 +, LNTV2 o')
      %title([param,' = ',num2str(P(vind).val)])
      
      Lind=Lind+1; %L22
      subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
      if makeplots==1
        if not(isempty(baseRun))
          plot(baseVals(pind),L(baseRun,2,2),[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),2,2),...
               [colors(srind),'+'])    
        end
      elseif makeplots==2
        if not(isempty(baseRun))
          plot(baseVals(pind),1,[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),2,2)/L(baseRun,2,2),...
               [colors(srind),'+'])    
        end      
      end
      hold off
      xlabel(convParams{pind})
      ylabel('L22')   
      
      Lind=Lind+1; %L13
      subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
      if makeplots==1
        if not(isempty(baseRun))
          plot(baseVals(pind),L(baseRun,1,3),[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),1,3),...
               [colors(srind),'+'])    
        end
      elseif makeplots==2
        if not(isempty(baseRun))
          plot(baseVals(pind),1,[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),1,3)/L(baseRun,1,3),...
               [colors(srind),'+'])    
        end
      end
      hold off
      xlabel(convParams{pind})
      ylabel('L13')
      
      Lind=Lind+1; %L23
      subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
      if makeplots==1
        if not(isempty(baseRun))
          plot(baseVals(pind),L(baseRun,2,3),[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),2,3),...
             [colors(srind),'+'])    
        end
      elseif makeplots==2
        if not(isempty(baseRun))
          plot(baseVals(pind),1,[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),2,3)/L(baseRun,2,3),...
             [colors(srind),'+'])    
        end
      end
      hold off
      hold off
      % [min(scanval),max(scanval)],-4.57*/runs{1}.iota*[1,1],'k:')
      %Matt's value 4.57 could be wrong here
      xlabel(convParams{pind})
      ylabel('L23')
      
      Lind=Lind+1; %L33
      subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
      if makeplots==1
        if not(isempty(baseRun))
          plot(baseVals(pind),L(baseRun,3,3),[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),3,3),...
               [colors(srind),'+'])    
        end
      elseif makeplots==2
        if not(isempty(baseRun))
          plot(baseVals(pind),1,[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),L(scanRuns{pind}(srind),3,3)/L(baseRun,3,3),...
               [colors(srind),'+'])    
        end
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
      if makeplots==1
        if not(isempty(baseRun))
          plot(baseVals(pind),LNTV(baseRun,3),[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),LNTV(scanRuns{pind}(srind),3),...
               [colors(srind),'+'])    
        end
      elseif makeplots==2
        if not(isempty(baseRun))
          plot(baseVals(pind),1,[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),LNTV(scanRuns{pind}(srind),3)/LNTV(baseRun,3),...
               [colors(srind),'+'])    
        end
      end
      hold off
      %plot(baseVals(pind),LNTV(baseRun,3),'r+',scanVals{pind},LNTV(scanRuns{pind},3),'g+')
      xlabel(convParams{pind})
      ylabel('LNTV_3')
    end
    
  else %%%%%%%% RHSMode=1
    
    for pind=1:NParam
      
      Nscanruns=length(scanRuns{pind});
      colors=[];
      for j=1:Nscanruns
        colors=[colors,'b'];
      end
      colors(find(finished(scanRuns{pind})~=1))='c';
      Lind=0;

      for spec=1:runs.Nspecies(1)
        Lind=Lind+1; %particleFlux
        subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
        
        if 0
          e=1.6022e-19;
          mp=1.6726e-27;
          vbar=sqrt(2*1e3*e/mp);
          a=0.6258; %%%%%% NOTE MANUAL INPUT
          showthis=runs.particleFlux_vm_psiN*a/2/runs.rN(1)*vbar*1e20; str='particleFlux /m^2/s';
        end
        
        %showthis=runs.NTVfromFlux; str='NTVfromFlux';
        showthis=runs.particleFlux_vm_psiN; str='particleFlux';
        
        if makeplots==1
          if not(isempty(baseRun))
            plot(baseVals(pind),showthis(baseRun,spec),[baseColor,'+']);hold on;
          end
          for srind=1:Nscanruns
            plot(scanVals{pind}(srind),showthis(scanRuns{pind}(srind),spec),...
                 [colors(srind),'+'])    
          end
        elseif makeplots==2
          if not(isempty(baseRun))
            plot(baseVals(pind),1,[baseColor,'+']);hold on;
          end
          for srind=1:Nscanruns
            plot(scanVals{pind}(srind),...
            showthis(scanRuns{pind}(srind),spec)/...
            showthis(baseRun,spec) ,...
                 [colors(srind),'+'])    
          end
        end
        xlabel(convParams{pind})
        ylabel([str,' spec. ',num2str(spec)])

        Lind=Lind+1;
        subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
        if makeplots==1
          if not(isempty(baseRun))
            plot(baseVals(pind),runs.heatFlux_vm_psiN(baseRun,spec),[baseColor,'+']);hold on;
          end
          for srind=1:Nscanruns
            plot(scanVals{pind}(srind),runs.heatFlux_vm_psiN(scanRuns{pind}(srind),spec),...
                 [colors(srind),'+'])    
          end
        elseif makeplots==2
          if not(isempty(baseRun))
            plot(baseVals(pind),1,[baseColor,'+']);hold on;
          end
          for srind=1:Nscanruns
            plot(scanVals{pind}(srind),...
                 runs.heatFlux_vm_psiN(scanRuns{pind}(srind),spec)/...
                 runs.heatFlux_vm_psiN(baseRun,spec),...
                 [colors(srind),'+'])    
          end
        end
        xlabel(convParams{pind})
        ylabel(['heatFlux spec. ',num2str(spec)])
  
        Lind=Lind+1;
        subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
        if makeplots==1
          if not(isempty(baseRun))
            plot(baseVals(pind),runs.FSABFlow(baseRun,spec),[baseColor,'+']);hold on;
          end
          for srind=1:Nscanruns
            plot(scanVals{pind}(srind),runs.FSABFlow(scanRuns{pind}(srind),spec),...
                 [colors(srind),'+'])    
          end
        elseif makeplots==2
          if not(isempty(baseRun))
            plot(baseVals(pind),1,[baseColor,'+']);hold on;
          end
          for srind=1:Nscanruns
            plot(scanVals{pind}(srind),...
                 runs.FSABFlow(scanRuns{pind}(srind),spec)/...
                 runs.FSABFlow(baseRun,spec),...
                 [colors(srind),'+'])    
          end
        end
        xlabel(convParams{pind})
        ylabel(['FSABFlow spec. ',num2str(spec)])

        Lind=Lind+1;
        subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
        if makeplots==1
          if not(isempty(baseRun))
            plot(baseVals(pind),runs.NTV(baseRun,spec),[baseColor,'+']);hold on;
          end
          for srind=1:Nscanruns
            plot(scanVals{pind}(srind),runs.NTV(scanRuns{pind}(srind),spec),...
                 [colors(srind),'+'])    
          end
        elseif makeplots==2
          if not(isempty(baseRun))
            plot(baseVals(pind),1,[baseColor,'+']);hold on;
          end
          for srind=1:Nscanruns
            plot(scanVals{pind}(srind),...
                 runs.NTV(scanRuns{pind}(srind),spec)/...
                 runs.NTV(baseRun,spec),...
                 [colors(srind),'+'])    
          end
        end
        xlabel(convParams{pind})
        ylabel(['NTV spec. ',num2str(spec)])

      end
      
      %factor=1.6022e-19*1e3*1e20; %conversion to SI
      factor=1;
      
      Lind=Lind+1;
      subplot(Nrow,Ncol,(pind-1)*Ncol+Lind)
      if makeplots==1
        if not(isempty(baseRun))
          plot(baseVals(pind),sum(runs.NTV(baseRun,:)),[baseColor,'+']);hold on;
          plot(baseVals(pind),sum(runs.NTVfromFlux(baseRun,:)),[baseColor,'o']);
        end        
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),...
               factor*sum(runs.NTV(scanRuns{pind}(srind),:)),[colors(srind),'+'])
          plot(scanVals{pind}(srind),...
               factor*sum(runs.NTVfromFlux(scanRuns{pind}(srind),:)),[colors(srind),'o'])    
          %sr=sum(runs.NTVfromFlux(scanRuns{pind}(srind),:))./sum(runs.NTV(scanRuns{pind}(srind),:))
        end
      elseif makeplots==2
        if not(isempty(baseRun))
          plot(baseVals(pind),1,[baseColor,'+']);hold on;
        end
        for srind=1:Nscanruns
          plot(scanVals{pind}(srind),...
               sum(runs.NTV(scanRuns{pind}(srind),:))/...
               sum(runs.NTV(baseRun,:)),...
               [colors(srind),'+'])    
        end
      end
      xlabel(convParams{pind})
      ylabel('NTV tot')

    end
  end
  annotation('textbox',[0 0.96 1 .07],...
             'HorizontalAlignment','center',...
            'Interpreter','none',...
             'VerticalAlignment','bottom',...
            'FontSize',12,'LineStyle','none',...
             'String',stringForTop);
end
end