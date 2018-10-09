function d11fits=dkes_ftd11(Geomdat,dkdata,varargin);

if nargin>=3
  fignr=varargin{1};
else
  fignr=1;
end
if nargin==4
  dkesoutputpath=varargin{2};
  if not(exist([dkesoutputpath,'ps'],'dir'))
    success=mkdir(dkesoutputpath,'ps');
    if not(success)
      error(['Could not create the directory ',dkesoutputpath,'ps'])
    end
  end
else
  dkesoutputpath='';
end
if ishandle(fignr)
  close(fignr)
end

for rind=1:dkdata.Nradii
  %fignr=1+mod(floor((rind-1)/2),4)+((-1)^rind+1)/2*4;
  %fignr=rind;
  
  %Check if it is a tokamak first
  if all(Geomdat.n{rind}==0)
      d11fits.eps_efh(rind)=-1;%NaN;
      d11fits.g11_ft(rind)=-1;%NaN;
      d11fits.g11_er(rind)=-1;%NaN;
      d11fits.ex_er(rind)=-1;%NaN;
      d11fits.er_u(rind)=-1;%NaN;    
  else

    nodk=dkdata.Nruns(rind); %number of dkes data
    
    %ind11=find(Geomdat.m{rind}==1 & Geomdat.n{rind}==1);
    %ind01=find(Geomdat.m{rind}==0 & Geomdat.n{rind}==1);
    ind10=find(Geomdat.m{rind}==1 & Geomdat.n{rind}==0);
    %b_11=Geomdat.Bnorm{rind}(ind11);
    %b_01=Geomdat.Bnorm{rind}(ind01);
    b_10=Geomdat.Bnorm{rind}(ind10);

    b0dk=1.0;
    
    R=Geomdat.R00(rind);
    epsilon=Geomdat.rnorm(rind)*Geomdat.minorradiusW7AS/R;
    xkn=abs(b_10/epsilon);
    aiota=abs(Geomdat.iota(rind));
    
    %del=1e-3;
    %x_l=1e10;
    %x_u=-1e10;
    %y1_l=1e10;
    %y1_u=-1e10;
    
    % sorting the data for efield = 0
    
    [~,E0inds]=find(dkdata.efld{rind}==0);
    if length(E0inds)>0
      [cmul_E0,srti]=sort(dkdata.cmul{rind}(E0inds),'descend'); %HM: anu
      d11_E0 =abs(dkdata.d11{rind}(E0inds(srti)));  %HM: a110  
      d11e_E0=dkdata.d11e{rind}(E0inds(srti)); %HM: a110e    
    else
      warning('There were no runs for Er = 0 !')
      disp('Using the lowest Er for extrapolation!')
      cmul_E0=dkdata.cmulfirst{rind}.cmul;
      d11_E0=[];
      d11e_E0=[];
      for cmulind=1:dkdata.cmulfirst{rind}.Ncmuls
        [~,lowEind]=min(abs(dkdata.cmulfirst{rind}.efld{cmulind}));
        d11_E0=[d11_E0,dkdata.cmulfirst{rind}.d11{cmulind}(lowEind)];
        d11e_E0=[d11e_E0,dkdata.cmulfirst{rind}.d11e{cmulind}(lowEind)];
      end
    end
    %x_l=min(cmul_E0);
    %x_u=max(cmul_E0);
    %y1_l=min(abs(d11_E0));
    %y1_u=max(abs(d11_E0));
    
    %scaling to the analytical plateau value
    
    fnrm11=pi/8/(R*aiota*b0dk^2);
    fnrmnu=R/aiota;
    fnrmEr=1/epsilon/aiota/b0dk;
    
    % convergence plot for D_11 in the 1/nu regime
    
    %xs = log10(max(1e-9,cmul_E0));
    g11=d11_E0.*cmul_E0;
    g11_l=(d11_E0-d11e_E0).*cmul_E0;
    g11_u=(d11_E0+d11e_E0).*cmul_E0;
    %eps_h=(4.9982*g11*R^2*b0dk^2).^(2/3);
    %xp=xs;
    %yp=eps_h;
    %ype(1,:)=xs;
    %ype(2,:)=xs;
    %ype(3,:)=(4.9982*g11_l*R^2*b0dk^2).^(2/3);
    %ype(4,:)=(4.9982*g11_u*R^2*b0dk^2).^(2/3);
    
    %y4g_l=min(1e5,g11_u);
    %y4e_l=min(ype(4,:));
    
    if length(cmul_E0)<2
      disp(['Warning: There was < 2 runs for Er = 0 at radius #',num2str(rind)])
      makeafit=0;
      acceptedfit=1;
    else
      ask_for_exp_fit=0;
      if ask_for_exp_fit
        disp('Please give the exp in the fit of g11 = cmul*D_11 with y0+y1*cmul^exp')
        exp_fit=input('exp=? (default is 0.5):');
        if isempty(exp_fit)
          exp_fit=0.5;
        end
      else
        exp_fit=0.5;
      end

      % try to estimate how many points that are needed in the fit 
      % by finding the inflection point
      Nconsider=min(9,length(cmul_E0)-1);
      cmulexp_10=cmul_E0(end-Nconsider+1:end).^exp_fit;
      g11_10=g11(end-Nconsider+1:end);
      dfdx=diff(g11_10)./diff(cmulexp_10);
      x1=(cmulexp_10(1:end-1)+cmulexp_10(2:end))/2;
      d2fdx2=diff(dfdx)./diff(x1);
      x2=(x1(1:end-1)+x1(2:end))/2;
      inflectionpoint=interp1(d2fdx2,x2,0);
      if not(isempty(inflectionpoint))
        indx=find(diff(sign(cmulexp_10-inflectionpoint(1))));
        Ndataforfit=max(4,Nconsider-indx(1)+1);
      else
        Ndataforfit=5;
      end
      
      fprintf(1,'\n')
      disp('----------------------------------------------')
      disp(['Radius ',num2str(rind),'/',num2str(dkdata.Nradii),...
            ': r=',num2str(Geomdat.minorradiusW7AS*Geomdat.rnorm(rind)),...
            ' m, Figure nr ',num2str(fignr)])
      disp('----------------------------------------------')
      fprintf(1,'\n')
      
      %set(0, 'DefaultFigureVisible', 'off');
      fig(fignr)
      makeafit=1;
      acceptedfit=0;
      
      iters=1;
      while not(acceptedfit)
        close(fignr)
        fig(fignr)    
        pos=get(gcf,'Position');
        pos(3)=840;
        set(gcf,'Position',pos);
        
        while 1
          fitinds=length(g11)-Ndataforfit+1:length(g11);
          g11_part=g11(fitinds);              %HM: awk(:,1)
          cmulexp_part=cmul_E0(fitinds).^exp_fit; %HM: awk(:,2)
          weights = min(1e4,(d11e_E0(fitinds)./d11_E0(fitinds)).^(-2));
          
          F0=fit(cmulexp_part',g11_part','poly1','weight',weights);
          y1=F0.p1;
          y0=F0.p2; %HM: g11_ft ???
          g11fit=F0.p1*cmul_E0.^exp_fit+F0.p2;
          %g11_part_fit=F0.p1*cmulexp_part+F0.p2;
          eps_efh=(4.9982*y0*R^2*b0dk^2).^(2/3);
          if iters<2 && g11(fitinds(1)-1)<g11fit(fitinds(1)-1)
            Ndataforfit=Ndataforfit+1;
            iters=iters+1; 
          else
            break
          end
        end
        
        subplot(121)    
        showinds=max(1,length(g11)-Ndataforfit-1):length(g11);
        %       cmul_E0(fitinds).^exp_fit,g11(fitinds),'rx',...
        hold off
        plot(cmul_E0(showinds).^exp_fit,g11fit(showinds),'r-',...
             cmul_E0(showinds).^exp_fit,g11(showinds),'rs',...
             cmul_E0(fitinds).^exp_fit,g11(fitinds),'r*',...
             cmul_E0(fitinds).^exp_fit,g11(fitinds),'rd',...
             cmul_E0(fitinds).^exp_fit,g11(fitinds),'r+',...
             cmul_E0(fitinds).^exp_fit,g11_l(fitinds),'r^',...
             cmul_E0(fitinds).^exp_fit,g11_u(fitinds),'rv')%,...
                                                           %     inflectionpoint*[1,1],[min(g11(showinds)),max(g11(showinds))],'r:')
                                                           %xlim([cmul_E0(end).^exp_fit,cmul_E0(end-Ndataforfit-1).^exp_fit])
                                                           %text(inflectionpoint,max(g11(showinds)),'inflection point');
        xlabel(['(\nu/v)^{',num2str(exp_fit),'}'])
        ylabel('g_{11} = (\nu/v) * |\Gamma_{11}|')
        title(['r=',num2str(100*Geomdat.minorradiusW7AS*Geomdat.rnorm(rind)),...
               'cm: ',num2str(rind),'/',num2str(dkdata.Nradii),...
               ', \epsilon_{eff}=',num2str(eps_efh,'%5.5f')])%,...
                                                             %', g_{11,ft}=',num2str(y0,'%5.2e')])
        legend([num2str(y0,'%5.2e'),'+y_1(\nu/v)^{',num2str(exp_fit),'}'],2)
        drawnow
        %set(0, 'DefaultFigureVisible', 'on');
        %figHandles = findall(0, 'Type', 'figure');
        %set(figHandles(:), 'visible', 'on')
        
        disp(['fit of g11 = cmul*D_11 with y0+y1*cmul^',num2str(exp_fit)])
        disp(['up to cmul =',num2str(cmul_E0(fitinds(1)))])
        disp(['--> effective helical ripple: ',...
              'eps_efh=(4.9982*y0*R^2*B0^2)^(2/3)= ',num2str(eps_efh)])
        if y0<0
          disp(' WARNING, y0 < 0! Choose another number of points!')
        end
        newNdataforfit=...
            input('If you accept the fit press return, otherwise give the number of points (0=make no fit):');
        if isempty(newNdataforfit)
          acceptedfit=1;
        elseif newNdataforfit==0
          acceptedfit=1;
          makeafit=0;
        else
          Ndataforfit=newNdataforfit;
        end
      end  
    end  
    
    if not(makeafit)
      d11fits.eps_efh(rind)=NaN;
      d11fits.g11_ft(rind)=NaN;
      d11fits.g11_er(rind)=NaN;
      d11fits.ex_er(rind)=NaN;
      d11fits.er_u(rind)=NaN;
    else
      d11fits.eps_efh(rind)=eps_efh;
      d11fits.g11_ft(rind)=y0;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%% Fit of sqrt(nu) regime %%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      %First make a plot so that the used can decide a range of cmuls and efields

      figure(fignr)
      subplot(122)
      colorsandstyle={'ro','go','bo','co','mo','ko',...
                      'r+','g+','b+','c+','m+','k+',...
                      'rx','gx','bx','cx','mx','kx',...
                      'rd','gd','bd','cd','md','kd'};
      disp('List of Efld:')
      disp('---------------------')
      for efldind=1:dkdata.efldfirst{rind}.Neflds
        theefld=dkdata.efldfirst{rind}.efld(efldind);
        loglog(dkdata.efldfirst{rind}.cmul{efldind},...
               abs(dkdata.efldfirst{rind}.d11{efldind}),...
               colorsandstyle{mod(efldind-1,length(colorsandstyle))+1})
        hold on
        xlabel('\nu/v');
        ylabel('|\Gamma_{11}|')
        disp(['Efld: ',num2str(theefld),', symbol: ',...
              colorsandstyle{mod(efldind-1,length(colorsandstyle))+1}])
      end
      hold off
      disp('---------------------')
      
      %if Ndataforfit<=2
      %  disp('Too little data to make a sqrt(nu) regime fit!')
      %  disp('----------------------------------------------')
      %  d11fits.g11_er(rind)=0;
      %  d11fits.ex_er(rind)=0;
      %  d11fits.er_u(rind)=0;
      %else
      
      
      cmul_ftl_def=cmul_E0(end); %take lowest cmul_E0 as default lower limit
      cmul_ftu_def=min(1e-4,cmul_E0(max(1,length(cmul_E0)-Ndataforfit))); %default uppr limit
      er_res=aiota*epsilon;
      er_u_def=0.05*er_res;
      er_l_def=max(3e-7,cmul_ftl_def);

      do_input=1;
      acceptedfit=0;
      iters=0;
      while not(acceptedfit)
        iters=iters+1;
        
        if not(do_input)
          cmul_ftl=[];cmul_ftu=[];er_l=[];er_u=[];      
        else
          if iters == 1 %use defaults the first time
            cmul_ftl=cmul_ftl_def;
            cmul_ftu=cmul_ftu_def;
            er_l=er_l_def;
            er_u=er_u_def;
          else
            problemreading=1;
            while problemreading
              cmul_ftl_str=input(['lower cmul limit for E_r fit (deafult=',...
                                  num2str(cmul_ftl_def),'): '],'s');
              if strcmp(cmul_ftl_str,'')
                problemreading=0;
                cmul_ftl=cmul_ftl_def;
              else
                cmul_ftl=str2num(cmul_ftl_str);
                if isempty(cmul_ftl)
                  disp('INPUT NOT UNDERSTANDABLE')
                else
                  problemreading=0;
                end
              end
            end
            
            problemreading=1;
            while problemreading
              cmul_ftu_str=input(['upper cmul limit for E_r fit (deafult=',...
                                  num2str(cmul_ftu_def),'): '],'s');
              if strcmp(cmul_ftu_str,'')
                problemreading=0;
                cmul_ftu=cmul_ftu_def;
              else
                cmul_ftu=str2num(cmul_ftu_str);
                if isempty(cmul_ftu)
                  disp('INPUT NOT UNDERSTANDABLE')
                else
                  problemreading=0;
                end
              end
            end
            
            problemreading=1;
            while problemreading
              er_l_str=input(['lower efield limit for E_r fit (deafult=',...
                              num2str(er_l_def),'): '],'s');
              if strcmp(er_l_str,'')
                problemreading=0;
                er_l=er_l_def;
              else
                er_l=str2num(er_l_str);
                if isempty(er_l)
                  Ercolorind=find(strcmp(colorsandstyle,er_l_str));
                  if isempty(Ercolorind)
                    disp('INPUT NOT UNDERSTANDABLE')
                  else
                    problemreading=0;
                    er_l=dkdata.efldfirst{rind}.efld(Ercolorind);
                  end
                else
                  problemreading=0;
                end
              end
            end
            
            problemreading=1;
            while problemreading
              er_u_str=input(['upper efield limit for E_r fit (deafult=',...
                              num2str(er_u_def),'): '],'s');
              if strcmp(er_u_str,'')
                problemreading=0;
                er_u=er_u_def;
              else
                er_u=str2num(er_u_str);
                if isempty(er_u)
                  if length(er_u_str)==1
                    er_u_str=[er_u_str,'o'];
                  end
                  Ercolorind=find(strcmp(colorsandstyle,er_u_str));
                  if isempty(Ercolorind)
                    disp('INPUT NOT UNDERSTANDABLE')
                  else
                    problemreading=0;
                    er_u=dkdata.efldfirst{rind}.efld(Ercolorind);
                  end
                else
                  problemreading=0;
                end
              end
            end
          end
        end

        good.ind=find(dkdata.efld{rind}>=er_l & ...
                      dkdata.efld{rind}<=er_u & ...
                      dkdata.cmul{rind}>=cmul_ftl & ...
                      dkdata.cmul{rind}<=cmul_ftu);
        good.efld=dkdata.efld{rind}(good.ind);
        good.cmul=dkdata.cmul{rind}(good.ind);
        good.d11=abs(dkdata.d11{rind}(good.ind));
        good.d11e=dkdata.d11e{rind}(good.ind);
        good.lgd11u=log(good.d11+good.d11e);
        good.lgd11l=log(good.d11-good.d11e);
        good.rel_err=good.d11e./good.d11;
        weights=1./good.rel_err;
        
        Nex_er=150;
        Ng11_er=150;
        ex_ers=linspace(0.5,8,Nex_er);
        g11_ers=logspace(-10,-2,Ng11_er);
        [ex_erM,g11_erM]=ndgrid(ex_ers,g11_ers);
        cost=zeros(size(ex_erM));
        for ex_er_ind=1:Nex_er
          for g11_er_ind=1:Ng11_er
            ex_er=ex_ers(ex_er_ind);
            g11_er=g11_ers(g11_er_ind);
            d11appr_sqrtnu=sqrt(good.cmul)*g11_er./good.efld.^(3/2);
            d11appr_nuinv=y0./good.cmul;
            lgd11appr=log((d11appr_nuinv.^(-ex_er)+d11appr_sqrtnu.^(-ex_er)).^(-1/ex_er));
            %cost(ex_er_ind,g11_er_ind)=sum((log(d11appr)-log(good.d11)).^2.*weights)./ ...
            %    sum(weights);
            cost(ex_er_ind,g11_er_ind)=...
                sum( (lgd11appr-good.lgd11u).^2.*(lgd11appr>good.lgd11u) +...
                     (good.lgd11l-lgd11appr).^2.*(lgd11appr<good.lgd11l));
          end
        end
        [~,totind]=min(cost(:));
        tmp=ex_erM(:);
        ex_er=tmp(totind);
        tmp=g11_erM(:);
        g11_er=tmp(totind);
        
        displayit=0;
        if displayit
          good_efld=good.efld
          good_cmul=good.cmul
          good_d11=good.d11
          good_d11e=good.d11e
          ex_er
          g11_er
          cmul_ftu_def
          error('stopppp')
        end
        
        %g11_er=sum(good.d11.*good.efld.^(3/2)./sqrt(good.cmul).*weights)/sum(weights)

        cmul_vec=linspace(cmul_ftl_def,cmul_ftu_def,100);
        
        showoutsidepoints=1;
        colors='rgbcmk';
        %styles={'+','o'};
        figure(fignr)
        subplot(122)
        hold off
        loglog(cmul_E0(showinds),g11fit(showinds)./cmul_E0(showinds),'r-.')
        title(['r=',num2str(Geomdat.minorradiusW7AS*Geomdat.rnorm(rind)),' m,',...
               ' ex_{Er}=',num2str(ex_er,'%4.2f'),...
               ', g11_{Er}=',num2str(g11_er,'%5.2e')])
        xlabel('\nu/v')
        ylabel('|\Gamma_{11}|')
        hold on
        for efldind=1:dkdata.efldfirst{rind}.Neflds
          theefld=dkdata.efldfirst{rind}.efld(efldind);
          if theefld>=er_l && theefld<=er_u
            d11appr_sqrtnu=sqrt(cmul_vec)*g11_er/theefld^(3/2);
            d11appr_nuinv=y0./cmul_vec;
            d11appr=(d11appr_nuinv.^(-ex_er)+d11appr_sqrtnu.^(-ex_er)).^(-1/ex_er);
            loglog(cmul_vec,d11appr,colors(mod(efldind-1,length(colors))+1))
          end
          for cmulind=1:dkdata.efldfirst{rind}.Ncmuls(efldind)
            thecmul=dkdata.efldfirst{rind}.cmul{efldind}(cmulind);
            
            if theefld>=er_l && theefld<=er_u && thecmul>=cmul_ftl && thecmul<=cmul_ftu
              outside=0;
              thecolor=colorsandstyle{mod(efldind-1,length(colorsandstyle))+1};
              %thecolor=[colors(mod(efldind,length(colors))+1),styles{2}];
              themaxcolor=[colors(mod(efldind-1,length(colors))+1),'v'];
              themincolor=[colors(mod(efldind-1,length(colors))+1),'^'];
            else
              outside=1;
              thecolor=colorsandstyle{mod(efldind-1,length(colorsandstyle))+1};
              %thecolor=[colors(mod(efldind,length(colors))+1),styles{1}];
            end
            thed11=abs(dkdata.efldfirst{rind}.d11{efldind}(cmulind));
            d11max=thed11+dkdata.efldfirst{rind}.d11e{efldind}(cmulind);
            d11min=thed11-dkdata.efldfirst{rind}.d11e{efldind}(cmulind);
            if not(outside && not(showoutsidepoints))
              loglog(thecmul,thed11,thecolor)
              hold on
            end
            if not(outside)
              loglog(thecmul,d11max,themaxcolor,...
                     thecmul,d11min,themincolor)
              hold on
            end
          end
        end
        hold off
        
        d11fits.g11_er(rind)=g11_er;
        d11fits.ex_er(rind)=ex_er;
        d11fits.er_u(rind)=er_u;
        
        askthequestion=0;
        if askthequestion
          acc=input(['Do you accept the fit (return=yes, n=no, ',...
                     'x=make no fit) ? '],'s');
        else
          acc='';
        end
        acceptedfit=isempty(acc);

        if not(acceptedfit)
          do_input=1;
        end
        if cmpstr(acc,'x')
          d11fits.ex_er(rind)=NaN;
          d11fits.er_u(rind)=NaN;  
          title(['r=',num2str(Geomdat.minorradiusW7AS*Geomdat.rnorm(rind)),' m,',...
                 ' ex_{Er}=',num2str(ex_er,'%4.2f'),', g11_{Er}=',...
                 num2str(g11_er,'%5.2e'),...
                 ' No Er fit made.'])

          acceptedfit=1;
        end
      end
      
      figurefile=['fits_for_radius',num2str(rind),'of',num2str(dkdata.Nradii),'.eps'];
      disp(['Printing figure file: ',figurefile])
      pwd
      if strcmp(dkesoutputpath,'')
        print('-depsc',figurefile)
      else
        print('-depsc',[dkesoutputpath,'ps/',figurefile])
      end
    end
  end
end

