function [tauP,tauE,kappa,varpi]=calcAlltau(runs,varargin)

if nargin==1
  Ntheta=runs.Ntheta(1);
  Nzeta=runs.Nzeta(1);
else
  Ntheta=varargin{1};
  Nzeta=varargin{2};%We might need more zeta resolution when we multiply the quantities in real space
end
  
if any(diff(runs.Nspecies)~=0)
  error('Strange!')
end

mndPhidpsi=NaN*zeros(runs.NumElements,1);
tauP=NaN*zeros(runs.NumElements,runs.Nspecies(1));
tauE=NaN*zeros(runs.NumElements,runs.Nspecies(1));
tauE2=NaN*zeros(runs.NumElements,runs.Nspecies(1));
varpi=NaN*zeros(runs.NumElements,runs.Nspecies(1));
kappa=NaN*zeros(runs.NumElements,runs.Nspecies(1));
NPeriods=runs.NPeriods(1);

for spec=1:runs.Nspecies(1)
  for ind=1:runs.NumElements
    if runs.finished(ind)
      e=1.6022e-19;
      mp=1.6726e-27;
      vbar=sqrt(1e3*e*2/mp);
      mbar=1.6726e-27;
      nbar=1e20;
      Tbar=1e3*e;
      Phibar=1e3;
      
      mndPhidpsi(ind)=runs.mHats(ind,spec)*mbar*...
          runs.nHats(ind,spec)*nbar*...
          runs.dPhiHatdpsiN(ind)*Phibar/runs.psiAHat(ind);
      
      vpar=ifftmn(fftmn(squeeze(runs.flow{ind})/runs.nHats(ind,spec)*vbar),...
                  NaN,Ntheta,Nzeta);
      vparoverB_mn=fftmn(squeeze(runs.flow{ind})/runs.nHats(ind,spec)*vbar...
                             ./runs.BHat{ind});
      vparoverB3_mn=fftmn(squeeze(runs.flow{ind})/runs.nHats(ind,spec)*vbar...
                             ./runs.BHat{ind}.^3);
      vparoverB=ifftmn(vparoverB_mn,...
                  NaN,Ntheta,Nzeta);
      [dvparoverBdtheta,dvparoverBdzeta]=ifftmn(mnmat(mnlistgrad(vparoverB_mn,NPeriods)),...
                  NaN,Ntheta,Nzeta);
      [dvparoverB3dtheta,dvparoverB3dzeta]=ifftmn(mnmat(mnlistgrad(vparoverB3_mn,NPeriods)),...
                  NaN,Ntheta,Nzeta);
      n1overn0=ifftmn(fftmn(squeeze(runs.densityPerturbation{ind})/runs.nHats(ind,spec)),NaN,Ntheta,Nzeta);
      
      
      ptilde=ifftmn(fftmn(runs.pressureAnisotropy{ind}*nbar*Tbar),NaN,Ntheta,Nzeta);
      B=ifftmn(fftmn(runs.BHat{ind}),NaN,Ntheta,Nzeta);
      dBdtheta_mn=fftmn(runs.dBHatdtheta{ind});
      dBdzeta_mn=fftmn(runs.dBHatdzeta{ind});
      dBdtheta=ifftmn(fftmn(runs.dBHatdtheta{ind}),NaN,Ntheta,Nzeta);
      dBdzeta=ifftmn(fftmn(runs.dBHatdzeta{ind}),NaN,Ntheta,Nzeta);

      [dBdtheta_amn,dBdzeta_amn]=grad(fftmn(runs.BHat{ind}),NPeriods);
      [dBdtheta_a,dBdzeta_a]=ifftmn(grad(fftmn(runs.BHat{ind}),NPeriods),...
                  NaN,Ntheta,Nzeta);
      
      %dBdzeta_mn.c
      %dBdzeta_amn.c
      %dBdzeta_mn.s
      %dBdzeta_amn.s
      %pause
      
      h=1./B.^2;
      G=runs.GHat(ind);
      I=runs.IHat(ind);
      iota=runs.iota(ind);
      FSAB2=runs.FSABHat2(ind);
      gamma=-G/FSAB2;
      
      u=calcu(fftmn(h),G,I,runs.iota(ind),runs.NPeriods(ind));
      [dudtheta,dudzeta]=ifftmn(grad(fftmn(u),NPeriods));
      
      
      %%%%%%%%%%%%%%% Calculate tauP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      NTVKernel=FSAB2/(4*pi^2)/(G+iota*I)*...
                1./B.^2.*(...
                  (gamma+u).*B.*(iota*dBdtheta+dBdzeta)+iota./B.*(G.*dBdtheta-I.*dBdzeta));
      
      tauP(ind,spec)=4*pi^2/Ntheta/Nzeta*sum(sum(NTVKernel.*ptilde));
      
      %%%%%%%%%%%%%%% Calculate varpi and kappa %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FSAu        =sum(sum(u.*h))/sum(sum(h));
      FSAvparB    =sum(sum(vpar.*B.*h))/sum(sum(h));
      FSAvparoverB=sum(sum(vpar./B.*h))/sum(sum(h));       
      
      %FSAvparB_Sfx=runs.FSABFlow(ind)/runs.nHats(ind,spec)*vbar
      
      %FSAvparB*runs.nHats(ind,spec)/vbar %=FSABFlow OK!
      varpi(ind,spec)=(FSAvparB/FSAB2-FSAvparoverB)/FSAu;
      kappa(ind,spec)=FSAvparB/FSAB2 - varpi(ind,spec)*G/FSAB2;
      
      grn=grad(fftmn(n1overn0),NPeriods);
      
      
      fac1=varpi(ind,spec)*iota*h;
      fac2=varpi(ind,spec)*(u-G./FSAB2)-kappa(ind,spec);
      
      rhs=grn(1).*( G*fac1+fac2*iota)...
          +grn(2).*(-I*fac1+fac2);
      
      %invJacBdotgrad(rhs,iota,NPeriods)

      Vpar1overB_try=ifftmn(invJacBdotgrad(rhs,iota,NPeriods));
      
      %%%%%%%%%%%%%%% Calculate tauE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      integrand=...
          (kappa(ind,spec)-varpi(ind,spec)*(u-G./FSAB2)).*...
          (u-G./FSAB2).*  ...
          2/(G+iota.*I).*B .* ...
          (G.*dBdtheta-I.*dBdzeta);
      integrand2=-(u-G./FSAB2).*B.^4/(G+iota.*I) .* ...
          (G.*dvparoverB3dtheta-I.*dvparoverB3dzeta);
      integrand3=-(u-G./FSAB2).*B.^2/(G+iota.*I) .* ...
          (G.*dvparoverBdtheta-I.*dvparoverBdzeta);
      %integrand4=-(-G./FSAB2).*B.^2/(G+iota.*I) .* ...
      %    (G.*dvparoverBdtheta-I.*dvparoverBdzeta);
      integrand4=-vparoverB...
          .*B.^2/(G+iota.*I) .* ...
          (G.*dudtheta-I.*dudzeta);
      integrand5=-(kappa(ind,spec)-varpi(ind,spec)*(u-G./FSAB2))...
          .*B.^2/(G+iota.*I) .* ...
          (G.*dudtheta-I.*dudzeta);
      
      
      tauE1(ind,spec) = -mndPhidpsi(ind)*sum(sum(integrand.*h))/sum(sum(h));
      tauE2(ind,spec) = -mndPhidpsi(ind)*sum(sum(integrand2.*h))/sum(sum(h));
      tauE3(ind,spec) = -mndPhidpsi(ind)*sum(sum(integrand3.*h))/sum(sum(h));
      tauE4(ind,spec) = -mndPhidpsi(ind)*sum(sum(integrand4.*h))/sum(sum(h));
      tauE5(ind,spec) = -mndPhidpsi(ind)*sum(sum(integrand5.*h))/sum(sum(h));
      
      tauE(ind,spec) = tauE2(ind,spec);
      
      doit=0;
      if doit
        
        %normth=max(max(max(abs(dBdtheta))), max(max(abs(dBdtheta_a))))
        %T1=dBdtheta+normth*2;
        %T2=dBdtheta_a+normth*2;
        %normz=max(max(max(abs(dBdzeta))), max(max(abs(dBdzeta_a))))
        %T1=dBdzeta+normz*2;
        %T2=dBdzeta_a+normz*2;
        %T1=(kappa(ind,spec)-varpi(ind,spec)*(runs.uHat{ind}-G./FSAB2));
        T1=(kappa(ind,spec)-varpi(ind,spec)*(u-G./FSAB2));
        T2=vparoverB;
        T3=squeeze(runs.flow{ind})/runs.nHats(ind,spec)*vbar./runs.BHat{ind};
        maxx=max(max(max(abs(T1))), max(max(abs(T2))));
        
        %fig(1)
        %surf(u)
        %shading flat;view(0,90);colorbar;title('u')
        fig(5)
        surf(T1)
        shading flat;view(0,90);colorbar;title(num2str(runs.rN(ind)))
        fig(2)
        surf((T2+2*maxx)./(T1+2*maxx))
        shading flat;view(0,90);colorbar;title(num2str(runs.rN(ind)))
        fig(3)
        surf(T2-T1)
        shading flat;view(0,90);colorbar;title(num2str(runs.rN(ind)))
        fig(4)
        surf(Vpar1overB_try)
        %surf(runs.densityPerturbation{ind})
        shading flat;view(0,90);colorbar;title(num2str(runs.rN(ind)))
        
        Umn=fftmn(u);
        USmn=fftmn(runs.uHat{ind});
        DUmn=fftmn(u-runs.uHat{ind});
        
        %DUmn.c
        %DUmn.s
        %Umn.c
        %Umn.s
        
        drawnow
        %pause
      end
    end
  end
end

%tauE'
%tauE2'
%tauE3'
%tauE4'
%tauE5'
%error('stopppppp')


