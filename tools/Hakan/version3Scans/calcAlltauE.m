function tauE=calcAlltauE(runs)

if any(diff(runs.Nspecies)~=0)
  error('Strange!')
end

for spec=1:runs.Nspecies(1)
  for ind=1:runs.NumElements
    
    e=1.6022e-19;
    mp=1.6726e-27;
    vbar=sqrt(1e3*e*2/mp);
    mbar=1.6726e-27;
    nbar=1e20;
    Phibar=1e3*e;
    
    mndPhidpsi(ind)=runs.mHats(ind,spec)*mbar*...
        runs.nHats(ind,spec)*nbar*...
        runs.dPhiHatdpsiN(ind)*Phibar/runs.psiAHat(ind);
    
    vpar=squeeze(runs.flow{ind})/runs.nHats(ind,spec)*vbar;
    
    
    BHat=runs.BHat{ind};
    dBdtheta=runs.dBHatdtheta{ind};
    dBdzeta=runs.dBHatdzeta{ind};
    h=1./BHat.^2;
    G=runs.GHat(ind);
    I=runs.IHat(ind);
    iota=runs.iota(ind);
    FSAB2=runs.FSABHat2(ind);
    
    u=calcu(BHat,G,I,runs.iota(ind),runs.NPeriods(ind));
            
    FSAu        =sum(sum(u.*h))/sum(sum(h));
    FSAvparB    =sum(sum(vpar.*BHat.*h))/sum(sum(h));
    FSAvparoverB=sum(sum(vpar./BHat.*h))/sum(sum(h)); 
    
    varpi=(FSAvparB/FSAB2-FSAvparoverB)/FSAu;
    kappa=FSAvparB/FSAB2 - varpi*runs.GHat(ind)/FSAB2;
    
    integrand=...
    (kappa-varpi*(u-G./FSAB2)).*...
        (u-G./FSAB2).*  ...
        2*(G+iota.*I)./BHat.^3 .* ...
        (G.*dBdtheta-I.*dBdzeta);
    
    tauE(ind,spec) = -mndPhidpsi(ind)*sum(sum(integrand.*h))/sum(sum(h));
  end
end


