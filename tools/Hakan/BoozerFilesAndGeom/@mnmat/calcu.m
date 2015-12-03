function [us,umnmats]=calcu(hs,Gs,Is,iotas,NPeriod)

if length(NPeriod)==1
  NPeriods=NPeriod*ones(size(iotas));
else
  NPeriods=NPeriod;
end

if length(iotas)==1
    %B=Bs;
    G=Gs;
    I=Is;
    iota=iotas;
    N=NPeriods;
    
    hmn=hs;
    %hmn=fftmn(1./B.^2);
    
    umn.c=iota*(G*hmn.m + I*hmn.n * N)./(hmn.n * N - iota*hmn.m).*hmn.c;
    umn.s=iota*(G*hmn.m + I*hmn.n * N)./(hmn.n * N - iota*hmn.m).*hmn.s;
    
    umn.c(hmn.m0ind,hmn.n0ind)=0;
    umn.s(hmn.m0ind,hmn.n0ind)=NaN;
    
    umn.m=hmn.m;
    umn.n=hmn.n;
    umn.m0ind=hmn.m0ind;
    umn.n0ind=hmn.n0ind;

    umnmats=mnmat(umn);
    us=ifftmn(umnmats);
else %inputs were given in vectors
  for ind=1:length(iotas)
    %B=Bs{ind};
    G=Gs(ind);
    I=Is(ind);
    iota=iotas(ind);
    N=NPeriods(ind);
    
    hmn=hs(ind);
    %hmn=fftmn(1./B.^2);
    
    umn.c=iota*(G*hmn.m + I*hmn.n * N)./(hmn.n * N - iota*hmn.m).*hmn.c;
    umn.s=iota*(G*hmn.m + I*hmn.n * N)./(hmn.n * N - iota*hmn.m).*hmn.s;
    
    umn.c(hmn.m0ind,hmn.n0ind)=0;
    umn.s(hmn.m0ind,hmn.n0ind)=NaN;

    umn.m=hmn.m;
    umn.n=hmn.n;
    umn.m0ind=hmn.m0ind;
    umn.n0ind=hmn.n0ind;

    umnmats{ind}=mnmat(umn);
    us{ind}=ifftmn(umnmats{ind});
  end
end