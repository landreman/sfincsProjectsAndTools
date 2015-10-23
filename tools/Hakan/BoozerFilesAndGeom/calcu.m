function us=calcu(Bs,Gs,Is,iotas,NPeriod)

if length(iotas)==1
    B=Bs;
    G=Gs;
    I=Is;
    iota=iotas;
     
    hmn=fftmn(1./B.^2);
    
    umn.c=iota*(G*hmn.m + I*hmn.n * NPeriod)./(hmn.n * NPeriod - iota*hmn.m).*hmn.c;
    umn.s=iota*(G*hmn.m + I*hmn.n * NPeriod)./(hmn.n * NPeriod - iota*hmn.m).*hmn.s;
    
    umn.c(hmn.m0ind,hmn.n0ind)=0;
    umn.s(hmn.m0ind,hmn.n0ind)=NaN;
    
    umn.m=hmn.m;
    umn.n=hmn.n;
    umn.m0ind=hmn.m0ind;
    umn.n0ind=hmn.n0ind;

    us=ifftmn(umn);
else %inputs were given in vectors
  for ind=1:length(iotas)
    B=Bs{ind};
    G=Gs(ind);
    I=Is(ind);
    iota=iotas(ind);
    
    hmn=fftmn(1./B.^2);
    
    umn.c=iota*(G*hmn.m + I*hmn.n * NPeriod)./(hmn.n * NPeriod - iota*hmn.m).*hmn.c;
    umn.s=iota*(G*hmn.m + I*hmn.n * NPeriod)./(hmn.n * NPeriod - iota*hmn.m).*hmn.s;
    
    umn.c(hmn.m0ind,hmn.n0ind)=0;
    umn.s(hmn.m0ind,hmn.n0ind)=NaN;

    umn.m=hmn.m;
    umn.n=hmn.n;
    umn.m0ind=hmn.m0ind;
    umn.n0ind=hmn.n0ind;

    us{ind}=ifftmn(umn);
  end
end