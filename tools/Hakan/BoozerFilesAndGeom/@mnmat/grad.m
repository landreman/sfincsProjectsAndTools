function varargout=grad(Fmn,Nperiods)

  Nu_even=isnan(Fmn.s(end,Fmn.n0ind));
  Nv_even=not(mod(size(Fmn.c,2),2));

  dduFmn.s = -Fmn.c.*Fmn.m;
  dduFmn.c =  Fmn.s.*Fmn.m;

  ddvFmn.s =  Fmn.c.*Fmn.n*Nperiods;
  ddvFmn.c = -Fmn.s.*Fmn.n*Nperiods;
  
  dduFmn.m=Fmn.m;
  ddvFmn.m=Fmn.m;
  dduFmn.n=Fmn.n;
  ddvFmn.n=Fmn.n;
  dduFmn.m0ind=Fmn.m0ind;
  ddvFmn.m0ind=Fmn.m0ind;
  dduFmn.n0ind=Fmn.n0ind;
  ddvFmn.n0ind=Fmn.n0ind;
  
  %Let us force the size of the data matrix to be the same
  %In principle one should increase it because there is more sine 
  %information than it can accommodate now, but that would be unpractical
  %if we want to transform back to real space.
  
  dduFmn.c(Fmn.m0ind,Fmn.n0ind)=0;
  dduFmn.s(Fmn.m0ind,Fmn.n0ind)=NaN;
  if Nu_even
    dduFmn.c(end,Fmn.n0ind)=0;
    dduFmn.s(end,Fmn.n0ind)=NaN;
    if Nv_even
      dduFmn.c(end,end)=0;
      dduFmn.s(end,end)=NaN;
    end
  end
  if Nv_even
    dduFmn.c(Fmn.m0ind,end)=0;
    dduFmn.s(Fmn.m0ind,end)=NaN;
  end
  
  ddvFmn.c(Fmn.m0ind,Fmn.n0ind)=0;
  ddvFmn.s(Fmn.m0ind,Fmn.n0ind)=NaN;
  if Nu_even
    ddvFmn.c(end,Fmn.n0ind)=0;
    ddvFmn.s(end,Fmn.n0ind)=NaN;
    if Nv_even
      ddvFmn.c(end,end)=0;
      ddvFmn.s(end,end)=NaN;
    end
  end
  if Nv_even
    ddvFmn.c(Fmn.m0ind,end)=0;
    ddvFmn.s(Fmn.m0ind,end)=NaN;
  end
  dduFmn=mnmat(dduFmn);
  ddvFmn=mnmat(ddvFmn);

if nargout==2
  varargout{1}=dduFmn;
  varargout{2}=ddvFmn;
elseif nargout==1
  varargout{1}=[dduFmn,ddvFmn];
end