function varargout=mngrad(Fmn,Nperiods)


if isfield(Fmn,'cosparity')  %Fmn is of list form
  
  cind=find(Fmn.cosparity & not(Fmn.m==0 & Fmn.n==0));
  sind=find(not(Fmn.cosparity) & not(Fmn.m==0 & Fmn.n==0));
  
  dduFmn.m        = [Fmn.m(cind), Fmn.m(sind)];
  dduFmn.n        = [Fmn.n(cind), Fmn.n(sind)];
  dduFmn.cosparity= not([Fmn.cosparity(cind), Fmn.cosparity(sind)]);
  dduFmn.data     = [-Fmn.m(cind).*Fmn.data(cind), ...
                      Fmn.m(sind).*Fmn.data(sind)];
  
  ddvFmn.m        = [Fmn.m(cind), Fmn.m(sind)];
  ddvFmn.n        = [Fmn.n(cind), Fmn.n(sind)];
  ddvFmn.cosparity= not([Fmn.cosparity(cind), Fmn.cosparity(sind)]);
  ddvFmn.data     = [-Fmn.n(cind)*Nperiods.*Fmn.data(cind),...
                      Fmn.n(sind)*Nperiods.*Fmn.data(sind)];  
  
else %Fmn is of matrix form
  dduFmn.s = -Fmn.c.*Fmn.m;
  dduFmn.c =  Fmn.s.*Fmn.m;

  ddvFmn.s = -Fmn.c.*Fmn.n*Nperiods;
  ddvFmn.c =  Fmn.s.*Fmn.n*Nperiods;
  
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
  dduFmn.c(end,Fmn.n0ind)=0;
  dduFmn.c(Fmn.m0ind,end)=0;
  dduFmn.c(end,end)=0;
  dduFmn.s(Fmn.m0ind,Fmn.n0ind)=NaN;
  dduFmn.s(end,Fmn.n0ind)=NaN;
  dduFmn.s(Fmn.m0ind,end)=NaN;
  dduFmn.s(end,end)=NaN;
  
  ddvFmn.c(Fmn.m0ind,Fmn.n0ind)=0;
  ddvFmn.c(end,Fmn.n0ind)=0;
  ddvFmn.c(Fmn.m0ind,end)=0;
  ddvFmn.c(end,end)=0;
  ddvFmn.s(Fmn.m0ind,Fmn.n0ind)=NaN;
  ddvFmn.s(end,Fmn.n0ind)=NaN;
  ddvFmn.s(Fmn.m0ind,end)=NaN;
  ddvFmn.s(end,end)=NaN;
  
end

if nargout==2
  varargout{1}=dduFmn;
  varargout{2}=ddvFmn;
elseif nargout==1
  varargout{1}=[dduFmn,ddvFmn];
end