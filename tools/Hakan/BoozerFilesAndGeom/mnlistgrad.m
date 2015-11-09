function varargout=mnlistgrad(Fmn,Nperiods)


if strcmp(class(Fmn),'mnmat')
  varargout{1}=grad(Fmn,Nperiods);
else
  
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
  ddvFmn.data     = [ Fmn.n(cind)*Nperiods.*Fmn.data(cind),...
                     -Fmn.n(sind)*Nperiods.*Fmn.data(sind)];  

  if nargout==2
    varargout{1}=dduFmn;
    varargout{2}=ddvFmn;
  elseif nargout==1
    varargout{1}=[dduFmn,ddvFmn];
  end
end