%function g4out=g4(Geom,rind,lambda,Ntheta,Nzeta)
function g4out=g4(B,lambda,NPeriods,I,G,iota)

% ---------------------------------------------------------------------------------------
% Calculate g4 from harmonics of 1/sqrt(1-lambda B).
% Define w as
% \nabla_\parallel w = \nabla B \times \vector{B} \cdot \nabla \psi 
%
% Then g4 =< sqrt(1-lambda B) w >
%
% ---------------------------------------------------------------------------------------

%First, load data from Geom:
%NPeriods=Geom.Nperiods;
%I=Geom.Btheta(rind);
%G=Geom.Bphi(rind);
%iota=Geom.iota(rind);
%Ntheta=151;
%Nzeta=151;
%B=ifftmn(mnmat(Geom,rind,'B',Ntheta,Nzeta,'forceSize'));


%now do the calculation
w = zeros(size(B));
h=1./sqrt(1-lambda*B);

hmn=unmnmat(fftmn(h));

wmn.c=(G*hmn.m + I*hmn.n * NPeriods)./(hmn.n * NPeriods - iota*hmn.m) .* hmn.c;
wmn.c(hmn.m0ind,hmn.n0ind)=0;
wmn.s=(G*hmn.m + I*hmn.n * NPeriods)./(hmn.n * NPeriods - iota*hmn.m) .* hmn.s;
wmn.s(hmn.m0ind,hmn.n0ind)=NaN; %Important to set. Indicates M odd/even
wmn.m=hmn.m;wmn.n=hmn.n;
wmn.m0ind=hmn.m0ind;wmn.n0ind=hmn.n0ind;
mnmats.w=mnmat(wmn);
w=ifftmn(mnmats.w);

if size(w,1)~=151
  h
  %fftmn(h)
  size(h)
  size(hmn.c)
  size(hmn.s)
  size(wmn.c)
  size(wmn.s)
  hmn.n0ind
  wmn.s(hmn.m0ind,hmn.n0ind)
end

%FSA:
g4out = sqrt(1-lambda*B).*w;
%g4out = sum(sum(B.^(-2).*sqrt(1-lambda*B).*w)) / sum(sum(B.^(-2)));

  