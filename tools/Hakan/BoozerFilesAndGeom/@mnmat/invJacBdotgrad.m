function out=invJacBdotgrad(arg,iota,Nperiods,varargin)

% Assumes Boozer coordinates
% solves arg = iota d(out)/dtheta +d(out)/dzeta

% The optional 4th argument zeroout_Deltaiota can be used 
% to zero out mn components 
% which have abs(n/m*NPeriods - iota) < zeroout_Deltaiota.
%

if nargin>3
  zeroout_Deltaiota=varargin{1};
else
  zeroout_Deltaiota=-1; %means that it will not be used
end

sz=size(arg);
if not(mod(sz(1),2)) || not(mod(sz(2),2))
  error('information loss if Nu or Nv even!')
end

if arg.c(out.m0ind,out.n0ind)~=0
  error('The 00 component of the input has to be = 0!')
end

out.m=arg.m;
out.n=arg.n;
out.m0ind=arg.m0ind;
out.n0ind=arg.n0ind;
out.c = -arg.s./(iota*arg.m-Nperiods*arg.n);
out.s = arg.c./(iota*arg.m-Nperiods*arg.n);

m_tmp=out.m;
m_tmp(out.m0ind,:)=10*eps; %to avoid division by zero
badinds=find(abs(out.n./m_tmp*NPeriods - iota)<zeroout_Deltaiota);
out.c(badinds)=0;
out.s(badinds)=0; 

out.c(out.m0ind,out.n0ind)=0;
out.s(out.m0ind,out.n0ind)=NaN;

out=mnmat(out);

