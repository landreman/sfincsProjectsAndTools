function [us,umnmats]=calcu(hs,Gs,Is,iotas,NPeriod,varargin)
% This solves the equation
%
%   B\cdot\nabla u = - B\times\nabla\chi\cdot\nabla h
%
% given that <B^2 u> = 0
%
% Very often we are interested to solve the equation with 
% h=1/B^2, in which case one puts fftmn(1./B.^2) as the 
% first argument.
%
% The optional 6th argument zeroout_Deltaiota can be used 
% to zero out mn components of wmn
% which have abs(n/m*NPeriods - iota) < zeroout_Deltaiota.
%
% The optional 7th argument can be used to change from the
% poloidal flux \chi as the flux surface label in the above
% equation to the toroidal flux \psi by giving the strings
% 'poloidal flux' (default) or 'toroidal flux'.


if nargin>5
  zeroout_Deltaiota=varargin{1};
else
  zeroout_Deltaiota=-1; %means that it will not be used
end
if nargin>6
  radialvariable=varargin{2};
else
  radialvariable='poloidal flux';
end
if strcmp(radialvariable,'poloidal flux')
  iotaexp=1;
elseif strcmp(radialvariable,'toroidal flux')
  iotaexp=0;
else
  error('Radial variable not recognised!')
end

if length(NPeriod)==1
  NPeriods=NPeriod*ones(size(iotas));
else
  NPeriods=NPeriod;
end

if length(iotas)==1
    G=Gs;
    I=Is;
    iota=iotas;
    N=NPeriods;
    
    hmn=hs; % often fftmn(1./B.^2);
    
    umn.c=iota^iotaexp*(G*hmn.m + I*hmn.n * N)./(hmn.n * N - iota*hmn.m).*hmn.c;
    umn.s=iota^iotaexp*(G*hmn.m + I*hmn.n * N)./(hmn.n * N - iota*hmn.m).*hmn.s;
    
    umn.c(hmn.m0ind,hmn.n0ind)=0;
    umn.s(hmn.m0ind,hmn.n0ind)=NaN;
    
    umn.m=hmn.m;
    umn.n=hmn.n;
    umn.m0ind=hmn.m0ind;
    umn.n0ind=hmn.n0ind;

    m_tmp=hmn.m;
    m_tmp(hmn.m0ind,:)=10*eps; %to avoid division by zero
    badinds=find(abs(hmn.n./m_tmp*NPeriods - iota)<zeroout_Deltaiota);
    umn.c(badinds)=0;
    umn.s(badinds)=0; 
    
    umnmats=mnmat(umn);
    us=ifftmn(umnmats);
else %inputs were given in vectors
  for ind=1:length(iotas)
    G=Gs(ind);
    I=Is(ind);
    iota=iotas(ind);
    N=NPeriods(ind);
    
    hmn=hs(ind); % often fftmn(1./B.^2);
    
    umn.c=iota^iotaexp*(G*hmn.m + I*hmn.n * N)./(hmn.n * N - iota*hmn.m).*hmn.c;
    umn.s=iota^iotaexp*(G*hmn.m + I*hmn.n * N)./(hmn.n * N - iota*hmn.m).*hmn.s;
    
    umn.c(hmn.m0ind,hmn.n0ind)=0;
    umn.s(hmn.m0ind,hmn.n0ind)=NaN;

    umn.m=hmn.m;
    umn.n=hmn.n;
    umn.m0ind=hmn.m0ind;
    umn.n0ind=hmn.n0ind;

    m_tmp=hmn.m;
    m_tmp(hmn.m0ind,:)=10*eps; %to avoid division by zero
    badinds=find(abs(hmn.n./m_tmp*NPeriods - iota)<zeroout_Deltaiota);
    umn.c(badinds)=0;
    umn.s(badinds)=0; 
    
    umnmats{ind}=mnmat(umn);
    us{ind}=ifftmn(umnmats{ind});
  end
end