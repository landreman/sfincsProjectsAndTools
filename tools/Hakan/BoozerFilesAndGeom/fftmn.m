function Fmn=fftmn(f)

%
% This function gives the cos and sin coefficients of f(u,v)
% written on the forms
%
% even M, even N:
% f = \sum_{n=0}^{N/2}   f^c_{0n} \cos (nv) 
%   + \sum_{n=1}^{N/2-1} f^s_{0n} \sin (nv)
%   + \sum_{m=1,n=-N/2+1}^{M/2-1,N/2} f^c_{mn} cos(mu - nv) +
%   + \sum_{m=1,n=-N/2+1}^{M/2-1,N/2} f^s_{mn} sin(mu - nv) +
%   + \sum_{n=0}^{N/2}   f^c_{M/2,n} cos(M/2*u - nv)
%   + \sum_{n=1}^{N/2-1} f^c_{M/2,n} sin(M/2*u - nv)
%
% even M, odd N:
% f = \sum_{n=0}^{(N-1)/2}   f^c_{0n} \cos (nv) 
%   + \sum_{n=1}^{(N-1)/2}   f^s_{0n} \sin (nv)
%   + \sum_{m=1,n=-(N-1)/2}^{M/2-1,(N-1)/2} f^c_{mn} cos(mu - nv) +
%   + \sum_{m=1,n=-(N-1)/2}^{M/2-1,(N-1)/2} f^s_{mn} sin(mu - nv) +
%   + \sum_{n=0}^{(N-1)/2}   f^c_{M/2,n} cos(M/2*u - nv)
%   + \sum_{n=1}^{(N-1)/2}   f^c_{M/2,n} sin(M/2*u - nv)
%
% odd M, even N:
% f = \sum_{n=0}^{N/2}   f^c_{0n} \cos (nv) 
%   + \sum_{n=1}^{N/2-1} f^s_{0n} \sin (nv)
%   + \sum_{m=1,n=-N/2+1}^{(M-1)/2,N/2} f^c_{mn} cos(mu - nv) +
%   + \sum_{m=1,n=-N/2+1}^{(M-1)/2,N/2} f^s_{mn} sin(mu - nv)
%
% odd M, odd N:
% f = \sum_{n=0}^{(N-1)/2}   f^c_{0n} \cos (nv) 
%   + \sum_{n=1}^{(N-1)/2}   f^s_{0n} \sin (nv)
%   + \sum_{m=1,n=-(N-1)/2}^{(M-1)/2,(N-1)/2} f^c_{mn} cos(mu - nv) +
%   + \sum_{m=1,n=-(N-1)/2}^{(M-1)/2,(N-1)/2} f^s_{mn} sin(mu - nv)
%
%
% It is assumed that the u and v are uniform vectors
% u=(0:M-1)/M*2*pi
% v=(0:N-1)/N*2*pi
% Because u=v=0 is included in the grid, there are more cos than sin coefficients.
%

[Nu,Nv]=size(f);

%if Nv/2-floor(Nv/2)~=0 || Nu/2-floor(Nu/2)~=0
%  error('Nv and Nu must be even numbers!')
%end
%if Nv~=Nu
%  error('f must be square!')
%end

F=fliplr(fftshift(fft2(f),2));

mmax=floor(Nu/2);
nmax=floor(Nv/2);
Nu_is_even=not(mod(Nu,2));
Nv_is_even=not(mod(Nv,2));
%nrange_old=-Nv/2+1:Nv/2 %only worked for even Nv
nrange=-nmax+Nv_is_even:nmax;
mrange=0:mmax;
[m,n]=ndgrid(mrange,nrange);
Fmn.c=NaN*zeros(size(m));
Fmn.s=NaN*zeros(size(m));

Fmn.c(1:end,:) =  2/(Nu*Nv) * real(F(1:mmax+1,:));
Fmn.s(1:end,:) = -2/(Nu*Nv) * imag(F(1:mmax+1,:));

m0ind=1;
n0ind=ceil(Nv/2);

Fmn.c(m0ind,1:n0ind-1)=NaN;
Fmn.s(m0ind,1:n0ind)=NaN;
Fmn.c(m0ind,n0ind)=Fmn.c(m0ind,n0ind)/2;
if Nv_is_even
  Fmn.s(m0ind,end)=NaN;
  Fmn.c(m0ind,end)=Fmn.c(m0ind,end)/2;
end

if Nu_is_even
  Fmn.c(end,1:n0ind-1)=NaN;
  Fmn.s(end,1:n0ind)=NaN;
  Fmn.c(end,n0ind)=Fmn.c(end,n0ind)/2;
  if Nv_is_even
    Fmn.s(end,end)=NaN;
    Fmn.c(end,end)=Fmn.c(end,end)/2;
  end
end

Fmn.m=m;
Fmn.n=n;
Fmn.m0ind=m0ind;
Fmn.n0ind=n0ind;
