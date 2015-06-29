function varargout=ifftmn(Fmns,varargin)
% if length(varargout)=length(Fmn)+2 then the two last output arguments are the
% poloidal and toroidal dicretisation matrices u and v.

% The input Nperiods=varargin{1} is required if the output
% toroidal discretisation matrix v is wanted

% The input 3 and 4 can be used to specify a different size of the output
% compared to the input. Beware of information loss if you use this option.

Ninput=length(Fmns);

if isfield(Fmns,'cosparity')
  %list input was given!
  Fmns=mnmatrix(Fmns);
end

if nargin>=2
  Nperiods=varargin{1};
end
if nargout==Ninput+2 && nargin < 4
  error(['Nperiods, Nu and Nv must be given as input if the ',...
         'toroidal discretisation matrix is requested as output!'])
end

for inputind=1:Ninput
  Fmn=Fmns(inputind);
    
  Nv=size(Fmn.c,2);
  Nu=(size(Fmn.c,1)-1)*2;
  if Nv~=size(Fmn.s,2) || Nu~=(size(Fmn.s,1)-1)*2;
    error('Sizes of Fmn.c and Fmn.s are inconsistent!')
  end

  if nargin==4
    Nu_in=Nu;
    Nv_in=Nv;
    Nu=varargin{2};
    Nv=varargin{3};
    if (Nu/2-floor(Nu/2)~=0) || (Nv/2-floor(Nv/2)~=0)
    error('The size must be an even number!')
    end

    nrange=-Nv/2+1:Nv/2;
    mrange=0:Nu/2;
    [m,n]=ndgrid(mrange,nrange);
    m0ind=1;
    n0ind=Nv/2;
    
    Fmn_in=Fmn;
    if Nu>Nu_in && Nv>Nv_in
      du=(Nu-Nu_in)/2;
      dv=(Nv-Nv_in)/2;
      Fmn.c=zeros(Nu/2+1,Nv);
      Fmn.s=zeros(Nu/2+1,Nv);
      Fmn_in.c(end,1:Nv_in/2-1)=0; %remove these NaNs before enlarging
      Fmn_in.s(end,1:Nv_in/2)=0; %remove these NaNs before enlarging
      Fmn_in.s(1,end)=0; %remove this NaN before enlarging
      Fmn_in.s(end,end)=0; %remove this NaN before enlarging
      
      Fmn.c(1:end-du,1+dv:end-dv)=Fmn_in.c;
      Fmn.s(1:end-du,1+dv:end-dv)=Fmn_in.s;
      
      Fmn.c(1,1:n0ind-1)=NaN;
      Fmn.c(end,1:n0ind-1)=NaN;
      Fmn.s(1,1:n0ind)=NaN;
      Fmn.s(end,1:n0ind)=NaN;
      Fmn.s(1,end)=NaN;
      Fmn.s(end,end)=NaN;
      
    elseif Nu<Nu_in && Nv<Nv_in
      disp('Warning: Requested spatial discretisation may not resolve all modes!')
      du=(Nu_in-Nu)/2;
      dv=(Nv_in-Nv)/2;
      Fmn.c=Fmn_in.c(1:end-du,1+dv:end-dv);
      Fmn.s=Fmn_in.s(1:end-du,1+dv:end-dv);
      Fmn.c(end,1:n0ind-1)=NaN;
      Fmn.s(end,1:n0ind)=NaN;
      Fmn.s(1,end)=NaN;
      Fmn.s(end,end)=NaN;
    elseif Nu==Nu_in && Nv==Nv_in
      Fmn=Fmn_in;
    else
      error('The case where Nu is increased and Nv decreased is not implemented yet!')
    end
  else
    nrange=-Nv/2+1:Nv/2;
    mrange=0:Nu/2;
    [m,n]=ndgrid(mrange,nrange);
    m0ind=1;
    n0ind=Nv/2;
  end
  
  uvec=(0:Nu-1)/Nu*2*pi;
  vvec=(0:Nv-1)/Nv*2*pi;
  [u,v]=ndgrid(uvec,vvec);

  %Normal assembly, is very slow
  if 0
    f=zeros(Nu);
    for m=1:Nu/2-1
      %disp([num2str(m),' / ',num2str(Nu/2-1)])
      for n= -Nv/2+1:Nv/2
        f=f+Fmn.c(m+1,n+Nv/2)*cos(m*u-n*v)+Fmn.s(m+1,n+Nv/2)*sin(m*u-n*v);
      end
    end
    
    for n=0:Nv/2
      f=f+Fmn.c(1,n+Nv/2)*cos(-n*v)+Fmn.c(end,n+Nv/2)*cos(Nu/2*u-n*v);
    end
    
    for n=1:Nv/2-1
      f=f+Fmn.s(1,n+Nv/2)*sin(-n*v)+Fmn.s(end,n+Nv/2)*sin(Nu/2*u-n*v);
    end
  else %ifft assembly, is much faster
    
    Fmn.c(m0ind,n0ind)=Fmn.c(m0ind,n0ind)*2;
    Fmn.c(end,n0ind)=Fmn.c(end,n0ind)*2;
    Fmn.c(m0ind,end)=Fmn.c(m0ind,end)*2;
    Fmn.c(end,end)=Fmn.c(end,end)*2;
    
    Fmn.c(m0ind,1:n0ind-1)=fliplr(Fmn.c(m0ind,n0ind+1:end-1));
    Fmn.s(m0ind,n0ind)=0;
    Fmn.s(m0ind,end)=0;
    Fmn.s(m0ind,1:n0ind-1)=-fliplr(Fmn.s(m0ind,n0ind+1:end-1));
    
    Fmn.c(end,1:n0ind-1)=fliplr(Fmn.c(end,n0ind+1:end-1));
    Fmn.s(end,n0ind)=0;
    Fmn.s(end,end)=0;
    Fmn.s(end,1:n0ind-1)=-fliplr(Fmn.s(end,n0ind+1:end-1));

    Fmncos=zeros(Nu,Nv);
    Fmnsin=zeros(Nu,Nv);
    Fmncos(1:Nu/2+1,:)=Fmn.c;
    Fmnsin(1:Nu/2+1,:)=Fmn.s;
    
    Fmncos(Nu/2+2:end,n0ind)=flipud(Fmncos(2:Nu/2,n0ind));
    Fmncos(Nu/2+2:end,end)  =flipud(Fmncos(2:Nu/2,end));
    Fmnsin(Nu/2+2:end,n0ind)=-flipud(Fmnsin(2:Nu/2,n0ind));
    Fmnsin(Nu/2+2:end,end)  =-flipud(Fmnsin(2:Nu/2,end));

    
    Fmncos(Nu/2+2:end,1:Nv/2-1)     =fliplr(flipud(Fmncos(2:Nu/2,Nv/2+1:end-1)));
    Fmncos(Nu/2+2:end,Nv/2+1:end-1) =fliplr(flipud(Fmncos(2:Nu/2,1:Nv/2-1)));
    Fmnsin(Nu/2+2:end,1:Nv/2-1)     =-fliplr(flipud(Fmnsin(2:Nu/2,Nv/2+1:end-1)));
    Fmnsin(Nu/2+2:end,Nv/2+1:end-1) =-fliplr(flipud(Fmnsin(2:Nu/2,1:Nv/2-1)));
    
    F=Nu*Nv/2*(Fmncos-i*Fmnsin);
    
    f=real(ifft2(fftshift(fliplr(F),2)));
  end

  varargout{inputind}=f;
end

if nargout==Ninput+2
  varargout{Ninput+1}=u;
  varargout{Ninput+2}=v/Nperiods;
end