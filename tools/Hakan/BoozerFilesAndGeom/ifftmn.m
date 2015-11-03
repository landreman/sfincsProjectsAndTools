function varargout=ifftmn(Fmns,varargin)
% This routine can process more than one Fmn input by giving a vector as input:
% Fmns=[Fmn1,Fmn2,Fmn3...]
% if length(varargout)=length(Fmns)+2 then the two last output arguments are the
% poloidal and toroidal dicretisation matrices u and v.

% The input Nperiods=varargin{1} is required if the output
% toroidal discretisation matrix v is wanted

% The input 3 and 4 can be used to specify a different size of the output
% compared to the input. Beware of information loss if you use this option.

Ninput=length(Fmns);

if nargin>=2
  Nperiods=varargin{1};
end
if nargout==Ninput+2 && nargin < 4
  error(['Nperiods, Nu and Nv must be given as input if the ',...
         'toroidal discretisation matrix is requested as output!'])
end

listinput=0;
if isfield(Fmns,'cosparity')
  %list input was given!
  listinput=1;
  if nargin==4
    Fmns=mnmatrix(Fmns,varargin{2},varargin{3});
  else
    Fmns=mnmatrix(Fmns);
  end
end


for inputind=1:Ninput
  Fmn=Fmns(inputind);
  make_normal_assembly=0; %0 is default. Only special cases will switch to 1 later
    
  if any(size(Fmn.c)~=size(Fmn.s))
    error('Sizes of Fmn.c and Fmn.s are inconsistent!')
  end
  Nv=size(Fmn.c,2);
  n0ind_tmp=ceil(Nv/2);
  % The boolean isnan(Fmn.s(end,Fmn.n0ind)) marks whether Nu is even
  Nu=(size(Fmn.c,1)-1)*2+not(isnan(Fmn.s(end,n0ind_tmp)));
  Nudiscr=Nu; %used for the assembly below
  Nvdiscr=Nv;  %used for the assembly below
  
  if nargin==4
    Nu_in=Nu;
    Nv_in=Nv;
    Nu=varargin{2};
    Nv=varargin{3};
    Nudiscr=Nu;
    Nvdiscr=Nv;
    if mod(Nu-Nu_in,2)~=0 || mod(Nv-Nv_in,2)~=0
      error('The new M must have the same even/odd parity as the old M. The same goes for N!')
    end

    mmax=floor(Nu/2);
    nmax=floor(Nv/2);
    Nu_even=not(mod(Nu,2));
    Nv_even=not(mod(Nv,2));
    nrange=-nmax+Nv_even:nmax;
    mrange=0:mmax;
    %nrange=-Nv/2+1:Nv/2;
    %mrange=0:Nu/2;
    [m,n]=ndgrid(mrange,nrange);
    m0ind=1;
    n0ind=ceil(Nv/2);
    
    Fmn_in=Fmn;
    if Nu>Nu_in && Nv>Nv_in
      du=(Nu-Nu_in)/2;
      dv=(Nv-Nv_in)/2;
      m0ind_in=1;
      n0ind_in=ceil(Nv_in/2);
      Fmn.c=zeros(mmax+1,Nv);
      Fmn.s=zeros(mmax+1,Nv);
      if Nu_even
        Fmn_in.c(end,1:n0ind_in-1)=0; %remove these NaNs before enlarging
        Fmn_in.s(end,1:n0ind_in)=0; %remove these NaNs before enlarging
        if Nv_even
          Fmn_in.s(end,end)=0; %remove this NaN before enlarging
        end
      end
      if Nv_even
        Fmn_in.s(1,end)=0; %remove this NaN before enlarging
      end
      
      Fmn.c(1:end-du,1+dv:end-dv)=Fmn_in.c;
      Fmn.s(1:end-du,1+dv:end-dv)=Fmn_in.s;
      
      Fmn.c(1,1:n0ind-1)=NaN;
      Fmn.s(1,1:n0ind)=NaN;
      if Nu_even
        Fmn.c(end,1:n0ind-1)=NaN;
        Fmn.s(end,1:n0ind)=NaN;
        if Nv_even
          Fmn.s(end,end)=NaN;
        end
      end
      if Nv_even
        Fmn.s(1,end)=NaN;
      end     
    elseif Nu<=Nu_in && Nv<=Nv_in
      if (Nu<Nu_in || Nv<Nv_in) && not(listinput)
        disp(['ifftmn.m: Warning: Requested spatial discretisation ',...
              'may not resolve all modes!'])
      end
      if (Nu<Nu_in || Nv<Nv_in) && listinput
        disp(['ifftmn.m: Warning: Requested spatial discretisation is small. ',...
              'Resorting to listwise ifft. Aliasing of non-resolved high modes may appear.'])
        make_normal_assembly=1;
        %Go back to input representation
        Nudiscr=Nu;
        Nvdiscr=Nv;
        Fmn=Fmn_in;
        Nu=Nu_in;
        Nv=Nv_in;
        mmax=floor(Nu/2);
        nmax=floor(Nv/2);
        nrange=-nmax+Nv_even:nmax;
        mrange=0:mmax;
        [m,n]=ndgrid(mrange,nrange);
        m0ind=1;
        n0ind=ceil(Nv/2);        
      else
        du=(Nu_in-Nu)/2;
        dv=(Nv_in-Nv)/2;
        Fmn.c=Fmn_in.c(1:end-du,1+dv:end-dv);
        Fmn.s=Fmn_in.s(1:end-du,1+dv:end-dv);
        if Nu_even
          Fmn.c(end,1:n0ind-1)=NaN;
          Fmn.s(end,1:n0ind)=NaN;
          if Nv_even
            Fmn.s(end,end)=NaN;
          end
        end
        if Nv_even
          Fmn.s(1,end)=NaN;
        end
      end
    elseif Nu==Nu_in && Nv==Nv_in %(already covered by the previous case)
      Fmn=Fmn_in;
    else
      Nu
      Nu_in
      Nv
      Nv_in
      error('The case where Nu is increased and Nv decreased is not implemented yet!')
    end
  else
    mmax=floor(Nu/2);
    nmax=floor(Nv/2);
    Nu_even=not(mod(Nu,2));
    Nv_even=not(mod(Nv,2));
    nrange=-nmax+Nv_even:nmax;
    mrange=0:mmax;
    [m,n]=ndgrid(mrange,nrange);
    m0ind=1;
    n0ind=ceil(Nv/2);
  end
  
  uvec=(0:Nudiscr-1)/Nudiscr*2*pi;
  vvec=(0:Nvdiscr-1)/Nvdiscr*2*pi;
  [u,v]=ndgrid(uvec,vvec);

  %Normal assembly, is very slow
  if make_normal_assembly
    f=zeros(Nudiscr,Nvdiscr);
    for m=1:mmax-Nu_even
      %disp([num2str(m),' / ',num2str(mmaxNu-Nu_even)])
      for n= -nmax+Nv_even:nmax
        f=f+Fmn.c(m+1,n+n0ind)*cos(m*u-n*v)+Fmn.s(m+1,n+n0ind)*sin(m*u-n*v);
      end
    end
    
    for n=0:nmax
      f=f+Fmn.c(1,n+n0ind)*cos(-n*v)+...
        Nu_even*Fmn.c(end,n+n0ind)*cos(Nu/2*u-n*v);
    end
    
    for n=1:nmax-Nv_even
      f=f+Fmn.s(1,n+n0ind)*sin(-n*v)+...
        Nu_even*Fmn.s(end,n+n0ind)*sin(Nu/2*u-n*v);
    end
  else %ifft assembly, is much faster
    
    Fmn.c(m0ind,n0ind)=Fmn.c(m0ind,n0ind)*2;
    Fmn.c(m0ind,1:n0ind-1)=fliplr(Fmn.c(m0ind,n0ind+1:end-Nv_even));
    Fmn.s(m0ind,n0ind)=0;
    Fmn.s(m0ind,1:n0ind-1)=-fliplr(Fmn.s(m0ind,n0ind+1:end-Nv_even));
    if Nv_even
      Fmn.c(m0ind,end)=Fmn.c(m0ind,end)*2;    
      Fmn.s(m0ind,end)=0;
    end
    
    if Nu_even
      Fmn.c(end,n0ind)=Fmn.c(end,n0ind)*2;
      Fmn.c(end,1:n0ind-1)=fliplr(Fmn.c(end,n0ind+1:end-Nv_even));
      Fmn.s(end,n0ind)=0;
      Fmn.s(end,1:n0ind-1)=-fliplr(Fmn.s(end,n0ind+1:end-Nv_even));
      if Nv_even
        Fmn.c(end,end)=Fmn.c(end,end)*2;
        Fmn.s(end,end)=0;
      end
    end
    
    Fmncos=zeros(Nu,Nv);
    Fmnsin=zeros(Nu,Nv);
    Fmncos(1:mmax+1,:)=Fmn.c;
    Fmnsin(1:mmax+1,:)=Fmn.s;
    
    Fmncos(mmax+2:end,n0ind)=flipud(Fmncos(2:mmax+1-Nu_even,n0ind));
    Fmnsin(mmax+2:end,n0ind)=-flipud(Fmnsin(2:mmax+1-Nu_even,n0ind));
    
    if Nv_even
      Fmncos(mmax+2:end,end)  =flipud(Fmncos(2:mmax+1-Nu_even,end));
      Fmnsin(mmax+2:end,end)  =-flipud(Fmnsin(2:mmax+1-Nu_even,end));
    end
    
    Fmncos(mmax+2:end,1:n0ind-1)           ...
        =fliplr(flipud(Fmncos(2:mmax+1-Nu_even,n0ind+1:end-Nv_even)));
    Fmncos(mmax+2:end,n0ind+1:end-Nv_even) ...
        =fliplr(flipud(Fmncos(2:mmax+1-Nu_even,1:n0ind-1)));
    Fmnsin(mmax+2:end,1:n0ind-1)     ...
        =-fliplr(flipud(Fmnsin(2:mmax+1-Nu_even,n0ind+1:end-Nv_even)));
    Fmnsin(mmax+2:end,n0ind+1:end-Nv_even)       ...
        =-fliplr(flipud(Fmnsin(2:mmax+1-Nu_even,1:n0ind-1)));
    
    F=Nu*Nv/2*(Fmncos-i*Fmnsin);
    f=real(ifft2(ifftshift(fliplr(F),2)));
  end

  varargout{inputind}=f;
end

if nargout==Ninput+2
  varargout{Ninput+1}=u;
  varargout{Ninput+2}=v/Nperiods;
end