function varargout=mnmat(lsts,varargin)
securityLevel=-1; %-1: no security, 0: warnings, 1:errors

for inputind=1:length(lsts)
  lst=lsts(inputind);
  clear Fmn
  
  if strcmp(class(lst),'mnmat')
    Fmn=lst;
  elseif not(isfield(lst,'cosparity')) &&...
        isfield(lst,'m') && ...
        isfield(lst,'n') && ...
        isfield(lst,'m0ind') && ...
        isfield(lst,'n0ind') && ...
        isfield(lst,'c') && ...
        isfield(lst,'s')
    Fmn.c=lst.c;
    Fmn.s=lst.s;
    Fmn.m=lst.m;
    Fmn.n=lst.n;
    Fmn.m0ind=lst.m0ind;
    Fmn.n0ind=lst.n0ind;

    Fmn=class(Fmn,'mnmat');
  elseif isfield(lst,'headertext')&&...
        isfield(lst,'StelSym')&&...
        isfield(lst,'Nperiods')&&...
        isfield(lst,'Bmn')&&...
        isfield(lst,'R')&&...
        isfield(lst,'Z')&&...
        isfield(lst,'Dphi')&&...
        isfield(lst,'m')&&...
        isfield(lst,'n')&&...
        isfield(lst,'rnorm')&&...
        isfield(lst,'s')&&...
        isfield(lst,'psi_a')
    
    %disp('mnmat input from a Geom struct.')
    %This is a Geometry struct
    %I require the input (Geom,rind,fieldtoextract)
    %or (Geom,rind,fieldtoextract,Nu,Nv,forceSize_option)
    if nargin<3
      error('Not enough inputs')
    end
    rind=varargin{1};
    objecttoextract=varargin{2};
    if strcmp(objecttoextract,'Bmn') || strcmp(objecttoextract,'B')
      fieldtoextract='Bmn';
    elseif strcmp(objecttoextract,'R')
      fieldtoextract='R';
    elseif strcmp(objecttoextract,'Z')
      fieldtoextract='Z';
    elseif strcmp(objecttoextract,'Dzetacylphi') || ...
          strcmp(objecttoextract,'Dphi')
      %Dzetacylphi=2*pi/Nperiods*Dphi
      fieldtoextract='Dphi';
    else
      fieldtoextract=objecttoextract;
    end    
    clear lista
    lista.m=lst.m{rind};
    lista.n=lst.n{rind};
    if lst.StelSym
      parity=ones(size(lista.m));
    else
      parity=lst.parity{rind};
    end
    if strcmp(fieldtoextract,'Z') || strcmp(fieldtoextract,'Dphi')
      lista.cosparity=not(parity); %only sin for StelSym 
    else
      lista.cosparity=parity; %only cos for StelSym 
    end
    tmp=getfield(lst,fieldtoextract);
    if strcmp(fieldtoextract,'Dphi')
      lista.data=2*pi/getfield(lst,'Nperiods')*tmp{rind};
    else
      lista.data=tmp{rind};
    end
    if nargin==3
      Fmn=mnmat(lista); %recursive call to default case below  
    elseif nargin==5
      Fmn=mnmat(lista,varargin{3},varargin{4});
    elseif nargin==6
      Fmn=mnmat(lista,varargin{3},varargin{4},varargin{5});
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else %Default case: list input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(lst.m,1)>size(lst.n,2)
      lst.m=lst.m';
    end
    if size(lst.n,1)>size(lst.n,2)
      lst.n=lst.n';
    end
    if size(lst.data,1)>size(lst.data,2)
      lst.data=lst.data';
    end
    if size(lst.cosparity,1)>size(lst.cosparity,2)
      lst.cosparity=lst.cosparity';
    end
    
    if any(lst.m<0)
      error('Negative m values found!')
    end
    if any(lst.n(find(lst.m==0))<0)
      if securityLevel<=0
        if securityLevel==0   
          disp(['Warning: Negative n values for m=0 found! ',...
                'Combining these components with the corresponding positive n data'])
        end
        m0inds=find(lst.m==0);
        problems=m0inds(lst.n(m0inds)<0);
        for pindi=1:length(problems)
          pind=problems(pindi);
          par=lst.cosparity(pind);
          n  =lst.n(pind);   
          dat=lst.data(pind);
          opposite=find(lst.cosparity==par & lst.n==-n & lst.m==0);
          if not(isempty(opposite))
            lst.data(opposite)=lst.data(opposite)+(-1)^(par+1)*dat;
          else
            lst.m=[lst.m,0];
            lst.n=[lst.n,-n];
            lst.cosparity=[lst.cosparity,par];
            lst.data=[lst.data,(-1)^(par+1)*dat];
          end
        end
      else
        error('Negative n values for m=0 found!')
      end
    end
    m0n0inds=find(lst.m==0 & lst.n==0);
    if not(isempty(m0n0inds))
      if any(lst.cosparity(m0n0inds)==0)
        pind=m0n0inds(find(lst.cosparity(m0n0inds)==0));
        if lst.data(pind)~=0
          if securityLevel<=0
            if securityLevel==0
              disp('Warning: A m=n=0 sine term is present! It will be ignored.')
            end
            lst.m   =[lst.m(1:pind-1),   lst.m(pind+1:end)];
            lst.n   =[lst.n(1:pind-1),   lst.n(pind+1:end)];
            lst.data=[lst.data(1:pind-1),lst.data(pind+1:end)];
            lst.cosparity=[lst.cosparity(1:pind-1),lst.cosparity(pind+1:end)];
          else
            error('A m=n=0 sine term is present!')
          end
        end
    end
    end
    
    if nargin>=3 %We know if the user wants even/odd Nu and Nv
      Nu_even=not(mod(varargin{1},2));
      Nv_even=not(mod(varargin{2},2));
    else
      Nu_even=1;%Default
      Nv_even=1;%Default
    end
    
    forceSize=0;
    if nargin==4
      forceSize=strcmp(varargin{3},'forceSize');
    end
    
    if forceSize
      Nu=varargin{1};
      Nv=varargin{2};
    else
      mmax=max(lst.m)+1; %add 1 not to loose sign information on the highest m components
      nmax=max(abs(lst.n))+1; %same for n
      
      nmax=nmax + (not(mod(nmax,2))~=Nv_even);
      Nv=2*nmax+not(Nv_even);
      Nu=2*mmax+not(Nu_even);
      
      if nargin>=3
        Nu=max(varargin{1},Nu);
        Nv=max(varargin{2},Nv);
      end
      
      if Nv_even
        inds=find(lst.n==-Nv/2);
        if not(isempty(inds))
          lst.n(inds)=-lst.n(inds);
        end
      end
      
      problem=0;
      ind1=[];ind2=[];ind3=[];inds4=[];
      if Nv_even
        ind1=find(lst.m==0 & lst.n==Nv/2);
        if Nu_even
          ind2=find(lst.m==Nu/2 & lst.n==Nv/2);
        end
      end
      if Nu_even
        ind3=find(lst.m==Nu/2 & lst.n==0);
        inds4=find(lst.m==Nu/2 & lst.n<0);
      end
      if not(isempty(ind1))
        if lst.cosparity(ind1)==0
          problem=1;
      end
      end
      if not(isempty(ind2))
        if lst.cosparity(ind2)==0
          problem=2;
        end
      end
      if not(isempty(ind3))
        if lst.cosparity(ind3)==0
          problem=3;
        end
      end
      if not(isempty(inds4))
        problem=3;
      end
      
      switch problem
       case 1
        Nv=Nv+2;
       case 2
        Nv=Nv+2;
        Nu=Nu+2;
       case 3
        Nu=Nu+2;
      end
    end
    
    mmax=floor(Nu/2);
    nmax=floor(Nv/2);
    nrange=-nmax+Nv_even:nmax;
    mrange=0:mmax;
    
    [m,n]=ndgrid(mrange,nrange);
    Fmn.c=zeros(size(m));
    Fmn.s=zeros(size(m));
    Fmn.m=m;
    Fmn.n=n;
    Fmn.m0ind=1;
    Fmn.n0ind=ceil(Nv/2);
    
    mind=lst.m+1;
    nind=lst.n+Fmn.n0ind;
    
    
    for li=1:length(lst.m)
      if lst.m(li)<=mmax && lst.n(li)<=nmax && lst.n(li)>=-nmax+Nv_even
        if lst.cosparity(li)
          Fmn.c(mind(li),nind(li))=lst.data(li);
        else
          Fmn.s(mind(li),nind(li))=lst.data(li);
        end
      end
    end
    
    Fmn.c(Fmn.m0ind,1:Fmn.n0ind-1)=NaN;
    Fmn.s(Fmn.m0ind,1:Fmn.n0ind)=NaN;
    if Nu_even
      Fmn.c(end,1:Fmn.n0ind-1)=NaN;
      Fmn.s(end,1:Fmn.n0ind)=NaN;
      if Nv_even
      Fmn.s(end,end)=NaN;
      end
    end
    if Nv_even
      Fmn.s(Fmn.m0ind,end)=NaN;
    end   
    Fmn=class(Fmn,'mnmat');
    
  end
  
  %Store the output
  if nargout==length(lsts)
    varargout{inputind}=Fmn;
  elseif nargout==1 || nargout==0
    varargout{1}(inputind)=Fmn;
  else
    error('Incorrect number of output!')
  end
end