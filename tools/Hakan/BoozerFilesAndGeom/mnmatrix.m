function Fmns=mnmatrix(lsts,varargin)
securityLevel=0; %-1: no security, 0: warnings, 1:errors

for inputind=1:length(lsts)
  
  lst=lsts(inputind);

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
      for pind=problems
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


  mmax=max(lst.m)+1; %add 1 not to loose sign information on the highest m components
  nmax=max(abs(lst.n))+1; %same for n

  %Nu=2^(ceil(log(max(mmax,nmax))/log(2))+1); %If we want powers of 2
  Nu=ceil(max(mmax,nmax))*2;
  Nv=Nu;

  if nargin==3
    Nu=max(varargin{1},Nu);
    Nv=max(varargin{2},Nv);
  end
  
  inds=find(lst.n==-Nv/2);
  if not(isempty(inds))
    lst.n(inds)=-lst.n(inds);
  end

  problem=0;
  ind1=find(lst.m==0 & lst.n==Nv/2);
  ind2=find(lst.m==Nu/2 & lst.n==Nv/2);
  ind3=find(lst.m==Nu/2 & lst.n==0);
  inds4=find(lst.m==Nu/2 & lst.n<0);
  if not(isempty(ind1))
    if lst.cosparity(ind1)==0
      problem=1;
    end
  end
  if not(isempty(ind2))
    if lst.cosparity(ind2)==0
      problem=1;
    end
  end
  if not(isempty(ind3))
    if lst.cosparity(ind3)==0
      problem=1;
    end
  end
  if not(isempty(inds4))
    problem=1;
  end

  if problem
    Nv=Nv+2;
    Nu=Nu+2;
  end

  nrange=-Nv/2+1:Nv/2;
  mrange=0:Nu/2;
  [m,n]=ndgrid(mrange,nrange);
  Fmn.c=zeros(size(m));
  Fmn.s=zeros(size(m));
  Fmn.m=m;
  Fmn.n=n;
  Fmn.m0ind=1;
  Fmn.n0ind=Nv/2;

  mind=lst.m+1;
  nind=lst.n+Nv/2;

  for li=1:length(lst.m)
    if lst.cosparity(li)
      Fmn.c(mind(li),nind(li))=lst.data(li);
    else
      Fmn.s(mind(li),nind(li))=lst.data(li);
    end
  end

  Fmn.c(Fmn.m0ind,1:Fmn.n0ind-1)=NaN;
  Fmn.s(Fmn.m0ind,1:Fmn.n0ind)=NaN;
  Fmn.c(end,1:Fmn.n0ind-1)=NaN;
  Fmn.s(end,1:Fmn.n0ind)=NaN;
  Fmn.s(Fmn.m0ind,end)=NaN;
  Fmn.s(end,end)=NaN;

  %Store the output
  Fmns(inputind)=Fmn;
end