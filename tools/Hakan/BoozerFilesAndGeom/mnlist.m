function Fmnlists=mnlist(Fmns,varargin)
securityLevel=-1; %-1: no security, 0: warnings, 1:errors

if nargin==1
  threshold=0;
  thresholdType='absolute';
elseif nargin==2
  threshold=varargin{1};
  thresholdType='relative';
elseif nargin==3
  threshold=varargin{1};
  thresholdType=varargin{2};
end

for inputind=1:length(Fmns)
  Fmn=Fmns(inputind);

  Nu_even=isnan(Fmn.s(end,Fmn.n0ind));
  Nv_even=not(mod(size(Fmn.c,2),2));
  
  securityLevel=1;
  if not(isnan(Fmn.s(Fmn.m0ind,Fmn.n0ind)))
    error('A m=n=0 sine term is present!')
  end
  if Nv_even
    if not(isnan(Fmn.s(Fmn.m0ind,end)))
      error('Unresolvable sine information included!')
    end
    if Nu_even
      if not(isnan(Fmn.s(end,end)))
        error('Unresolvable sine information included!')
      end
    end
  end
  if not(all(isnan(Fmn.s(Fmn.m0ind,1:Fmn.n0ind-1))))
    if securityLevel<=0
      if securityLevel==0
        warning(['Negative n values for m=0 sine terms found! Subtracting these components ', ...
                 'from the corresponding positive n data'])
      end
      Fmn.s(Fmn.m0ind,Fmn.n0ind+1:end-1)=Fmn.s(Fmn.m0ind,Fmn.n0ind+1:end-1) ...
          -fliplr(Fmn.s(Fmn.m0ind,2:Fmn.n0ind-1));
    else % securityLevel==1
      error('Negative n values for m=0 sine terms found!')
    end
  end
  if not(all(isnan(Fmn.c(Fmn.m0ind,1:Fmn.n0ind-1))))
    if securityLevel<=0
      if securityLevel==0
        warning(['Negative n values for m=0 cos terms found! Adding these components ', ...
                 'to the corresponding positive n data'])
      end
      Fmn.c(Fmn.m0ind,Fmn.n0ind+1:end-1)=Fmn.c(Fmn.m0ind,Fmn.n0ind+1:end-1) ...
          +fliplr(Fmn.c(Fmn.m0ind,2:Fmn.n0ind-1));
    else % securityLevel==1
      error('Negative n values for m=0 cos terms found!')
    end
  end
  if Nu_even
    if not(all(isnan(Fmn.s(end,1:Fmn.n0ind-1))))
      if securityLevel==0
        warning(['Negative n values for m=mmax sine terms found! Subtracting these components ', ...
                 'from the corresponding positive n data'])
        Fmn.s(end,Fmn.n0ind+1:end-1)=Fmn.s(end,Fmn.n0ind+1:end-1) ...
            -fliplr(Fmn.s(end,2:Fmn.n0ind-1));
      else % securityLevel==1
        error('Negative n values for m=mmax sine terms found!')
      end
    end
    if not(all(isnan(Fmn.c(end,1:Fmn.n0ind-1))))
      if securityLevel==0
        warning(['Negative n values for m=mmax cos terms found! Adding these components ', ...
                 'to the corresponding positive n data'])
        Fmn.c(end,Fmn.n0ind+1:end-1)=Fmn.c(end,Fmn.n0ind+1:end-1) ...
            +fliplr(Fmn.c(end,2:Fmn.n0ind-1));
      else % securityLevel==1
        error('Negative n values for m=mmax cos terms found!')
      end
    end
  end

  if strcmp(thresholdType,'relative')
    threshold=threshold*Fmn.c(Fmn.m0ind,Fmn.n0ind);
  elseif not(strcmp(thresholdType,'absolute'))
    error('thresholdType not recognized!')
  end

  cinds=find(abs(Fmn.c)>threshold)';
  sinds=find(abs(Fmn.s)>threshold)';

  Fmnlist.data=[Fmn.c(cinds),Fmn.s(sinds)];
  Fmnlist.m=[Fmn.m(cinds),Fmn.m(sinds)];
  Fmnlist.n=[Fmn.n(cinds),Fmn.n(sinds)];
  Fmnlist.cosparity=[ones(1,length(cinds)),zeros(1,length(sinds))];

  %Store the output
  Fmnlists(inputind)=Fmnlist;
end