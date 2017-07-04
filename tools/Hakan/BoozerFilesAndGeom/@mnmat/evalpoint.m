function vals=evalpoint(Fmns,u,v,Nperiods)

Ninput=length(Fmns);
vals=zeros(1,Ninput);

if Ninput>1
  Fmn2=Fmns(1);
  for inputind=2:Ninput
    Fmn1=Fmn2;
    Fmn2=Fmns(inputind);
    if not(compatible(Fmn1,Fmn2))
      error('The input must have the same format!')
    end
  end
end


make_lists_first=0;
if make_lists_first
  lists=mnlist(Fmns,-eps,'absolute');
  if Ninput>1
    for ind=1:length(lists(1).m)
      c=cos(lists(1).m(ind) * u - lists(1).n(ind) * Nperiods * v);
      s=sin(lists(1).m(ind) * u - lists(1).n(ind) * Nperiods * v);
      if lists(1).cosparity(ind)
        for lind=1:Ninput
          vals(lind)=vals(lind)+lists(lind).data(ind) * c;
        end
      else
        for lind=1:Ninput
          vals(lind)=vals(lind)+lists(lind).data(ind) * s;
        end      
      end
    end
  else
    for ind=1:length(lists.m)
      c=cos(lists.m(ind) * theta - lists.n(ind) * Nperiods * zeta);
      s=sin(lists.m(ind) * theta - lists.n(ind) * Nperiods * zeta);
      if lists.cosparity(ind)
        vals=vals+lists.data(ind) * c;
      else
        vals=vals+lists.data(ind) * s;
      end
    end
  end
else
  if Ninput>1
    c=cos(Fmns(1).m * u - Fmns(1).n * Nperiods * v);
    s=sin(Fmns(1).m * u - Fmns(1).n * Nperiods * v);
    cnan=find(isnan(Fmns(1).c));
    snan=find(isnan(Fmns(1).s));
    for lind=1:Ninput
      Fmns(lind).c(cnan)=0;
      Fmns(lind).s(snan)=0;
      vals(lind)=sum(sum(Fmns(lind).c.*c+...
                         Fmns(lind).s.*s));
    end
  else
    if length(u)==1 && length(v)==1    
      c=cos(Fmns.m * u - Fmns.n * Nperiods * v);
      s=sin(Fmns.m * u - Fmns.n * Nperiods * v);
      cnan=find(isnan(Fmns.c));
      snan=find(isnan(Fmns.s));
      Fmns.c(cnan)=0;
      Fmns.s(snan)=0;
      vals=sum(sum(Fmns.c.*c+Fmns.s.*s));
    else %the case of vector input u and/or v
      su=size(u);sv=size(v);
      if su(1)>su(2)
        u=u';
      end
      if sv(1)>sv(2)
        v=v';
      end
      c=cos(Fmns.m(:) * u - Fmns.n(:) * v * Nperiods);
      s=sin(Fmns.m(:) * u - Fmns.n(:) * v * Nperiods);
      cnan=find(isnan(Fmns.c));
      snan=find(isnan(Fmns.s));
      Fmns.c(cnan)=0;
      Fmns.s(snan)=0;
      vals=Fmns.c(:)'*c+Fmns.s(:)'*s;
    end
  end
end