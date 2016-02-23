function varargout=invgrad(Fmn,varargin)

Nu_even=isnan(Fmn.s(end,Fmn.n0ind));
 Nv_even=not(mod(size(Fmn.c,2),2));

if nargin==2
  %Only one argument given. Calculate two outputs whose u resp v derivatives equal Fmn.
  Nperiods=varargin{1};
  
  if get00(F)~=0
    error('There is a 00 component in an mnmat to be integrated!')
  end
  
  Fmndu.s =  Fmn.c./Fmn.m;
  Fmndu.c = -Fmn.s./Fmn.m;

  Fmndv.s = -Fmn.c./(Fmn.n*Nperiods);
  Fmndv.c =  Fmn.s./(Fmn.n*Nperiods);

  Fmndu.m=Fmn.m;
  Fmndv.m=Fmn.m;
  Fmndu.n=Fmn.n;
  Fmndv.n=Fmn.n;
  Fmndu.m0ind=Fmn.m0ind;
  Fmndv.m0ind=Fmn.m0ind;
  Fmndu.n0ind=Fmn.n0ind;
  Fmndv.n0ind=Fmn.n0ind;

  %Let us force the size of the data matrix to be the same
  %In principle one should increase it because there is more sine 
  %information than it can accommodate now, but that would be unpractical
  %if we want to transform back to real space.

  Fmndu.c(Fmn.m0ind,Fmn.n0ind)=0;
  Fmndu.s(Fmn.m0ind,Fmn.n0ind)=NaN;
  if Nu_even
    Fmndu.c(end,Fmn.n0ind)=0;
    Fmndu.s(end,Fmn.n0ind)=NaN;
    if Nv_even
      Fmndu.c(end,end)=0;
      Fmndu.s(end,end)=NaN;
    end
  end
  if Nv_even
    Fmndu.c(Fmn.m0ind,end)=0;
    Fmndu.s(Fmn.m0ind,end)=NaN;
  end

  Fmndv.c(Fmn.m0ind,Fmn.n0ind)=0;
  Fmndv.s(Fmn.m0ind,Fmn.n0ind)=NaN;
  if Nu_even
    Fmndv.c(end,Fmn.n0ind)=0;
    Fmndv.s(end,Fmn.n0ind)=NaN;
    if Nv_even
      Fmndv.c(end,end)=0;
      Fmndv.s(end,end)=NaN;
    end
  end
  if Nv_even
    Fmndv.c(Fmn.m0ind,end)=0;
    Fmndv.s(Fmn.m0ind,end)=NaN;
  end
  Fmndu=mnmat(Fmndu);
  Fmndv=mnmat(Fmndv);

  if nargout==2
    varargout{1}=Fmndu;
    varargout{2}=Fmndv;
  elseif nargout==1
    varargout{1}=[Fmndu,Fmndv];
  end
  
elseif nargin>=3 %Both dduGmn and ddvGmn are given and Gmn is sought
  dduGmn=Fmn;
  ddvGmn=varargin{1};
  Nperiods=varargin{2};
  method=1;
  if nargin==4
    method=varargin{3};
  end
  
  if not(compatible(dduGmn,ddvGmn))
    error('Non-compatible inputs!')
  end

  Gmn.m=dduGmn.m;
  Gmn.n=dduGmn.n;
  Gmn.m0ind=dduGmn.m0ind;
  Gmn.n0ind=dduGmn.n0ind;

  Gmn.c=zeros(size(dduGmn.c));
  Gmn.s=zeros(size(dduGmn.s));
  
  nnon0=Gmn.n0ind+1:size(Gmn.c,2);
  
  %Two options that should be equivalent
  if method
    Gmn.s(2:end,:) =  dduGmn.c(2:end,:)./dduGmn.m(2:end,:);
    Gmn.c(2:end,:) = -dduGmn.s(2:end,:)./dduGmn.m(2:end,:);

    Gmn.s(1,nnon0) = -ddvGmn.c(1,nnon0)./(ddvGmn.n(1,nnon0)*Nperiods);
    Gmn.c(1,nnon0) =  ddvGmn.s(1,nnon0)./(ddvGmn.n(1,nnon0)*Nperiods);
  else
    Gmn.s(:,nnon0) = -ddvGmn.c(:,nnon0)./(ddvGmn.n(:,nnon0)*Nperiods);
    Gmn.c(:,nnon0) =  ddvGmn.s(:,nnon0)./(ddvGmn.n(:,nnon0)*Nperiods);
    
    Gmn.s(2:end,Gmn.n0ind) =  dduGmn.c(2:end,Gmn.n0ind)./dduGmn.m(2:end,Gmn.n0ind);
    Gmn.c(2:end,Gmn.n0ind) = -dduGmn.s(2:end,Gmn.n0ind)./dduGmn.m(2:end,Gmn.n0ind);    
  end
  
  Gmn.c(Gmn.m0ind,Gmn.n0ind)=0;
  Gmn.s(Gmn.m0ind,Gmn.n0ind)=NaN;
  if Nu_even
    Gmn.c(end,Gmn.n0ind)=0;
    Gmn.s(end,Gmn.n0ind)=NaN;
    if Nv_even
      Gmn.c(end,end)=0;
      Gmn.s(end,end)=NaN;
    end
  end
  if Nv_even
    Gmn.c(Gmn.m0ind,end)=0;
    Gmn.s(Gmn.m0ind,end)=NaN;
  end
  
  Gmn=mnmat(Gmn);
  
  if nargout==1
    varargout{1}=Gmn;
  else
    error('Too many outputs!')
  end
end