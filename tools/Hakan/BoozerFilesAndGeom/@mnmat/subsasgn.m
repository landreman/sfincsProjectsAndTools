function A=subsasgn(Ain,S,B)

A=Ain;
if strcmp(S(1).type,'{}')
  if not iscell(A)
    error('Indexing an mnmat with {} brackets is undefined!')
  else
    if length(S)==1
      if length(S(1).subs)==1
        A{S(1).subs{1}}=B;
      elseif length(S(1).subs)==2
        A{S(1).subs{1},S(1).subs{2}}=B;
      else
        error('3D tensors of mnmats are not supported')
      end
    else
      error('indexing like Bmn{3}.cos(3,-2) is not supported')
    end
  end
elseif strcmp(S(1).type,'()')
  if length(S(1).subs)>1
    error('Matrices of mnmats are not supported!')
  end
  if length(S)==1
    A(S(1).subs{1})=B;
  else
    error('indexing like Amns(3).c=... is not supported')
  end
else %it starts witn a '.'
  if length(S)==1
    if strcmp(S.subs,'c')
      out=A.c;
    elseif strcmp(S(1).subs,'s')
      out=A.s;
    else
      error(['Indexing like Amn.',S.subs,'=... is not supported'])
    end
  else
    if length(S)~=2
      error(['The reference to elements of an mnmat is made with the notation ',...
             'A.c(mind,nind), A.s(mind,nind), A.cos(m,n) or A.sin(m,n). (The keyword ',...
              '''end'' cannot be used with .cos or .sin.)'])
    end
    if ~strcmp(S(2).type,'()')
      error(['The reference to elements of an mnmat is made with the notation ',...
             'A.c(mind,nind), A.s(mind,nind), A.cos(m,n) or A.sin(m,n).'])
    end

    if length(S(2).subs)~=2
      error(['The reference to elements of an mnmat with the notation ',...
             'A.c(mind,nind), A.s(mind,nind), A.cos(m,n) or A.sin(m,n) needs two input arguments.'])
    end
    %Now we know that length(S.subs)=2.
 
    if strcmp(S(1).subs,'c') || strcmp(S(1).subs,'s')
      minds=S(2).subs{1};
    elseif strcmp(S(1).subs,'cos') || strcmp(S(1).subs,'sin')
      if isnumeric(S(2).subs{1})
        minds=A.m0ind+S(2).subs{1};
      else
        error('The input m in the notation A.cos(m,n) or A.sin(m,n) must be numeric')
      end
    end

    if strcmp(S(1).subs,'c') || strcmp(S(1).subs,'s')
      ninds=S(2).subs{2};
    elseif strcmp(S(1).subs,'cos') || strcmp(S(1).subs,'sin')
      if isnumeric(S(2).subs{2})
        ninds=A.n0ind+S(2).subs{2};
      else
        error('The input n in the notation A.cos(m,n) or A.sin(m,n) must be numeric')
      end
    end
    
    if strcmp(S(1).subs,'cos') || strcmp(S(1).subs,'sin')
      if (any(minds<1) || any(minds>size(A.c,1))) && (any(ninds<1) || any(ninds>size(A.c,2)))
        error(['The inputs m and n are both out of range! Allowed is 0<=m<=',...
               num2str(size(A.c,1)-1),', ',...
               num2str(1-A.n0ind),'<=n<=',num2str(size(A.c,2)-A.n0ind)])
      end
      if any(minds<1) || any(minds>size(A.c,1)) 
        error(['The input m is out of range! Allowed is 0<=m<=',num2str(size(A.c,1)-1)])
      end
      if any(ninds<1) || any(ninds>size(A.c,2)) 
        error(['The input n is out of range! Allowed is ',num2str(1-A.n0ind),'<=n<=',num2str(size(A.c,2)-A.n0ind)])
      end
    end
    
    Snew.type='()';
    Snew.subs{1}=minds;
    Snew.subs{2}=ninds;
    if strcmp(S(1).subs,'cos') || strcmp(S(1).subs,'c')
      if all(all(isnan(subsref(A.c,Snew))==isnan(B)))
        A.c=subsasgn(A.c,Snew,B); %A.c(minds,ninds)=B;
      else
        error('Assignment not allowed! NaN pattern does not coincide with original.')
      end
    elseif strcmp(S(1).subs,'sin') || strcmp(S(1).subs,'s')
      if all(all(isnan(subsref(A.s,Snew))==isnan(B)))
        A.s=subsasgn(A.s,Snew,B); %A.s(minds,ninds)=B;
      else
        error('Assignment not allowed! NaN pattern does not coincide with original.')
      end
    end

  
  end
end