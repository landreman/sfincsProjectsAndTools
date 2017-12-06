function out=subsref(A,S)
%
% This function is used to retrieve single m,n components of an mnmat A.
% Usage:
% A.cos(m,n) or A.sin(m,n)
%
% It is also called when indexing a cell or a vector of mnmats
%

if strcmp(S(1).type,'{}')
  if not iscell(A)
    error('Indexing an mnmat with {} brackets is undefined!')
  else
    if length(S(1).subs)==1
      out=A{S(1).subs{1}};
    elseif length(S(1).subs)==2
      out=A{S(1).subs{1},S(1).subs{1}};
    else
      error('3D tensors of mnmats are not supported')
    end
  end
elseif strcmp(S(1).type,'()')
  if length(S(1).subs)>1
    error('Matrices of mnmats are not supported!')
  end
  if length(S)==1
    out=A(S(1).subs{1});
  else
    A=A(S(1).subs{1});
    S=S(2:end);
    out=subsref(A,S);
  end
else

  if length(S)~=2
    error(['The reference to elements of an mnmat is made with the notation ',...
           'A.cos(m,n) or A.sin(m,n). The keyword ''end'' cannot be used.'])
  end
  if ~strcmp(S(2).type,'()')
    error(['The reference to elements of an mnmat is made with the notation ',...
           'A.cos(m,n) or A.sin(m,n).'])
  end

  if length(S(2).subs)~=2
    error(['The reference to elements of an mnmat with the notation ',...
           'A.cos(m,n) or A.sin(m,n) needs two input arguments, m and n.'])
  end
  %Now we know that length(S.subs)=2.
  
  if isnumeric(S(2).subs{1})
    minds=A.m0ind+S(2).subs{1};
  else
    error('The input m in the notation A.cos(m,n) or A.sin(m,n) must be numeric')
  end

  if isnumeric(S(2).subs{2})
    ninds=A.n0ind+S(2).subs{2};
  else
    error('The input n in the notation A.cos(m,n) or A.sin(m,n) must be numeric')
  end
  
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
  
  if strcmp(S(1).subs,'cos')
    out=A.c(minds,ninds);
  elseif strcmp(S(1).subs,'sin')
    out=A.s(minds,ninds);
  end

end