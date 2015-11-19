function prod=mtimes(A,B)
% the * opteration
% various ways of taking a product in the Fourier domain

if isa(A,'mnmat') && isa(B,'mnmat')
  %convolution
  if compatible(A,B)
    prod=A;
    prod.c=A.c.*B.c;
    prod.s=A.s.*B.s;
    warning(['This is a multiplication in Fourier space, '...
             'i.e. a convolution in real space!'])
  else
    error('A and B not compatible!')
  end
elseif isa(B,'double')
  if length(B)==1
    prod=A;
    prod.c=A.c*B;
    prod.s=A.s*B;
  elseif all(size(A)==size(B))
    Bmn=fftmn(B);
    if compatible(A,B)
      prod=A.*B;
    end
  end
elseif isa(A,'double') && isa(B,'mnmat')
  prod=mtimes(B,A); %The above case
else
  error('not implemented case')
end
  