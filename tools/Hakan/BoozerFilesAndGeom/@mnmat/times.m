function prod=times(A,B)
% the .* opteration
%various ways of taking a product in the spatial domain

if isa(A,'mnmat') && isa(B,'mnmat')
  %convolution
  if compatible(A,B)
    prod=fftmn(ifftmn(A).*ifftmn(B));
  else
    error('A and B not compatible!')
  end
elseif isa(B,'double')
  if length(B)==1
    prod=A;
    prod.c=A.c*B;
    prod.s=A.s*B;
  elseif all(size(A)==size(B))
    prod=fftmn(ifftmn(A).*B);
  end
else
  error('not implemented case')
end
  