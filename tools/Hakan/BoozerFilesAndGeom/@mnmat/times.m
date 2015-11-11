function prod=times(A,B)
% the .* opteration
%various ways of taking a product in the spatial domain

if isa(B,'double')
  if length(B)==1
    prod=A;
    prod.c=A.c*B;
    prod.s=A.s*B;
  elseif all(size(A)==size(B))
    prod=fftmn(ifftmn(A).*B);
  else
    size_A=size(A)
    size_B=size(B)
    error('A and B cannot be multiplied')
  end
elseif  isa(B,'mnmat')
  %convolution
  if compatible(A,B)
    prod=fftmn(ifftmn(A).*ifftmn(B));
  else
    error('A and B not compatible!')
  end
else
  error('not implemented case')
end
  