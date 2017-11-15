function result=mrdivide(A,B)
% the .* opteration
%various ways of taking a product in the spatial domain

if isa(B,'double')
  if length(B)==1
    result=A;
    result.c=A.c/B;
    result.s=A.s/B;
  elseif all(size(A)==size(B))
    result=fftmn(ifftmn(A)./B);
  else
    size_A=size(A)
    size_B=size(B)
    error('A and B cannot be multiplied')
  end
elseif  isa(B,'mnmat')
  if compatible(A,B)
    result=fftmn(ifftmn(A)./ifftmn(B));
  else
    error('A and B not compatible!')
  end 
  warning('Division Amn/Bmn was done by inverse transforming, and then transforming A/B back!')
else
  error('not implemented case')
end
  