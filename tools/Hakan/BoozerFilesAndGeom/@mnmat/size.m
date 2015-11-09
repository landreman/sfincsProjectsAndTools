function sz=size(A)

if length(A)==1
  mmaxp1=size(A.m,1);
  isodd=not(isnan(A.s(A.m0ind,end)));
  sz=[(mmaxp1-1)*2+isodd,size(A.m,2)];
else
  sz=[numel(A),size(A(1))];
end