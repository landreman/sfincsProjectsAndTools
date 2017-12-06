function sz=size(A)
%
% This function returns the size of the corresponding spatial discretisation for A
%

if length(A)==1
  mmax=size(A.m,1)-1;
  %isodd=not(isnan(A.s(A.m0ind,end))); %old wrong
  isodd=not(isnan(A.s(end,A.n0ind))); %corrected 20171205
  sz=[mmax*2+isodd,size(A.m,2)];
else
  sz=[numel(A),size(A(1))];
end
