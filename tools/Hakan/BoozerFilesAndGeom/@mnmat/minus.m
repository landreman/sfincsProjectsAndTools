function out=minus(A,B)

if compatible(A,B)
  out.m=A.m;
  out.n=A.n;
  out.m0ind=A.m0ind;
  out.n0ind=A.n0ind;
  out.c=A.c-B.c;
  out.s=A.s-B.s;
  out=mnmat(out);
else
  error('A and B not compatible')
end