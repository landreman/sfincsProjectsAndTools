function summa=plus(A,B)

if compatible(A,B)
  summa.m=A.m;
  summa.n=A.n;
  summa.m0ind=A.m0ind;
  summa.n0ind=A.n0ind;
  summa.c=A.c+B.c;
  summa.s=A.s+B.s;
  summa=mnmat(summa);
else
  error('A and B not compatible')
end