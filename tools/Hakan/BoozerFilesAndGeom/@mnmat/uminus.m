function out=uminus(A)

out.m=A.m;
out.n=A.n;
out.m0ind=A.m0ind;
out.n0ind=A.n0ind;
out.c=-A.c;
out.s=-A.s;

out=mnmat(out);