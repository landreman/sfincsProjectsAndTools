function out=invJacBdotgrad(arg,iota,Nperiods)

%Assumes Boozer coordinates

sz=size(arg);
if not(mod(sz(1),2)) || not(mod(sz(2),2))
  error('information loss if Nu or Nv even!')
end

out.m=arg.m;
out.n=arg.n;
out.m0ind=arg.m0ind;
out.n0ind=arg.n0ind;
out.c = -arg.s./(iota*arg.m-Nperiods*arg.n);
out.s = arg.c./(iota*arg.m-Nperiods*arg.n);

out.c(out.m0ind,out.n0ind)=0;
out.s(out.m0ind,out.n0ind)=NaN;

out=mnmat(out);

