function [F,u,v,Nu,Nv,Du,Dv]=griddatacyclic(unoneq,vnoneq,Fnoneq,N)
% u must be 2*pi periodic.
% v must be 2*pi/N periodic.
% (unoneq, vnoneq) is the non-equidistant grid on which F is given.
% (u,v) is the equidistant grid generated with ndgrid which to interpolate to.

Nu=size(Fnoneq,1);
Nv=size(Fnoneq,2);

Dv=2*pi/Nv/N;
Du=2*pi/Nu;
vvec=(0:Nv-1)*Dv;
uvec=(0:Nu-1)'*Du;
[u,v] = ndgrid(uvec,vvec);

Mleftadd  = max(0,ceil(max(vnoneq(:,1))/(2*pi/N)));
Mrightadd = max(0,ceil(1-min(vnoneq(:,end))/(2*pi/N)));
Mv=Mleftadd+1+Mrightadd;

Mbottomadd  = max(0,ceil(max(unoneq(1,:))/(2*pi)));
Mtopadd     = max(0,ceil(1-min(unoneq(1,:))/(2*pi)));
Mu=Mtopadd+1+Mbottomadd;

unoneqBig=NaN*zeros(Nu*Mu,Nv*Mv);
vnoneqBig=NaN*zeros(Nu*Mu,Nv*Mv);
FnoneqBig=NaN*zeros(Nu*Mu,Nv*Mv);

for mui=1:Mu
  for mvi=1:Mv
    vadd=(mui-(Mleftadd+1))*2*pi/N;
    uadd=(mvi-(Mbottomadd+1))*2*pi;
    unoneqBig((mui-1)*Nu+1:mui*Nu, (mvi-1)*Nv+1:mvi*Nv)=unoneq+uadd;
    vnoneqBig((mui-1)*Nu+1:mui*Nu, (mvi-1)*Nv+1:mvi*Nv)=vnoneq+vadd;
    FnoneqBig((mui-1)*Nu+1:mui*Nu, (mvi-1)*Nv+1:mvi*Nv)=Fnoneq;
  end
end

F=griddata(unoneqBig,vnoneqBig,FnoneqBig,u,v);
