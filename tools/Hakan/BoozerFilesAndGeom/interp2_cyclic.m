function Fout=interp2_cyclic(u,v,F,uout,vout,N)

% u must be 2*pi periodic.
% v must be 2*pi/N periodic.

vDirection=sign(v(1,end)-v(1,1));
uDirection=sign(u(end,1)-u(1,1));

vbig=[v,v(:,1)+vDirection*2*pi/N];
vbig=[vbig;vbig(1,:)];
ubig=[u;u(1,:)+uDirection*2*pi];
ubig=[ubig,ubig(:,1)];
Fbig=[F,F(:,1)];
Fbig=[Fbig;Fbig(1,:)];

Fout=interp2(vbig,ubig,Fbig,mod(vout,2*pi/N),mod(uout,2*pi));