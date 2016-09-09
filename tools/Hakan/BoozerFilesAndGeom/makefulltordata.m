function out=makefulltordata(in,Ntor,isTorCoord,varargin)
% This function is used by the flux surface plot routine CoordPlots
% to add one extra line poloidally and toroidally to a matrix "in"
% which has the dimensions (poloidal coord, toroidal coord). This is to make 
% the plotted surface close on itself.

%First index is theta and second is zeta

if nargin==2
  isTorCoord='';
end
if strcmp(isTorCoord,'toroidal coordinate') 
  torPeriod=2*pi/Ntor;
elseif strcmp(isTorCoord,'aligned toroidal coordinate')
  torPeriod=0;
else
  torPeriod=0;
end
if strcmp(isTorCoord,'poloidal coordinate')
  polPeriod=2*pi;
else
  polPeriod=0;
end
  

out = zeros(size(in,1),Ntor*size(in,2)+1);
for ind=0:Ntor-1
  out(:,1+ind*size(in,2):(ind+1)*size(in,2))=in+ind*torPeriod;
end
out(:,Ntor*size(in,2)+1)=in(:,1)+Ntor*torPeriod;

if strcmp(isTorCoord,'aligned toroidal coordinate') 
  %The discretisation should not be aligned, but the coordinate "in" should be.
  out(find(out(:,end)==0),end)=2*pi;
  iota=varargin{1};
  Deltazeta=2*pi/iota;
  vect0_2pi=linspace(0,2*pi,size(out,2));
  torv1=mod(vect0_2pi-Deltazeta,2*pi);
  outend=interp1(vect0_2pi,out(1,:),torv1);
  out=[out;outend];
else
  out=[out;out(1,:)+polPeriod]; %add an extra line to complete the poloidal turn
end