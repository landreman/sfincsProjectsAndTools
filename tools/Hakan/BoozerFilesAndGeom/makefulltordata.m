function out=makefulltordata(in,Ntor,isTorCoord)
% This function is used by the flux surface plot routine CoordPlots
% to add one extra line poloidally and toroidally to a matrix "in"
% which has the dimensions (pol. coord, tor. coord). This is to make 
% the plotted surface close on itself.

%First index is theta and second is zeta

if nargin==2
  isTorCoord='';
end
if strcmp(isTorCoord,'toroidal coordinate')
  %torPeriod=2*pi/Ntor*sign(in(1,end)-in(1,1)); %risky
  torPeriod=2*pi/Ntor;
else
  torPeriod=0;
end
if strcmp(isTorCoord,'poloidal coordinate')
  %polPeriod=2*pi*sign(in(end,1)-in(1,1)); %risky
  polPeriod=2*pi;
else
  polPeriod=0;
end
  
in=[in;in(1,:)+polPeriod]; %add an extra line to complete the poloidal turn

out = zeros(size(in,1),Ntor*size(in,2)+1);
for ind=0:Ntor-1
  out(:,1+ind*size(in,2):(ind+1)*size(in,2))=in+ind*torPeriod;
end
out(:,Ntor*size(in,2)+1)=in(:,1)+Ntor*torPeriod;