function plotonfluxsurf(Discr,toplot,varargin)

if nargin==2
  samp_pol=1;
  samp_tor=1;
elseif nargin==3
  samp_pol=varargin{1}; %For subsampling
  samp_tor=varargin{1}; %For subsampling
elseif nargin==4
  samp_pol=varargin{1}; %For subsampling
  samp_tor=varargin{2}; %For subsampling
end

Nperiods=Discr.Nperiods;
cylphi=Discr.cylphi(1:samp_pol:end,1:samp_tor:end);
cylR=Discr.cylR(1:samp_pol:end,1:samp_tor:end);
cylZ=Discr.cylZ(1:samp_pol:end,1:samp_tor:end);

Ntheta=size(cylR,1);
Nzeta=size(cylR,2);

cylphif=makefulltordata(cylphi,Nperiods,'toroidal coordinate');

geomangf=-cylphif; %cylphi is minus the geometrical angle
Rf=makefulltordata(cylR,Nperiods);
Xf=Rf.*cos(geomangf);
Yf=Rf.*sin(geomangf);
Zf=makefulltordata(cylZ,Nperiods);

if ndims(toplot)==2

  toplotf=makefulltordata(toplot(1:samp_pol:end,1:samp_tor:end),Nperiods);
  
  shading flat
  surf(Xf,Yf,Zf,toplotf)
  shading flat
  axis equal

elseif ndims(toplot)==3
  if size(toplot,1)~=3
    error('I do not understand the format!')
  end
  
  vX=squeeze(toplot(1,1:samp_pol:end,1:samp_tor:end));
  vY=squeeze(toplot(2,1:samp_pol:end,1:samp_tor:end));
  vZ=squeeze(toplot(3,1:samp_pol:end,1:samp_tor:end));
  
  geomang=-cylphi;
  vR=vX.*cos(geomang)+vY.*sin(geomang);
  vg=-vX.*sin(geomang)+vY.*cos(geomang);
  
  vXf=zeros(size(Xf));
  vYf=zeros(size(Xf));
  vZf=zeros(size(Xf));
  
  vXf(:,1:Nzeta)=[vX;vX(1,:)];
  vYf(:,1:Nzeta)=[vY;vY(1,:)];
  vZf(:,1:Nzeta)=[vZ;vZ(1,:)];
  vXf(:,end)=vXf(:,1);
  vYf(:,end)=vYf(:,1);
  vZf(:,end)=vZf(:,1);
  
  for ind=1:Nperiods-1
    thisang=geomang+2*pi/Nperiods*ind;
    vXloc=vR.*cos(thisang)-vg.*sin(thisang);
    vYloc=vR.*sin(thisang)+vg.*cos(thisang);
    
    vXf(:,Nzeta*ind+1:Nzeta*(ind+1))=[vXloc;vXloc(1,:)];
    vYf(:,Nzeta*ind+1:Nzeta*(ind+1))=[vYloc;vYloc(1,:)];
    vZf(:,Nzeta*ind+1:Nzeta*(ind+1))=[vZ;vZ(1,:)];
  end
  
  quiver3(Xf,Yf,Zf,vXf,vYf,vZf)
  axis equal
  
else
  error('Strange number of dimensions to plot')
end