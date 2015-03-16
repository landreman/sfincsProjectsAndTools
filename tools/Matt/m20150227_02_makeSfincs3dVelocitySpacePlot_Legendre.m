clear

% xMax is the maximum normalized speed used in each plane of the figure.
xMax = 3;

% N is the number of grid points used in each plane of the figure.
N = 300;
% A good value is 30-100 for testing things out, but you might want to go
% up to 300 for a "production" figure. This script starts to take a long
% time as this number increases.

filename = 'C:\Users\landreman\Documents\Ubuntu\sfincsOutput_velocitySpacePlotLegendres.h5';

% Pick which species to plot, which (export_f_theta,export_f_zeta) location to plot, and which iteration to plot:
ispecies = 1;
itheta=1;
izeta=1;
iteration=1;

whichf = 1;
% 0 = full f
% 1 = delta f

switch whichf
    case 0
        data = h5read(filename,'/full_f');
    case 1
        data = h5read(filename,'/delta_f');
    otherwise
        error('Invalid whichf')
end

theta = h5read(filename,'/export_f_theta');
zeta = h5read(filename,'/export_f_zeta');
export_f_xi_option = h5read(filename,'/export_f_xi_option');
if export_f_xi_option ~= 0
    error('export_f_xi_option must be 0')
end
x = h5read(filename,'/export_f_x');
Nspecies = h5read(filename,'/Nspecies');
Niterations = size(data,1);
Nxi = double(h5read(filename,'/Nxi'));

fprintf('Number of available iterations: %d\n',Niterations)
fprintf('Number of available species: %d\n',Nspecies)
fprintf('Available values of theta: ')
theta
fprintf('Using theta = %g for plot.\n',theta(itheta))
fprintf('Available values of zeta: ')
zeta
fprintf('Using zeta = %g for plot.\n',zeta(izeta))

figure(1)
clf
set(gcf,'Color','w')

vx = linspace(-xMax,xMax,N);
vy = vx;
vz = vx;
[vx3D, vy3D, vz3D] = meshgrid(vx,vy,vz);
v3D = sqrt(vx3D.^2 + vy3D.^2 + vz3D.^2);
xi3D = vz3D ./ v3D;
data3D = zeros(size(vx3D));
fprintf('Processing L = ')
for L = 0:(Nxi-1)
    fprintf('%d ',L)
    if mod(L,10)==0 && L>0
        fprintf('\n')
    end
    if L==0
        LegendreP = ones(size(vx3D));
    elseif L==1
        LegendreP = xi3D;
        LegendreP_last = ones(size(vx3D));
    else
        % Recursion relation for Legendres
        LegendreP_lastlast = LegendreP;
        LegendreP = ((2*L-1)*xi3D .* LegendreP - (L-1)*LegendreP_last)/L;
        LegendreP_last = LegendreP_lastlast;
    end
    data3D = data3D + LegendreP .* interp1(x,squeeze(data(iteration,ispecies,itheta,izeta,L+1,:)), v3D, 'linear',0);
end
fprintf('\n')

maxval = max(max(max(abs(data3D))));

slice(vx3D,vy3D,vz3D,data3D, 0,0,0)
hold on
set(findobj(gca,'Type','Surface'),'EdgeColor','none')
%set(gca,'CLim',[-maxval,maxval])
daspect([1,1,1])
camproj('perspective')
camup([-0.9554    0.2711   -0.1168])
campos([ -6.5433  -19.4499    8.3788])
camva(19)
axis off
light('Position',[-5,-4,3])
set(findobj(gca,'type','surface'),...
    'AmbientStrength',.5,'DiffuseStrength',.8,...
    'SpecularStrength',.9,'SpecularExponent',25)

arrowMax = 4;
arrowTipMax = 4.5;
r=0.05;
arrow_r=0.1;
gray=0.5;
phi = linspace(0,2*pi,20);
w = [xMax,arrowMax];

% Add vz axis
[phi2D,w2D]=meshgrid(phi,w);
surf(r*cos(phi2D),r*sin(phi2D),w2D,'FaceColor',[gray,gray,gray],'EdgeColor','none','FaceLighting','gouraud')
% Add vz arrowhead
surf([arrow_r*cos(phi);zeros(size(phi))],[arrow_r*sin(phi);zeros(size(phi))],[arrowMax*ones(size(phi)); arrowTipMax*ones(size(phi))],...
    'FaceColor',[gray,gray,gray],'EdgeColor','none','FaceLighting','gouraud')
% Add vx axis
[phi2D,w2D]=meshgrid(phi,w);
surf(-w2D,r*cos(phi2D),r*sin(phi2D),'FaceColor',[gray,gray,gray],'EdgeColor','none','FaceLighting','gouraud')
% Add vx arrowhead
surf(-[arrowMax*ones(size(phi)); arrowTipMax*ones(size(phi))],[arrow_r*cos(phi);zeros(size(phi))],[arrow_r*sin(phi);zeros(size(phi))],...
    'FaceColor',[gray,gray,gray],'EdgeColor','none','FaceLighting','gouraud')
% Add vy axis
[phi2D,w2D]=meshgrid(phi,w);
surf(r*cos(phi2D),-w2D,r*sin(phi2D),'FaceColor',[gray,gray,gray],'EdgeColor','none','FaceLighting','gouraud')
% Add vy arrowhead
surf([arrow_r*cos(phi);zeros(size(phi))],-[arrowMax*ones(size(phi)); arrowTipMax*ones(size(phi))],[arrow_r*sin(phi);zeros(size(phi))],...
    'FaceColor',[gray,gray,gray],'EdgeColor','none','FaceLighting','gouraud')
