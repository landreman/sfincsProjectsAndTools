% This example file shows how to make 3D plots on a flux surface

%input
Boozerfile='/afs/ipp-garching.mpg.de/home/s/smithh/Forskning/Stellarator/sfincs/gitsfincs/equilibria/w7x-sc1.bc';
rnormwish=0.9;
Ntheta=41; %must be odd
Nzeta=41;  %must be odd
min_Bmn=0; %only load Bmn's bigger than this
max_m=inf;
maxabs_n=inf;

%Do this only once, if it takes a lot of time
Geom=readBoozerfile(Boozerfile,min_Bmn,max_m,maxabs_n,'StelSym','signcorr2');

rind=findnearest(Geom.rnorm,rnormwish);

[Ham,Booz,Cyl,Pest]=makeHamada(Geom,rind,Ntheta,Nzeta,'forceSize');

%Ham contains everything disretized on a uniform grid in Hamada coordinates
%Booz contains everything disretized on a uniform grid in Booz coordinates
%Pest contains everything disretized on a uniform grid in Pest coordinates
%type help makeHamada to get more info

%Try these out if you want all coordinates plotted in 3D using one of the discretizations 
%CoordPlots(Ham)
%CoordPlots(Booz);
%CoordPlots(Pest)
%CoordPlots(Cyl)


if 1
  fig(1)
  surf(Booz.zeta,Booz.theta,Booz.B);view(0,90);shading flat;colorbar
  xlabel('Boozer \zeta');ylabel('Boozer \theta')
  title('B in Boozer discretisation')
  
  fig(2)
  surf(Pest.zeta,Pest.theta,Pest.B);view(0,90);shading flat;colorbar
  xlabel('Boozer \zeta');ylabel('Boozer \theta')
  title('B in Pest discretisation')

  fig(3)
  surf(Booz.pzeta,Booz.ptheta,Booz.B);view(0,90);shading flat;colorbar
  xlabel('Pest \zeta');ylabel('Pest \theta')
  title('B in Boozer discretisation')
  
  fig(4)
  surf(Pest.pzeta,Pest.ptheta,Pest.B);view(0,90);shading flat;colorbar
  xlabel('Pest \zeta');ylabel('Pest \theta')
  title('B in Pest discretisation')
end

%Plot a quantiB on a 3D surface using the Pest discretisation:
fig(15)
subsampling_pol=1; %for plot efficiency. 1 means all points, 2 every second etc.
subsampling_tor=1;
plotonfluxsurf(Pest,Pest.B,subsampling_pol,subsampling_tor)
