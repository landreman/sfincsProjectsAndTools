% This example file shows how to use CoordPlots to make 3D
% plots of the flux surfaces

%input
Boozerfile='~/Forskning/Stellarator/sfincs/gitsfincs/equilibria/w7x-sc1.bc';
rnormwish=0.9;
Ntheta=100;
Nzeta=100;

%Do this only once, if it takes a lot of time
Geom=readBoozerfile(Boozerfile);

rind=findnearest(Geom.rnorm,rnormwish);

[Ham,Booz,Cyl]=makeHamada(Geom,rind,Ntheta,Nzeta);

%CoordPlots(Ham)
CoordPlots(Booz)
%CoordPlots(Cyl)

if 0
  fig(1)
  %surf(Booz.zeta,Booz.theta,Booz.cylth);view(0,90);shading flat;colorbar
  %xlabel('\zeta');ylabel('\theta')
  %surf(Cyl.zeta,Cyl.theta,Cyl.cylth);view(0,90);shading flat;colorbar
  fig(1);surf(Cyl.zeta,Cyl.theta,Cyl.Dthetacylth);view(0,90);shading flat;colorbar
  fig(2);surf(Booz.zeta,Booz.theta,Booz.Dthetacylth);view(0,90);shading flat;colorbar

  fig(3);surf(Cyl.zeta,Cyl.theta,Cyl.Dzetacylphi);view(0,90);shading flat;colorbar
  fig(4);surf(Booz.zeta,Booz.theta,Booz.Dzetacylphi);view(0,90);shading flat;colorbar


  fig(2)
  surf(Booz.zeta,Booz.theta,Booz.cylr);view(0,90);shading flat;colorbar
  xlabel('\zeta');ylabel('\theta')

  fig(15)
  plot(Booz.cylR,Booz.cylZ)
end