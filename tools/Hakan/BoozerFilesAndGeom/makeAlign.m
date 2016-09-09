function [Align,Discr]=makeAlign(Discr)
%
% This creates a discretisation in a filed aligned coordinate system.
%
% Discr could be Booz or Ham or some other discretisation with straight field lines.
% Discr in the output must be the same as the input.
%
% The coordinates are called apol,ator
% where 
% apol (from 0 to 2*pi) equals the input poloidal coordinate and 
% ator (from 0 to 2*pi) is the point where the magnetic field line crosses the outboard midplane apol=0
%
% Note that the Align coordinates are not periodic in the poloidal angle (which
% measures distance along the magnetic field line)

if Discr.Nperiods ~=1
  error(['Aligned coordinates only works for N=1. Otherwise they would depend on which ', ...
         'toroidal period one is in!'])
end

if strcmp(Discr.name,'Boozer')
  bpol=Discr.theta;
  btor=Discr.zeta;
  Align.name='Aligned Boozer';
elseif strcmp(Discr.name,'Hamada')
  bpol=Discr.vthet;
  btor=Discr.phi;
  Align.name='Aligned Hamada';
elseif strcmp(Discr.name,'Pest')
  bpol=Discr.ptheta;
  btor=Discr.pzeta;
  Align.name='Aligned Pest';
else
  error('Not a valid input coordinate discretisation to makeAlign.m!')
end

Discr_apol=bpol;
Discr_ator=mod(btor-bpol/Discr.iota,2*pi/Discr.Nperiods);

Align.Nperiods=Discr.Nperiods;
Align.iota=Discr.iota;

Align_apol=bpol; %Just copy the uniform discretisation matix from 0 to 2*pi
Align_ator=btor; %Just copy the uniform discretisation matix from 0 to 2*pi

Align_bpol=Align_apol;
Align_btor=Align_ator+Align_apol/Discr.iota;


Align.B=interp2_cyclic(bpol,btor,Discr.B,Align_bpol,Align_btor,Discr.Nperiods);


if strcmp(Discr.name,'Boozer')
  Align.theta=Align_bpol;
  Align.zeta=Align_btor;
  Align.altheta=Align_apol;
  Align.alzeta=Align_ator;
  Discr.altheta=Discr_apol;
  Discr.alzeta=Discr_ator;
elseif strcmp(Discr.name,'Hamada')
  Align.vthet=Align_bpol;
  Align.phi=Align_btor;
  Align.alvthet=Align_apol;
  Align.alphi=Align_ator;
  Discr.alvthet=Discr_apol;
  Discr.alphi=Discr_ator;
elseif strcmp(Discr.name,'Pest')
  Align.pthet=Align_bpol;
  Align.pzeta=Align_btor;
  Align.alptheta=Align_apol;
  Align.alpzeta=Align_ator;
  Discr.alptheta=Discr_apol;
  Discr.alpzeta=Discr_ator;
end

%fig(1)
%surf(Align.apol,Align.ator,Align.btor)
%surf(Discr.theta,Discr.zeta,Discr_ator)
%view(0,90)