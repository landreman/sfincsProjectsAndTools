addpath ~/Forskning/Stellarator/sfincs/gitsfincs/equilibria
addpath ~/Forskning/Stellarator/sfincs/sfincsProjectsAndTools/tools/Hakan/BoozerFilesAndGeom

[P.rnorm,P.r,P.ne20,P.ni20,P.TekeV,P.TikeV,P.ErkVm,P.Zeff,P.flux21,P.QMW,comments]...
    =loadprofiles();


boozerfile='w7x-lim-op1_1.bc';
%only run once
%Geom=readBoozerfile(boozerfile);
G=Geom.Bphi;
I=Geom.Btheta;


drnorm=[0,diff(rnorm)]+[diff(rnorm),0];

P.dne20drnorm=([0,diff(P.ne20)]+[diff(P.ne20),0])./drnorm;
P.dni20drnorm=([0,diff(P.ni20)]+[diff(P.ni20),0])./drnorm;
P.dTekeVdrnorm=([0,diff(P.TekeV)]+[diff(P.TekeV),0])./drnorm;
P.dTikeVdrnorm=([0,diff(P.TikeV)]+[diff(P.TikeV),0])./drnorm;

eps0=8.8542e-12;
mp=1.6726e-27;
me=9.1094e-31;
e=1.6022e-19;
mi=1*mp; %HS assumption


P.lnLambda=25.3-1.15*log10(P.ni20*1e14)+2.30*log10(P.TikeV*1e3);
P.nuii=1/(4*pi*eps0)^2*4*sqrt(2*pi)/3*...
     e^4*P.Zeff.^4.*(P.ni20*1e20)/sqrt(mi)...
     ./(P.TikeV*1e3*e).^(3/2).*P.lnLambda;

S.rnorm=interp1(Geom.rnorm,Geom.rnorm,[0.15:0.05:0.95],'nearest');
%S.rnorm=Geom.rnorm;

S.G=interp1(Geom.rnorm,Geom.Bphi,S.rnorm);
S.I=interp1(Geom.rnorm,Geom.Btheta,S.rnorm);
S.iota=interp1(Geom.rnorm,Geom.iota,S.rnorm);
S.B00=interp1(Geom.rnorm,Geom.B00,S.rnorm);

S.ne20=interp1(P.rnorm,P.ne20,S.rnorm);
S.ni20=interp1(P.rnorm,P.ni20,S.rnorm);
S.TekeV=interp1(P.rnorm,P.TekeV,S.rnorm);
S.TikeV=interp1(P.rnorm,P.TikeV,S.rnorm);
S.dne20drnorm=interp1(P.rnorm,P.dne20drnorm,S.rnorm);
S.dni20drnorm=interp1(P.rnorm,P.dni20drnorm,S.rnorm);
S.dTekeVdrnorm=interp1(P.rnorm,P.dTekeVdrnorm,S.rnorm);
S.dTikeVdrnorm=interp1(P.rnorm,P.dTikeVdrnorm,S.rnorm);
S.dne20dPsiN=S.dne20drnorm/2./S.rnorm;
S.dni20dPsiN=S.dni20drnorm/2./S.rnorm;
S.dTekeVdPsiN=S.dTekeVdrnorm/2./S.rnorm;
S.dTikeVdPsiN=S.dTikeVdrnorm/2./S.rnorm;

S.ErkVm=interp1(P.rnorm,P.ErkVm,S.rnorm);
S.dPotentialkVdPsiN = -S.ErkVm*Geom.minorradiusW7AS/2./S.rnorm;

S.vTi=sqrt(S.TikeV*1e3*e/mi*2);
S.nuii=interp1(P.rnorm,P.nuii,S.rnorm);
S.nuPrime=S.nuii.*(S.G+S.iota.*S.I)./S.vTi./S.B00;

fig(5)
plot(S.rnorm,S.nuPrime)
title('\nu prime')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Make runspecfile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if length(S.rnorm)==1
%  runspecfilename=['w7x_20160309.010_0.32_r',num2str(round(S.rnorm*1000)/1000),'.runspec'];
%else
runspecfilename='w7x_20160309.010_0.32.runspec';
%end

fid=fopen(runspecfilename,'w');
fprintf(fid,'! ---------------- Hennings comments: -----------------\n');
for ind=1:length(comments)
  fprintf(fid,['! ',comments{ind},'\n']);
end
fprintf(fid,'! ---------------- Boozer file comments: ----------------\n');
bccom=strread(Geom.headertext.maincomment, '%s', 'delimiter', sprintf('\n'));
for ind=1:length(bccom)
  fprintf(fid,['! ',bccom{ind},'\n']);
end
fprintf(fid,'!-----------------------------------------------\n');
fprintf(fid,['!   rN_wish',...
                 '     nHats_1',...
                 '     nHats_2',...
                 '     THats_1',...
                 '     THats_2',...
                 ' dNHatdpsiNs_1',...
                 ' dNHatdpsiNs_2',...
                 ' dTHatdpsiNs_1',...
                 ' dTHatdpsiNs_2',...
                 '  dPhiHatdpsiN',...
                 '  EParallelHat\n']);

for ind=1:length(S.rnorm)
  fprintf(fid,'%11.5e %11.5e %11.5e %11.5e %11.5e % 12.6e % 12.6e % 12.6e % 12.6e % 12.6e % 12.6e\n',...
          S.rnorm(ind),S.ne20(ind),S.ni20(ind),S.TekeV(ind),S.TikeV(ind), ...
          S.dne20dPsiN(ind),S.dni20dPsiN(ind), ...
          S.dTekeVdPsiN(ind),S.dTikeVdPsiN(ind), ...
          S.dPotentialkVdPsiN(ind),0);
end


fclose(fid);