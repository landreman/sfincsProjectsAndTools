%%%%%%%%%% Equilibrium data %%%%%%%%%%%%%%%%%%%%%
equilibriumfile='~/sfincs/sfincs/equilibria/w7x-sc1.bc';
addpath ~/sfincs/sfincsProjectsAndTools/tools/Hakan/BoozerFilesAndGeom
Geom=readBoozerfile(equilibriumfile);
rind=24;
%Geom.rnorm(rind)

G=Geom.Bphi(rind);
I=Geom.Btheta(rind);
iota=Geom.iota(rind);
B00=Geom.B00(rind);

dPsidr=Geom.rnorm(rind)*Geom.minorradius*B00;

Chandra1=(erf(sqrt(1))-sqrt(1)*2/sqrt(pi).*exp(-1))./2./1;
cmulovernuPrime=...
   B00/(G+iota*I)*3*sqrt(pi)/4*(erf(1)-Chandra1)

g11overL11s=sqrt(pi)/8*G/(G+iota*I)*G/B00/dPsidr^2;
g31overL31s=sqrt(pi)/4*G/dPsidr;
g33overL33s=-sqrt(pi)/2*(G+iota*I)*B00;


%%%%%%%%%%%%%%%%%% DKES data %%%%%%%%%%%%%%%%%%%%%%%5
dk=read_dkes_dkfile('~/sfincs/sfincsProjectsAndTools/tools/Hakan/dkes/w7x-sc1-ecb2.dk');
data=dk.data{4};


%Look at Er=0 values

Er0ind=find(data.EovervB==0);


cmul=data.cmul(Er0ind);
g11_i=data.g11_i(Er0ind)/B00^2;
g13_i=data.g13_i(Er0ind);
g33_i=data.g33_i(Er0ind)*B00^2;
g11_e=data.g11_e(Er0ind)/B00^2;
g13_e=data.g13_e(Er0ind);
g33_e=data.g33_e(Er0ind)*B00^2;

fig(1)
loglog(cmul,-g11_i,'k',cmul,-(g11_i+g11_e),'g',cmul,-(g11_i-g11_e),'g')
xlabel('cmul')
ylabel('g_{11}')

fig(2)
semilogx(cmul,g13_i,'k',cmul,g13_i+g13_e,'g',cmul,g13_i-g13_e,'g')
xlabel('cmul')
ylabel('g_{13}')

fig(3)
loglog(cmul,-g33_i,'k',cmul,-(g33_i+g33_e),'g',cmul,-(g33_i-g33_e),'g')
xlabel('cmul')
ylabel('g_{33}')


%%%%%%%%%%%%%%% Monoenergetic sfincs data

%nuPrime = 3e-5
nu3em5dir='~/sfincs/sfincs/fortran/version3/runs/monoenergetic/w7x-sc1-ecb2/r0p5Er0nu3e-5/41_220_220/try16_4_4';
[runs,missing]=getconvergencescanresults(nu3em5dir,0);

runind=runs.NumElements;%the last one is the base case
transportCoeffs(1,:,:)=squeeze(runs.transportCoeffs(runind,:,:));
nuPrime(1)=runs.nuPrime(runind);

%nuPrime = 1e-4
nu1em4dir='~/sfincs/sfincs/fortran/version3/runs/monoenergetic/w7x-sc1-ecb2/r0p5Er0nu1e-4/41_155_155/try16_4_4';
[runs,missing]=getconvergencescanresults(nu1em4dir,0);

runind=runs.NumElements;%the last one is the base case
transportCoeffs(2,:,:)=squeeze(runs.transportCoeffs(runind,:,:));
nuPrime(2)=runs.nuPrime(runind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
loglog(cmulovernuPrime*nuPrime,-g11overL11s*transportCoeffs(:,1,1),'r+')
hold off

figure(2)
hold on
semilogx(cmulovernuPrime*nuPrime,...
         g31overL31s*transportCoeffs(:,1,2),'r+',...
         cmulovernuPrime*nuPrime,...
         g31overL31s*transportCoeffs(:,2,1),'r+')
hold off

figure(3)
hold on
loglog(cmulovernuPrime*nuPrime,-g33overL33s*transportCoeffs(:,2,2),'r+')
hold off

