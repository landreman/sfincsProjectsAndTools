

%add the paths to
addpath([getenv('SFINCS_PROJECTSANDTOOLS_HOME'),'/tools/Hakan/BoozerFilesAndGeom']) %This is where readBoozerfile.m resides.
addpath([getenv('SFINCS_PROJECTSANDTOOLS_HOME'),'/tools/Hakan/version3Scans']) 
addpath([getenv('SFINCS_HOME'),'/equilibria/']) 

Geom=readBoozerfile('w7x-sc1.bc');

rwishs=0.1:0.1:1;
for ii=1:length(rwishs)
  rinds(ii)=findnearest(Geom.rnorm,rwishs(ii));
end

rnorms=Geom.rnorm(rinds);
ss=Geom.s(rinds);

n0=0.7; %($10^{20} m^{-3}$),
n1=0.7;
n2=2.5;
n3=3.0;

T0=1.0; %(keV)
T1=0.7;
T2=1.5;
T3=3.0;


n=n0*(1-n1*ss.^n2).^n3;
T=T0*(1-T1*ss.^T2).^T3;
dnds=n0*n3.*(1-n1*ss.^n2).^(n3-1).*(-n1*n2*ss.^(n2-1));
dTds=T0*T3.*(1-T1*ss.^T2).^(T3-1).*(-T1*T2*ss.^(T2-1));



outputonscreenonly=0;
 
runspecfile='runspec.dat';
if outputonscreenonly
  runspec_fid=1;
else
  runspec_fid=fopen(runspecfile,'w');
end


%HEADER
fprintf(runspec_fid,'! This file was generated by make_SatakeReganaBenchmarkrunspecfile.m\n');
fprintf(runspec_fid,'!\n');
fprintf(runspec_fid,['!      rN_wish',...
                    '       nHats_1',...
                    '       nHats_2',...
                    '       THats_1',...
                    '       THats_2',...
                    ' dNHatdpsiNs_1',...
                    ' dNHatdpsiNs_2',...
                    ' dTHatdpsiNs_1',...
                    ' dTHatdpsiNs_2\n']);
          
for rind=1:length(rnorms)
  rnorm=rnorms(rind);

  fprintf(runspec_fid,...
        ['%14.6e%14.6e%14.6e%14.6e%14.6e',...
         '%14.6e%14.6e%14.6e%14.6e\n'],...
        rnorm,n(rind),n(rind),T(rind),T(rind),...
         dnds(rind),dnds(rind),dTds(rind),dTds(rind) );
end

if not(outputonscreenonly)
  fclose(runspec_fid);
end
