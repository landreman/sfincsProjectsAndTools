function writeBoozerfile(Geom,filename)

if not(filename(end-2:end)=='.bc')
  filename=[filename,'.bc'];
end

fid = fopen(filename,'wt');
if fid<0
 error(['Unable to open file: ',filename])
end

fprintf(fid,'%s\n',Geom.headertext.maincomment);
fprintf(fid,'CC The file has beeen interpolated to a new set of radii\n');

fprintf(fid,'%s\n',Geom.headertext.globalvars);
fprintf(fid,'%3d%5d%6d%4d%16.6E%10.5f%10.5f\n',...
        Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
        Geom.minorradius,Geom.majorradius);

for rind=1:Geom.nsurf
  fprintf(fid,'%s\n',Geom.headertext.surfvars);
  fprintf(fid,'%s\n',Geom.headertext.surfvarunits);
  fprintf(fid,'%12.4E%12.4E%12.4E%12.4E%12.4E%12.4E\n',...
          Geom.s(rind),Geom.iota(rind),...
          Geom.Bphi(rind)/(Geom.Nperiods/2/pi*(4*pi*1e-7)),...
          Geom.Btheta(rind)/(1/2/pi*(4*pi*1e-7)),...
          Geom.dpds(rind),Geom.dVdsoverNper(rind));
  fprintf(fid,'%s\n',Geom.headertext.datavars);
  for ind=1:Geom.nmodes(rind)
    fprintf(fid,'%5d%5d%16.8E%16.8E%16.8E%16.8E\n',...
            Geom.m{rind}(ind),Geom.n{rind}(ind),...
            Geom.R{rind}(ind),Geom.Z{rind}(ind),...
            Geom.Dphi{rind}(ind),Geom.Bmn{rind}(ind));
  end
end
fclose(fid);