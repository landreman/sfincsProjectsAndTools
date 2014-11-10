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
  if Geom.StelSym
    for ind=1:Geom.nmodes(rind)
      fprintf(fid,'%5d%5d%16.8E%16.8E%16.8E%16.8E\n',...
              Geom.m{rind}(ind),Geom.n{rind}(ind),...
              Geom.R{rind}(ind),Geom.Z{rind}(ind),...
              Geom.Dphi{rind}(ind),Geom.Bmn{rind}(ind));
    end
  else     
    allns=sort(Geom.n{rind});
    ns=allns([1,find(diff(allns))+1]);
    Nn=length(ns);
    for nind=1:Nn
      n=ns(nind);
      allms=sort(Geom.m{rind}(find(Geom.n{rind}==n)));
      ms=allms([1,find(diff(allms))+1]);
      Nm=length(ms);
      for mind=1:Nm
        m=ms(mind);
        mninds=find((Geom.n{rind}==n)&(Geom.m{rind}=m));
        if length(mninds)>2 || isempty(mninds)
          error('This is impossible!')
        end
        rmnc=0;rmns=0;zmnc=0;zmns=0;
        vmnc=0;vmns=0;bmnc=0;bmns=0;
        for iind=1:length(mninds)
          if Geom.parity{rind}(mninds(iind))==1 %1<=>cosinus Bmn (stelsym part)
            rmnc=Geom.R{rind}(mninds(iind));
            zmns=Geom.Z{rind}(mninds(iind));
            vmns=Geom.Dphi{rind}(mninds(iind));
            bmnc=Geom.Bmn{rind}(mninds(iind));
          else %0<=>sinus Bmn (non-stelsym part)
            rmns=Geom.R{rind}(mninds(iind));
            zmnc=Geom.Z{rind}(mninds(iind));
            vmnc=Geom.Dphi{rind}(mninds(iind));
            bmns=Geom.Bmn{rind}(mninds(iind));            
          end
        end
      fprintf(fid,'%5d%5d%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E\n',...
              m,n,rmnc,rmns,zmnc,zmns,vmnc,vmns,bmnc,bmns);        
      end
    end
  end
end
fclose(fid);