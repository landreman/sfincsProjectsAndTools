function writeBoozerfile(Geom,filename)

if not(filename(end-2:end)=='.bc')
  filename=[filename,'.bc'];
end

fid = fopen(filename,'wt');
if fid<0
 error(['Unable to open file: ',filename])
end

% Change the signs from right to left handed (r,poloidal,toroidal)
signchange=-1; %sign changer
Geom.psi_a=Geom.psi_a*signchange;
Geom.torfluxtot=Geom.torfluxtot*signchange;
Geom.Bphi=Geom.Bphi*signchange^2;
Geom.Btheta=Geom.Btheta*signchange;
Geom.iota=Geom.iota*signchange;
Geom.dVdsoverNper=Geom.dVdsoverNper*signchange;
for tmpind=1:length(Geom.n)
  Geom.n{tmpind}=-Geom.n{tmpind};
  Geom.Dphi{tmpind}=-Geom.Dphi{tmpind};
end


%Begin writing to the file

fprintf(fid,'%s\n',Geom.headertext.maincomment);

if isfield(Geom.headertext,'globalvars')
  fprintf(fid,'%s\n',Geom.headertext.globalvars);
else
  if Geom.StelSym %JG includes a little more info in this case
    fprintf(fid,[' m0b  n0b nsurf nper  flux/[Tm^2]     a/[m]     R/[m]   avol/[m]  ' ...
                 'Rvol/[m] Vol/[m^3]\n']);
  else %and ES a little less
    fprintf(fid,' m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]\n');
  end
end

minorradiusVMEC_exists=0;
if isfield(Geom,'minorradiusVMEC')
  if not(isnan(Geom.minorradiusVMEC))
    minorradiusVMEC_exists=1;
  end
end
majorradiusVMEC_exists=0;
if isfield(Geom,'majorradiusVMEC')
  if not(isnan(Geom.majorradiusVMEC))
    majorradiusVMEC_exists=1;
  end
end
VolumeVMEC_exists=0;
if isfield(Geom,'VolumeVMEC')
  if not(isnan(Geom.VolumeVMEC))
    VolumeVMEC_exists=1;
  end
end


if isnan(Geom.majorradiusLastbcR00)
  Geom.majorradiusLastbcR00=0;
end


if Geom.StelSym
  if minorradiusVMEC_exists && not(majorradiusVMEC_exists) && not(VolumeVMEC_exists)
    fprintf(fid,'%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f\n',...
            Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
            Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,Geom.minorradiusVMEC);
  elseif minorradiusVMEC_exists && majorradiusVMEC_exists && not(VolumeVMEC_exists)
    fprintf(fid,'%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n',...
            Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
            Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,...
            Geom.minorradiusVMEC,Geom.majorradiusVMEC,...
            pi*Geom.minorradiusVMEC^2*2*pi*Geom.majorradiusVMEC);    
  elseif minorradiusVMEC_exists && majorradiusVMEC_exists && VolumeVMEC_exists
    fprintf(fid,'%3d%5d%6d%4d%16.6E%10.5f%10.5f%11.5f%11.5f%11.5f\n',...
            Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
            Geom.minorradiusW7AS,Geom.majorradiusLastbcR00,...
            Geom.minorradiusVMEC,Geom.majorradiusVMEC,Geom.VolumeVMEC);        
  else
    fprintf(fid,'%3d%5d%6d%4d%16.6E%10.5f%10.5f\n',...
            Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
            Geom.minorradiusW7AS,Geom.majorradiusLastbcR00);
  end
else
  %Difficult! I Think Erika uses VMEC quantities here but I am not sure!
  if minorradiusVMEC_exists && majorradiusVMEC_exists
    fprintf(fid,'%4d%6d%7d%6d%15.8f%15.8f%15.8f\n',...
            Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
            Geom.minorradiusVMEC,Geom.majorradiusVMEC);
  else
    error('Did not find the minorradiusVMEC and majorradiusVMEC!')
    %fprintf(fid,'%4d%6d%7d%6d%15.8f%15.8f%15.8f\n',...
    %        Geom.m0b,Geom.n0b,Geom.nsurf,Geom.Nperiods,Geom.torfluxtot,...
    %        Geom.minorradiusW7AS,Geom.majorradiusLastbcR00);  
  end
end
for rind=1:Geom.nsurf
  if isfield(Geom.headertext,'surfvars')
    fprintf(fid,'%s\n',Geom.headertext.surfvars);
  else
    if Geom.StelSym
      fprintf(fid,['       s         iota  curr_pol/nper    curr_tor    pprime   sqrt ' ...
                   'g(0,0)\n']);
    else
      fprintf(fid,['        s               iota           Jpol/nper          Itor',...
                   '            pprime         sqrt g(0,0)\n']);
    end
  end
  if isfield(Geom.headertext,'surfvarunits')
    fprintf(fid,'%s\n',Geom.headertext.surfvarunits);
  else
    if Geom.StelSym
      fprintf(fid,['                            [A]            [A]   dp/ds,[Pa] ' ...
                   '(dV/ds)/nper\n']);
    else
      fprintf(fid,['                                          [A]           [A] ',...
                   '            [Pa]         (dV/ds)/nper\n']);
    end
  end
  if Geom.StelSym
    fprintf(fid,'%12.4E%12.4E%12.4E%12.4E%12.4E%12.4E\n',...
          Geom.s(rind),Geom.iota(rind),...
          Geom.Bphi(rind)/(Geom.Nperiods/2/pi*(4*pi*1e-7)),...
          Geom.Btheta(rind)/(1/2/pi*(4*pi*1e-7)),...
          Geom.dpds(rind),Geom.dVdsoverNper(rind));
  else
    fprintf(fid,'%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E\n',...
          Geom.s(rind),Geom.iota(rind),...
          Geom.Bphi(rind)/(Geom.Nperiods/2/pi*(4*pi*1e-7)),...
          Geom.Btheta(rind)/(1/2/pi*(4*pi*1e-7)),...
          Geom.dpds(rind),Geom.dVdsoverNper(rind));
  end
  if isfield(Geom.headertext,'datavars')
    fprintf(fid,'%s\n',Geom.headertext.datavars);
  else
    if Geom.StelSym
      fprintf(fid,['    m    n        r/[m]           z/[m] (phib-phi)*nper/twopi     ',...
                   'bmn/[T]\n']);
    else
      fprintf(fid,['    m    n      rmnc [m]         rmns [m]         zmnc [m]         ',...
                   'zmns [m]         vmnc [ ]         vmns [ ]         ',...
                   'bmnc [T]         bmns [T]\n']);
    end
  end
  if Geom.StelSym
    for ind=1:Geom.nmodes(rind)
      fprintf(fid,'%5d%5d%16.8E%16.8E%16.8E%16.8E\n',...
              Geom.m{rind}(ind),Geom.n{rind}(ind),...
              Geom.R{rind}(ind),Geom.Z{rind}(ind),...
              Geom.Dphi{rind}(ind),Geom.Bmn{rind}(ind));
    end
  else     
    allms=sort(Geom.m{rind});
    ms=allms([1,find(diff(allms))+1]);
    Nm=length(ms);
    for mind=1:Nm
      m=ms(mind);
      allns=sort(Geom.n{rind}(find(Geom.m{rind}==m)));
      ns=allns([1,find(diff(allns))+1]);
      Nn=length(ns);
      for nind=1:Nn
        n=ns(nind);
        mninds=find((Geom.n{rind}==n)&(Geom.m{rind}==m));
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