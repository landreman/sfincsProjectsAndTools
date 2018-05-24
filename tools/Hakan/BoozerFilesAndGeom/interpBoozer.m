function out=interpBoozer(Geom,invarname,invar,interptype)
%
% This function interpolates the Boozer data struct Geom to a new
% set of radii. The chosen new values invar of the flux label 
% given by the string invarname ('s', 'r', 'rnorm') are interpolated.
% The interptype can be 'linear' or 'nearest', default is 'linear'.

% invarname can also be 'index', in which case no interpolation takes place,
% the output is just Geom at the indices given in invar.

if nargin<4
  interptype='linear';
end

if strcmp(invarname,'s')
  s=invar;
elseif strcmp(invarname,'rnorm')
  s=interp1(Geom.rnorm,Geom.s,invar,'pchip');
  if any(isnan(s))
    error('Chosen rnorm out of range!')
  end
elseif strcmp(invarname,'r')
  s=interp1(Geom.rnorm*Geom.minorradiusW7AS,Geom.s,invar,'pchip');
  if any(isnan(s))
    error('Chosen r out of range!')
  end  
elseif strcmp(invarname,'index')
  s=Geom.s(invar);
  interptype='nearest';
else
  error('not implemented')
end

%if Geom.StelSym==0
%  error('not implemented')
%end
if size(s,1)>size(s,2)
  s=s';
end

out.StelSym=Geom.StelSym;
out.newsigncorr=Geom.newsigncorr;
if isfield(Geom,'headertext')
  out.headertext=Geom.headertext;
end
if isfield(Geom,'m0b')
  out.m0b=Geom.m0b;
  out.n0b=Geom.n0b;
end
out.nsurf=length(s);
out.Nperiods=Geom.Nperiods;
if isfield(Geom,'psi_a')
  out.psi_a=Geom.psi_a;
end
if isfield(Geom,'torfluxtot')
  out.torfluxtot=Geom.torfluxtot;
end
out.minorradiusW7AS=Geom.minorradiusW7AS;
out.majorradiusLastbcR00=Geom.majorradiusLastbcR00;
if isfield(Geom,'minorradiusVMEC')
  out.minorradiusVMEC=Geom.minorradiusVMEC;
end
if isfield(Geom,'majorradiusVMEC')
  out.majorradiusVMEC=Geom.majorradiusVMEC;
end
if isfield(Geom,'VolumeVMEC')
  out.VolumeVMEC=Geom.VolumeVMEC;
end
out.Bfilter=Geom.Bfilter;


out.s=s;

if strcmp(interptype,'linear') %This is the default case
  if isfield(Geom,'headertext')
    out.headertext.maincomment=sprintf([out.headertext.maincomment,'\n',...
                        'CC The file has been linearly interpolated to a new set of radii']);
  end
  corrindedge=find((Geom.s(end)<s) & (s<1.01*Geom.s(end)));
  if not(isempty(corrindedge))
    s(corrindedge)=Geom.s(end);
  end
  corrindcentre=find((Geom.s(1)>s) & (s>0.99*Geom.s(1)));
  if not(isempty(corrindcentre))
    s(corrindcentre)=Geom.s(1);
  end

  out.rnorm=interp1(Geom.s,Geom.rnorm,s);
  out.iota=interp1(Geom.s,Geom.iota,s);
  out.Bphi=interp1(Geom.s,Geom.Bphi,s);
  out.Btheta=interp1(Geom.s,Geom.Btheta,s);
  if isfield(Geom,'dpds')
    out.dpds=interp1(Geom.s,Geom.dpds,s);
  end
  if isfield(Geom,'dVdsoverNper')
    out.dVdsoverNper=interp1(Geom.s,Geom.dVdsoverNper,s);
  end
  out.B00=interp1(Geom.s,Geom.B00,s);
  out.R00=interp1(Geom.s,Geom.R00,s);


  for surfind=1:out.nsurf
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b%5i/%5i',surfind,out.nsurf)
    if isempty(find(Geom.s==s(surfind)))
      Lind=find(diff(sign(Geom.s-s(surfind))));
      if isempty(Lind)
        error('s is out of range!')
      end
    else
      exactind=find(Geom.s==s(surfind));
      if exactind==length(Geom.s)
        Lind=length(Geom.s)-1;
      else
        Lind=exactind;
      end
    end
    Rind=Lind+1;

    if not(Geom.StelSym)
      mode=0;
      for par=0:1

        pindsL=find(Geom.parity{Lind}==par);
        pindsR=find(Geom.parity{Rind}==par);
        
        mL=Geom.m{Lind}(pindsL);
        nL=Geom.n{Lind}(pindsL);
        mR=Geom.m{Rind}(pindsR);
        nR=Geom.n{Rind}(pindsR);
        
        allns=sort([nL,nR]);
        ns=allns([1,find(diff(allns))+1]);
        Nn=length(ns);
        
        for nind=1:Nn
          n=ns(nind);
          allms=sort([mL(find(nL==n)),mR(find(nR==n))]);
          ms=allms([1,find(diff(allms))+1]);
          Nm=length(ms);
          for mind=1:Nm       
            mode=mode+1;
            m=ms(mind);
            Lmnind=find((Geom.m{Lind}==m)&(Geom.n{Lind}==n)&(Geom.parity{Lind}==par));
            Lexists=not(isempty(Lmnind));
            Rmnind=find((Geom.m{Rind}==m)&(Geom.n{Rind}==n)&(Geom.parity{Rind}==par));
            Rexists=not(isempty(Rmnind));
            out.m{surfind}(mode)=m;
            out.n{surfind}(mode)=n;
            out.parity{surfind}(mode)=par;
            if Lexists && Rexists
              out.Bmn{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                             [Geom.Bmn{Lind}(Lmnind),Geom.Bmn{Rind}(Rmnind)],...
                                             s(surfind));
              out.dBmnds{surfind}(mode)=...
                  (Geom.Bmn{Rind}(Rmnind)-Geom.Bmn{Lind}(Lmnind))/(Geom.s(Rind)-Geom.s(Lind));
              out.R{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                           [Geom.R{Lind}(Lmnind),Geom.R{Rind}(Rmnind)],...
                                           s(surfind));
              out.Z{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                           [Geom.Z{Lind}(Lmnind),Geom.Z{Rind}(Rmnind)],...
                                           s(surfind));
              out.Dphi{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                              [Geom.Dphi{Lind}(Lmnind),Geom.Dphi{Rind}(Rmnind)],...
                                              s(surfind));
            elseif Lexists && not(Rexists)
              out.Bmn{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                             [Geom.Bmn{Lind}(Lmnind),0],...
                                             s(surfind));
              out.dBmnds{surfind}(mode)=...
                  (0-Geom.Bmn{Lind}(Lmnind))/(Geom.s(Rind)-Geom.s(Lind));
              out.R{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                           [Geom.R{Lind}(Lmnind),0],...
                                           s(surfind));
              out.Z{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                           [Geom.Z{Lind}(Lmnind),0],...
                                           s(surfind));
              out.Dphi{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                              [Geom.Dphi{Lind}(Lmnind),0],...
                                              s(surfind));        
            elseif not(Lexists) && Rexists
              out.Bmn{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                             [0,Geom.Bmn{Rind}(Rmnind)],...
                                             s(surfind));
              out.dBmnds{surfind}(mode)=...
                  (Geom.Bmn{Rind}(Rmnind)-0)/(Geom.s(Rind)-Geom.s(Lind));
              out.R{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                           [0,Geom.R{Rind}(Rmnind)],...
                                           s(surfind));
              out.Z{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                           [0,Geom.Z{Rind}(Rmnind)],...
                                           s(surfind));
              out.Dphi{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                              [0,Geom.Dphi{Rind}(Rmnind)],...
                                              s(surfind));        
            else
              error('this is impossible')
            end
          end    
        end
      end
    else %Stellarator symmetric case
      mL=Geom.m{Lind};
      nL=Geom.n{Lind};
      mR=Geom.m{Rind};
      nR=Geom.n{Rind};

      allns=sort([nL,nR]);
      ns=allns([1,find(diff(allns))+1]);
      Nn=length(ns);
      
      mode=0;
      for nind=1:Nn
        n=ns(nind);
        allms=sort([mL(find(nL==n)),mR(find(nR==n))]);
        ms=allms([1,find(diff(allms))+1]);
        Nm=length(ms);
        for mind=1:Nm
          mode=mode+1;
          m=ms(mind);
          Lmnind=find((mL==m)&(nL==n));
          Lexists=not(isempty(Lmnind));
          Rmnind=find((mR==m)&(nR==n));
          Rexists=not(isempty(Rmnind));
          out.m{surfind}(mode)=m;
          out.n{surfind}(mode)=n;
          if Lexists && Rexists
            out.Bmn{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                           [Geom.Bmn{Lind}(Lmnind),Geom.Bmn{Rind}(Rmnind)],...
                                           s(surfind));
            out.dBmnds{surfind}(mode)=...
                  (Geom.Bmn{Rind}(Rmnind)-Geom.Bmn{Lind}(Lmnind))/(Geom.s(Rind)-Geom.s(Lind));
            out.R{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                         [Geom.R{Lind}(Lmnind),Geom.R{Rind}(Rmnind)],...
                                         s(surfind));
            out.Z{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                         [Geom.Z{Lind}(Lmnind),Geom.Z{Rind}(Rmnind)],...
                                         s(surfind));
            if isfield(Geom,'Dphi')
              out.Dphi{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                              [Geom.Dphi{Lind}(Lmnind),Geom.Dphi{Rind}(Rmnind)],...
                                              s(surfind));
            end
          elseif Lexists && not(Rexists)
            out.Bmn{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                           [Geom.Bmn{Lind}(Lmnind),0],...
                                           s(surfind));
            out.dBmnds{surfind}(mode)=...
                  (0-Geom.Bmn{Lind}(Lmnind))/(Geom.s(Rind)-Geom.s(Lind));
            out.R{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                         [Geom.R{Lind}(Lmnind),0],...
                                         s(surfind));
            out.Z{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                         [Geom.Z{Lind}(Lmnind),0],...
                                         s(surfind));
            if isfield(Geom,'Dphi')
              out.Dphi{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                              [Geom.Dphi{Lind}(Lmnind),0],...
                                              s(surfind));    
            end
          elseif not(Lexists) && Rexists
            out.Bmn{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                           [0,Geom.Bmn{Rind}(Rmnind)],...
                                           s(surfind));
            out.dBmnds{surfind}(mode)=...
                  (Geom.Bmn{Rind}(Rmnind)-0)/(Geom.s(Rind)-Geom.s(Lind));
            out.R{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                         [0,Geom.R{Rind}(Rmnind)],...
                                         s(surfind));
            out.Z{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                         [0,Geom.Z{Rind}(Rmnind)],...
                                         s(surfind));
            if isfield(Geom,'Dphi')
              out.Dphi{surfind}(mode)=interp1(Geom.s(Lind:Rind)',...
                                              [0,Geom.Dphi{Rind}(Rmnind)],...
                                              s(surfind));   
            end
          else
            error('this is impossible')
          end
        end %for mind
      end %for nind
    end %if not(StelSym)
    out.Bnorm{surfind}=out.Bmn{surfind}/out.B00(surfind);
    out.nmodes(surfind)=mode;
  end %surfind loop
  fprintf(1,'\n')

  
elseif strcmp(interptype,'nearest')
  if isfield(Geom,'headertext')
    out.headertext.maincomment=sprintf([out.headertext.maincomment,'\n',...
                        'CC This file has only a subset of the flux surfaces of the original file.']);
  end
  if strcmp(invarname,'index')
    inds=invar;
    if size(inds,1)>size(inds,2)
      inds=inds';
    end
  else
    inds=interp1(Geom.s,1:length(Geom.s),max(min(s,Geom.s(end)),Geom.s(1)),'nearest');
  end
  
  out.rnorm=Geom.rnorm(inds);
  out.s=Geom.s(inds);
  out.iota=Geom.iota(inds);
  out.Bphi=Geom.Bphi(inds);
  out.Btheta=Geom.Btheta(inds);
  if isfield(Geom,'dpds')
    out.dpds=Geom.dpds(inds);
  end
  if isfield(Geom,'dVdsoverNper')
    out.dVdsoverNper=Geom.dVdsoverNper(inds);
  end
  out.B00=Geom.B00(inds);
  out.R00=Geom.R00(inds);
  out.nmodes=Geom.nmodes(inds);
  for surfind=1:out.nsurf
    ind=inds(surfind);
    out.Bmn{surfind}  =Geom.Bmn{ind};
    out.R{surfind}    =Geom.R{ind};
    out.Z{surfind}    =Geom.Z{ind};
    if isfield(Geom,'Dphi')
      out.Dphi{surfind} =Geom.Dphi{ind};
    end
    out.Bnorm{surfind}=Geom.Bnorm{ind};
    out.m{surfind}    =Geom.m{ind};
    out.n{surfind}    =Geom.n{ind};
    if not(Geom.StelSym)
      out.parity{surfind}=Geom.parity{ind};
    end
  end
else
  error('Unknown interpolation type!')
end