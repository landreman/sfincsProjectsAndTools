function out=interpBoozer(Geom,invarname,invar)

if strcmp(invarname,'s')
  s=invar;
elseif strcmp(invarname,'rnorm')
  s=interp1(Geom.rnorm,Geom.s,invar);
else
  error('not implemented')
end

%if Geom.StelSym==0
%  error('not implemented')
%end

out.StelSym=Geom.StelSym;
out.headertext=Geom.headertext;
out.m0b=Geom.m0b;
out.n0b=Geom.n0b;
out.nsurf=length(s);
out.Nperiods=Geom.Nperiods;
out.torfluxtot=Geom.torfluxtot;
out.minorradius=Geom.minorradius;
out.majorradius=Geom.majorradius;
out.Bfilter=Geom.Bfilter;


out.rnorm=interp1(Geom.s,Geom.rnorm,s);
out.s=s;
out.iota=interp1(Geom.s,Geom.iota,s);
out.Bphi=interp1(Geom.s,Geom.Bphi,s);
out.Btheta=interp1(Geom.s,Geom.Btheta,s);
out.dpds=interp1(Geom.s,Geom.dpds,s);
out.dVdsoverNper=interp1(Geom.s,Geom.dVdsoverNper,s);
out.B00=interp1(Geom.s,Geom.B00,s);
out.R00=interp1(Geom.s,Geom.R00,s);

for surfind=1:out.nsurf

  Lind=find(diff(sign(Geom.s-s(surfind))));
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
  out.Bnorm{surfind}=out.Bmn{surfind}/out.B00(surfind);
  out.nmodes(surfind)=mode;
end
