function out=reduceBoozertoLCD_Nperiod(Geom)
%Divide all n in the Boozer data with the lowest common denominator.

if Geom.StelSym==0
  error('not implemented yet')
end

out=Geom;

ns=[];

for rind=1:Geom.nsurf
  ns=[ns,Geom.n{rind}];
end

ns=abs(ns(find(ns~=0)));

found=0;
Ntry=max(ns)+1;
while Ntry>1 && found==0
  Ntry=Ntry-1;
  found=all(mod(ns,Ntry)==0);
end
Nper=Ntry;

out.Nperiods=Geom.Nperiods*Nper;
for rind=1:Geom.nsurf
  out.n{rind}=Geom.n{rind}/Nper;
  %out.Bphi(rind)=Geom.Bphi(rind)/Nper;
  %out.iota(rind)=Geom.iota(rind)/Nper;
  
  out.dVdsoverNper(rind)=Geom.dVdsoverNper(rind)/Nper;
  out.Dphi{rind}=out.Dphi{rind}*Nper;
end