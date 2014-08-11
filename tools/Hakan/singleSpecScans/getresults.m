function out=getresults(directory,sortafter)
% Retrieves the output from the sfincs runs in the numbered subdirectories in
% "directory". Only succesful runs are loaded. The output struct "out"
% is sorted after the vaiable with the name in the input "sortafter"

H=loadallh5(directory);

ind=0;
for hind=1:length(H)
  if isstruct(H{hind})
    ind=ind+1;
    goodhinds(ind)=hind;
    out.run(ind).dir       =[H{hind}.dirpath,H{hind}.rundir];
    out.run(ind).time      =runtime(out.run(ind).dir);
    out.RHSMode(ind)      =H{hind}.run1.RHSMode;
      
    out.normradius(ind)   =H{hind}.run1.normradius;
    out.nuN(ind)          =H{hind}.run1.nuN;
    out.nuPrime(ind)      =H{hind}.run1.nuPrime;
    out.EStar(ind)        =H{hind}.run1.EStar;
    out.dPhiHatdpsi(ind)  =H{hind}.run1.dPhiHatdpsi;
    out.nHat(ind)         =H{hind}.run1.nHat;
    out.THat(ind)         =H{hind}.run1.THat;
    
    out.Ntheta(ind)            =H{hind}.run1.Ntheta;
    out.Nzeta(ind)             =H{hind}.run1.Nzeta;
    out.Nxi(ind)               =H{hind}.run1.Nxi;
    out.NL(ind)                =H{hind}.run1.NL;
    out.Nx(ind)                =H{hind}.run1.Nx;
    out.NxPotentialsPerVth(ind)=H{hind}.run1.NxPotentialsPerVth;
    out.xMax(ind)              =H{hind}.run1.xMax;
    out.solverTolerance(ind)   =H{hind}.run1.solverTolerance;
    
    if isfield(H{hind}.run1,'NTV') %Only for backward compability
      out.tauhat_s(ind)     =H{hind}.run1.NTV;
    else
      out.tauhat_s(ind)     =H{hind}.run1.NTVsingle;
      out.tauhat_m(ind)     =H{hind}.run1.NTVmulti;
    end  
    out.particleFlux(ind) =H{hind}.run1.particleFlux;
    out.heatFlux(ind)     =H{hind}.run1.heatFlux;
    out.momentumFlux(ind) =H{hind}.run1.momentumFlux;
    out.FSAFlow(ind)      =H{hind}.run1.FSAFlow;
  end
end
out.NumElements=ind;

if out.NumElements>0
  if all((out.RHSMode==2)|(out.RHSMode==0)); %The 0 is for a certain type of error
    out.transportMatrix=zeros(out.NumElements,3,3);
    for ind=1:out.NumElements
      out.transportMatrix(ind,:,:)=H{goodhinds(ind)}.run1.transportMatrix;
    end
  end
end

if out.NumElements>1
  if nargin>1 %sort the results
    [dummy,orderedInds]=sort(getfield(out,sortafter));
    fnames=fieldnames(out);
    for find=1:length(fnames)
      if strcmp(fnames{find},'transportMatrix')
        out.transportMatrix=out.transportMatrix(orderedInds,:,:);
      elseif not(strcmp(fnames{find},'NumElements'))
        tmp=getfield(out,fnames{find});
        out=setfield(out,fnames{find},tmp(orderedInds));
      end
    end
  end
end

