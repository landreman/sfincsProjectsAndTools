function [out,missing]=getresults(directory,sortafter)
% Retrieves the output from the sfincs runs with the same discretisation
% in the numbered subdirectories in "directory". 
% Only succesful runs are loaded. The output struct "out"
% is sorted after the vaiable with the name in the input "sortafter"
%directory
H=loadallh5(directory);

ind=0;  %Index for successful runs
mind=0; %index for missing runs
missing=[];
%for backward compatibility, check if x is an output
xoutputexists=0;
for hind=1:length(H)
  if isstruct(H{hind})
    if isfield(H{hind}.run1,'fNormIsotropicBeforeSurfaceIntegral')
      xoutputexists=1;
      %xoutputlength=length(H{hind}.run1.x);
    end
  end
end

for hind=1:length(H)
  if not(isstruct(H{hind})) %this simulation is missing in action
    mind=mind+1;
    if directory(end)=='/'
      missing(mind).dir =[directory,H{hind}(1:2)];
    else
      missing(mind).dir =[directory,'/',H{hind}(1:2)];
    end
    missing(mind).message=H{hind};
  else
    ind=ind+1;
    goodhinds(ind)=hind;
    out.run(ind).dir       =[H{hind}.dirpath,H{hind}.rundir];
    [out.run(ind).time,out.run(ind).proc]=runtime(out.run(ind).dir);

    %find out if it is single or multi species version
    multi=isfield(H{hind}.run1,'Nspecies');
      
    if not(multi) %only in single species version
      out.RHSMode(ind)      =H{hind}.run1.RHSMode;
      out.normradius(ind)   =H{hind}.run1.normradius;
      out.nuN(ind)          =H{hind}.run1.nuN;
      out.nuPrime(ind)      =H{hind}.run1.nuPrime;
      out.EStar(ind)        =H{hind}.run1.EStar;
      out.dPhiHatdpsi(ind)  =H{hind}.run1.dPhiHatdpsi;
      out.EHat(ind)         =H{hind}.run1.EHat;
      out.nHat(ind)         =H{hind}.run1.nHat;
      out.THat(ind)         =H{hind}.run1.THat;
      if all(out.RHSMode==1)
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

        out.particleFluxBeforeSurfaceIntegral(ind,:,:)=...
            H{hind}.run1.particleFluxBeforeSurfaceIntegral;
        if xoutputexists %for backward compatibility
          if isfield(H{hind}.run1,'fNormIsotropicBeforeSurfaceIntegral')  %for backward compatibility
            out.x{ind}                 =H{hind}.run1.x;
            out.fNormIsotropic{ind}    =H{hind}.run1.fNormIsotropic;
            out.fNormIsotropicBeforeSurfaceIntegral{ind}=H{hind}.run1.fNormIsotropicBeforeSurfaceIntegral;
          else
            out.x{ind}              =NaN*ones(1,out.Nx(ind));
            out.fNormIsotropic{ind} =NaN*ones(1,out.Nx(ind));
            out.fNormIsotropicBeforeSurfaceIntegral{ind}=NaN*ones(out.Ntheta(ind),out.Nzeta(ind),out.Nx(ind));
          end
        end
      end
    else
      out.nu_n(ind)          =H{hind}.run1.nu_n;
      out.dPhiHatdpsi_N(ind) =H{hind}.run1.dPhiHatdpsi_N;
      out.nHats(ind,:)       =H{hind}.run1.nHats;
      out.THats(ind,:)       =H{hind}.run1.THats;
      if all(out.RHSMode==1)
        out.tauhat_m(ind,:)    =H{hind}.run1.NTV;
        out.particleFlux(ind,:) =H{hind}.run1.particleFlux;
        out.heatFlux(ind,:)     =H{hind}.run1.heatFlux;
        out.momentumFlux(ind,:) =H{hind}.run1.momentumFlux;
        out.FSABFlow(ind,:)      =H{hind}.run1.FSABFlow;
      end    
    end    
    
    out.NPeriods(ind)          =double(H{hind}.run1.NPeriods);
    out.Ntheta(ind)            =double(H{hind}.run1.Ntheta);
    out.Nzeta(ind)             =double(H{hind}.run1.Nzeta);
    out.Nxi(ind)               =double(H{hind}.run1.Nxi);
    out.NL(ind)                =double(H{hind}.run1.NL);
    out.Nx(ind)                =double(H{hind}.run1.Nx);
    out.NxPotentialsPerVth(ind)=double(H{hind}.run1.NxPotentialsPerVth);
    out.xMax(ind)              =H{hind}.run1.xMax;
    out.solverTolerance(ind)   =H{hind}.run1.solverTolerance;
    out.theta{ind}             =H{hind}.run1.theta;
    out.theta{ind}             =H{hind}.run1.theta;
    out.zeta{ind}              =H{hind}.run1.zeta;
    out.GHat(ind)          =H{hind}.run1.GHat;
    out.IHat(ind)          =H{hind}.run1.IHat;
    out.iota(ind)          =H{hind}.run1.iota;
    out.B0OverBBar(ind)    =H{hind}.run1.B0OverBBar;
    out.FSABHat2(ind)      =H{hind}.run1.FSABHat2;
    out.didItConverge(ind) =H{hind}.run1.didItConverge;
  end
end
out.NumElements=ind;

if out.NumElements>0
  if not(multi)
    if out.NumElements>0
      if all((out.RHSMode==2)|(out.RHSMode==0)); %The 0 is for a certain type of error
        out.transportMatrix=zeros(out.NumElements,3,3);
        for ind=1:out.NumElements
          out.transportMatrix(ind,:,:)=H{goodhinds(ind)}.run1.transportMatrix;
        end
        if isfield(H{goodhinds(1)}.run1,'NTVMatrix') %For backward compability
          out.NTVMatrix=zeros(out.NumElements,3);
          for ind=1:out.NumElements
            out.NTVMatrix(ind,:)=H{goodhinds(ind)}.run1.NTVMatrix;
          end
        end
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
        elseif strcmp(fnames{find},'NTVMatrix')
          out.NTVMatrix=out.NTVMatrix(orderedInds,:);
        elseif strcmp(fnames{find},'theta')
          out.theta=out.theta{orderedInds};
        elseif strcmp(fnames{find},'zeta')
          out.zeta=out.zeta{orderedInds};
        elseif strcmp(fnames{find},'x')
          out.x=out.x{orderedInds};
        elseif strcmp(fnames{find},'fNormIsotropic')
          out.fNormIsotropic=out.fNormIsotropic{orderedInds};
        elseif strcmp(fnames{find},'fNormIsotropic')
          out.fNormIsotropicBeforeSurfaceIntegral=out.fNormIsotropicBeforeSurfaceIntegral{orderedInds};
        elseif not(strcmp(fnames{find},'NumElements'))
          tmp=getfield(out,fnames{find});
          if size(tmp,1)>1 && size(tmp,2)>1
            out=setfield(out,fnames{find},tmp(orderedInds,:));
          else
            out=setfield(out,fnames{find},tmp(orderedInds));
          end
        end
      end
    end
  end
end
