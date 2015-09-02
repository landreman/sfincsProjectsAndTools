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

for hind=1:length(H)
  if not(isstruct(H{hind})) %this simulation is missing in action
    mind=mind+1;
    if not(iscell(directory))
      if directory(end)=='/'
        missing(mind).dir =[directory,H{hind}(1:2)];
      else
        missing(mind).dir =[directory,'/',H{hind}(1:2)];
      end
      missing(mind).message=H{hind};
    else
      missing(mind).dir =H{hind}(1:end-6);
      missing(mind).message=H{hind};
    end
  elseif H{hind}.RHSMode==0 || ...
        (H{hind}.RHSMode==3 && not(isfield(H{hind},'transportMatrix')))
    %(e.g., or if some other quantity=0. Add more checks here if 
    % necessary)
     mind=mind+1;
    if not(iscell(directory))
      if directory(end)=='/'
        missing(mind).dir =[directory,H{hind}.rundir];
      else
        missing(mind).dir =[directory,'/',H{hind}.rundir];
      end
      missing(mind).message=[H{hind}.rundir,' allval0'];
    else
      missing(mind).dir=H{hind}.rundir;
      missing(mind).message=[H{hind}.rundir,' allval0'];
    end
  else
    ind=ind+1;
    goodhinds(ind)=hind;
    out.run(ind).dir       =[H{hind}.dirpath,H{hind}.rundir];
    out.run(ind).input_namelist  =(H{hind}.input_namelist)';
    %[out.run(ind).time,out.run(ind).proc]=runtime(out.run(ind).dir);

    out.RHSMode(ind)      =double(H{hind}.RHSMode);
    out.Nspecies(ind)     =double(H{hind}.Nspecies);
    Nsp=out.Nspecies(ind);
    out.inputRadialCoordinate(ind)=H{hind}.inputRadialCoordinate;
    out.inputRadialCoordinateForGradients(ind)=H{hind}.inputRadialCoordinateForGradients;
    
    if ind>1
      if out.RHSMode(ind)~=out.RHSMode(ind-1)
        out
        error('Not all have the same RHSMode')
      end
      if out.Nspecies(ind)~=out.Nspecies(ind-1)
        out
        error('Not all have the same Nspecies')
      end
      if out.inputRadialCoordinate(ind)~=out.inputRadialCoordinate(ind-1)
        out
        error('Not all have the same inputRadialCoordinate')
      end
      if out.inputRadialCoordinateForGradients(ind)~=out.inputRadialCoordinateForGradients(ind-1)
        out
        error('Not all have the same inputRadialCoordinateForGradients')
      end
    end
    
    out.theta{ind}             =H{hind}.theta;
    out.theta{ind}             =H{hind}.theta;
    out.zeta{ind}              =H{hind}.zeta;
    out.psiAHat(ind)           =H{hind}.psiAHat;
    out.psiHat(ind)            =H{hind}.psiHat;
    out.rN(ind)            =H{hind}.rN;
    out.GHat(ind)          =H{hind}.GHat;
    out.IHat(ind)          =H{hind}.IHat;
    out.iota(ind)          =H{hind}.iota;
    out.B0OverBBar(ind)    =H{hind}.B0OverBBar;
    out.VPrimeHat(ind)     =H{hind}.VPrimeHat;
    out.FSABHat2(ind)      =H{hind}.FSABHat2;
    out.alpha(ind)         =H{hind}.alpha;
    out.Delta(ind)         =H{hind}.Delta;

    doCorrection=1;
    if doCorrection
      if out.psiAHat(ind)*out.GHat(ind)<0
        %This is to correct an old bug in geometry.F90 for psiAHat
        out.psiAHat(ind)=abs(out.psiAHat(ind))*sign(out.GHat(ind));
        if ind==1
        disp(' WARNING: Signs of G and psiAHat are inconsistent !!!')
        disp('          Changing sign of psiAHat.')  
        end
      end
      if out.psiHat(ind)*out.GHat(ind)<0 
        %This is to correct an old bug in geometry.F90 for psiHat (above it was phiAHat)
        out.psiHat(ind)=abs(out.psiHat(ind))*sign(out.GHat(ind));
        if ind==1
          disp(' WARNING: Signs of G and psiHat are inconsistent !!!')
          disp('          Changing sign of psiHat.')  
        end
      end
    end
    
    if out.RHSMode(ind)==3 %Monoenergetic
      if isfield(H{hind},'nuPrime')
        out.nuPrime(ind)=H{hind}.nuPrime; 
        out.EStar(ind)=H{hind}.EStar;
        out.transportCoeffs(ind,:,:)=H{hind}.transportMatrix;
      else
        warning('nuPrime not stored!!')
      end
    else
      out.Zs(ind,:)          =H{hind}.Zs;
      out.mHats(ind,:)       =H{hind}.mHats;
      out.nHats(ind,:)       =H{hind}.nHats;
      out.THats(ind,:)       =H{hind}.THats;
      out.vTHats(ind,:)      =sqrt(out.THats(ind,:)./out.mHats(ind,:));
      out.dnHatdpsiN(ind,:)  =H{hind}.dnHatdpsiN;
      out.dTHatdpsiN(ind,:)  =H{hind}.dTHatdpsiN;      
      out.dPhiHatdpsiN(ind)  =H{hind}.dPhiHatdpsiN;
      out.EParallelHat(ind)  =H{hind}.EParallelHat;
    end
    if out.RHSMode(ind)==2
       out.transportMatrix(ind,:,:)=H{hind}.transportMatrix;
       out.nu_n(ind)=H{hind}.nu_n;
       if not(isfield(H{hind},'nuPrime'))
         out.nuPrime(ind)=out.nu_n(ind)*...
             (out.GHat(ind)+ out.iota(ind)*out.IHat(ind)) ...
             /sqrt(out.THats(ind,1))/out.B0OverBBar(ind);
       end
       if not(isfield(H{hind},'EStar'))
         out.EStar(ind) = out.GHat(ind)/out.iota(ind)/out.vTHats(ind,1)...
             /out.B0OverBBar(ind)* out.dPhiHatdpsiN(ind)...
             *out.alpha(ind)*out.Delta(ind)/2/out.psiAHat(ind);
       end
    end
    

    if isfield(H{hind},'finished')
      out.finished(ind)               =H{hind}.finished;
    else    
      out.finished(ind)               =0;
    end
    
    if out.RHSMode(ind)==1
      if out.finished(ind)
        if isfield(H{hind},'NTV')
          out.NTV(ind,:)                  =H{hind}.NTV';
          out.particleFlux_vm_psiN(ind,:) =H{hind}.particleFlux_vm_psiN';
          out.heatFlux_vm_psiN(ind,:)     =H{hind}.heatFlux_vm_psiN';
          out.momentumFlux_vm_psiN(ind,:) =H{hind}.momentumFlux_vm_psiN';
          out.FSABFlow(ind,:)             =H{hind}.FSABFlow';
        else
          out.finished(ind)=-1;
          warning(['Run number ',num2str(ind),...
                   ' in the directory ',out.run(ind).dir ,...
                   ' says ''finished'' but output results are missing']);
          out.NTV(ind,:)                  =NaN*zeros(Nsp,1);
          out.particleFlux_vm_psiN(ind,:) =NaN*zeros(Nsp,1);
          out.heatFlux_vm_psiN(ind,:)     =NaN*zeros(Nsp,1);
          out.momentumFlux_vm_psiN(ind,:) =NaN*zeros(Nsp,1);
          out.FSABFlow(ind,:)             =NaN*zeros(Nsp,1);
        end
      else
        out.NTV(ind,:)                  =NaN*ones(out.Nspecies(ind),1);
        out.particleFlux_vm_psiN(ind,:) =NaN*ones(out.Nspecies(ind),1);
        out.heatFlux_vm_psiN(ind,:)     =NaN*ones(out.Nspecies(ind),1);
        out.momentumFlux_vm_psiN(ind,:) =NaN*ones(out.Nspecies(ind),1);
        out.FSABFlow(ind,:)             =NaN*ones(out.Nspecies(ind),1);
      end
    end    
    
    
    out.NPeriods(ind)          =double(H{hind}.NPeriods);
    out.Ntheta(ind)            =double(H{hind}.Ntheta);
    out.Nzeta(ind)             =double(H{hind}.Nzeta);
    out.Nxi(ind)               =double(H{hind}.Nxi);
    out.NL(ind)                =double(H{hind}.NL);
    if out.RHSMode(ind)~=3
      out.Nx(ind)                =double(H{hind}.Nx);
      out.NxPotentialsPerVth(ind)=double(H{hind}.NxPotentialsPerVth);
      out.xMax(ind)              =H{hind}.xMax;
      out.x{ind}                 =H{hind}.x;
    end
    out.solverTolerance(ind)   =H{hind}.solverTolerance;
    %out.didItConverge(ind) =H{hind}.didItConverge; %old name
    %out.finished(ind) =H{hind}.finished; %new name?
    if out.RHSMode(ind)==1 &&  out.finished(ind)==1
      %out.NTVtot(ind)=sum(out.NTV(ind,:));
      %out.NTVfromFlux(ind,:)=1e-3*out.iota(ind)*...
      %      sqrt(2*1.6022e-19*1e3/1.6726e-27)* out.psiAHat(ind)*...
      %      out.particleFlux_vm_psiN(ind,:).*out.Zs(ind,:);%is equal to:
      out.NTVfromFlux(ind,:)=out.psiAHat(ind)*out.iota(ind)*...
          out.Zs(ind,:).*out.particleFlux_vm_psiN(ind,:)/(out.Delta(ind)/2);
    else
      out.NTVfromFlux(ind,1:out.Nspecies(ind))=NaN;
    end
  end
end
    
out.NumElements=ind;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%{    
    
    
    
    %find out if it is single or multi species version
    multi=isfield(H{hind},'Nspecies');
      
    % The following is sadly not stored as output other than in 
    % H{hind}.input_namelist'
    %out.geometryScheme(ind) =H{hind}.geometryScheme;
    %if out.geometryScheme(ind)==11
    %  out.JGboozer_file{ind} = H{hind}.JGboozer_file;
    %elseif out.geometryScheme(ind)==12
    %  out.JGboozer_file_NonStelSym{ind}=H{hind}.JGboozer_file_NonStelSym;
    %end
    if not(multi) %only in single species version
      out.RHSMode(ind)      =H{hind}.RHSMode;
      out.normradius(ind)   =H{hind}.normradius;
      out.nuN(ind)          =H{hind}.nuN;
      out.nuPrime(ind)      =H{hind}.nuPrime;
      out.EStar(ind)        =H{hind}.EStar;
      out.dPhiHatdpsi(ind)  =H{hind}.dPhiHatdpsi;
      out.EHat(ind)         =H{hind}.EHat;
      out.nHat(ind)         =H{hind}.nHat;
      out.THat(ind)         =H{hind}.THat;
      if all(out.RHSMode==1)
        if isfield(H{hind},'NTV') %Only for backward compability
          out.tauhat_s(ind)     =H{hind}.NTV;
        else
          out.tauhat_s(ind)     =H{hind}.NTVsingle;
          out.tauhat_m(ind)     =H{hind}.NTVmulti;
        end  
        out.particleFlux(ind) =H{hind}.particleFlux;
        out.heatFlux(ind)     =H{hind}.heatFlux;
        out.momentumFlux(ind) =H{hind}.momentumFlux;
        out.FSAFlow(ind)      =H{hind}.FSAFlow;

        out.particleFluxBeforeSurfaceIntegral(ind,:,:)=...
            H{hind}.particleFluxBeforeSurfaceIntegral;
        if xoutputexists %for backward compatibility
          if isfield(H{hind},'fNormIsotropicBeforeSurfaceIntegral')  %for backward compatibility
            out.x{ind}                 =H{hind}.x;
            out.fNormIsotropic{ind}    =H{hind}.fNormIsotropic;
            out.fNormIsotropicBeforeSurfaceIntegral{ind}=H{hind}.fNormIsotropicBeforeSurfaceIntegral;
          else
            out.x{ind}              =NaN*ones(1,out.Nx(ind));
            out.fNormIsotropic{ind} =NaN*ones(1,out.Nx(ind));
            out.fNormIsotropicBeforeSurfaceIntegral{ind}=NaN*ones(out.Ntheta(ind),out.Nzeta(ind),out.Nx(ind));
          end
        end
      end
    else
      out.nu_n(ind)          =H{hind}.nu_n;
      out.dPhiHatdpsi_N(ind) =H{hind}.dPhiHatdpsi_N;
      out.nHats(ind,:)       =H{hind}.nHats;
      out.THats(ind,:)       =H{hind}.THats;
      if all(out.RHSMode==1)
        out.tauhat_m(ind,:)    =H{hind}.NTV;
        out.particleFlux(ind,:) =H{hind}.particleFlux;
        out.heatFlux(ind,:)     =H{hind}.heatFlux;
        out.momentumFlux(ind,:) =H{hind}.momentumFlux;
        out.FSABFlow(ind,:)      =H{hind}.FSABFlow;
      end    
    end    
    
    out.NPeriods(ind)          =double(H{hind}.NPeriods);
    out.Ntheta(ind)            =double(H{hind}.Ntheta);
    out.Nzeta(ind)             =double(H{hind}.Nzeta);
    out.Nxi(ind)               =double(H{hind}.Nxi);
    out.NL(ind)                =double(H{hind}.NL);
    out.Nx(ind)                =double(H{hind}.Nx);
    out.NxPotentialsPerVth(ind)=double(H{hind}.NxPotentialsPerVth);
    out.xMax(ind)              =H{hind}.xMax;
    out.solverTolerance(ind)   =H{hind}.solverTolerance;
    out.theta{ind}             =H{hind}.theta;
    out.theta{ind}             =H{hind}.theta;
    out.zeta{ind}              =H{hind}.zeta;
    out.GHat(ind)          =H{hind}.GHat;
    out.IHat(ind)          =H{hind}.IHat;
    out.iota(ind)          =H{hind}.iota;
    out.B0OverBBar(ind)    =H{hind}.B0OverBBar;
    out.FSABHat2(ind)      =H{hind}.FSABHat2;
    out.didItConverge(ind) =H{hind}.didItConverge;
  end
end
out.NumElements=ind;

if out.NumElements>0
  if not(multi)
      if all((out.RHSMode==2)|(out.RHSMode==0)); %The 0 is for a certain type of error
        out.transportMatrix=zeros(out.NumElements,3,3);
        for ind=1:out.NumElements
          out.transportMatrix(ind,:,:)=H{goodhinds(ind)}.transportMatrix;
        end
        if isfield(H{goodhinds(1)},'NTVMatrix') %For backward compability
          out.NTVMatrix=zeros(out.NumElements,3);
          for ind=1:out.NumElements
            out.NTVMatrix(ind,:)=H{goodhinds(ind)}.NTVMatrix;
          end
        end
      elseif all(out.RHSMode==3)
        out.transportCoeffs=zeros(out.NumElements,2,2);        
        for ind=1:out.NumElements
          out.transportCoeffs(ind,:,:)=H{goodhinds(ind)}.transportCoeffs;
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
        elseif strcmp(fnames{find},'transportCoeffs')
          out.transportCoeffs=out.transportCoeffs(orderedInds,:,:);
        elseif strcmp(fnames{find},'theta')
          out.theta={out.theta{orderedInds}};
        elseif strcmp(fnames{find},'zeta')
          out.zeta={out.zeta{orderedInds}};
        elseif strcmp(fnames{find},'x')
          out.x={out.x{orderedInds}};
        elseif strcmp(fnames{find},'fNormIsotropic')
          out.fNormIsotropic={out.fNormIsotropic{orderedInds}};
        elseif strcmp(fnames{find},'fNormIsotropic')
          out.fNormIsotropicBeforeSurfaceIntegral={out.fNormIsotropicBeforeSurfaceIntegral{orderedInds}};
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


%}