function [out,missing]=getresults(directory,sortafter,varargin)
% Retrieves the output from the sfincs runs with the same discretisation
% in the numbered subdirectories in "directory". 
% Only succesful runs are loaded. The output struct "out"
% is sorted after the vaiable with the name in the input "sortafter"
%directory

load_unfinished_too=0;
if nargin>2
  if strcmp(varargin{1},'load unfinished too')
    load_unfinished_too=1;
  end
end

if nargin==2
  H=loadallh5(directory,sortafter);
else
  H=loadallh5(directory);
end
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
    
    out.NPeriods(ind)          =double(H{hind}.NPeriods);
    out.Ntheta(ind)            =double(H{hind}.Ntheta);
    out.Nzeta(ind)             =double(H{hind}.Nzeta);
    out.Nxi(ind)               =double(H{hind}.Nxi);
    out.NL(ind)                =double(H{hind}.NL);
    out.includePhi1(ind)                =double(H{hind}.includePhi1);

    out.theta{ind}             =H{hind}.theta;
    out.theta{ind}             =H{hind}.theta;
    out.zeta{ind}              =H{hind}.zeta;
    out.psiAHat(ind)           =H{hind}.psiAHat;
    out.psiHat(ind)            =H{hind}.psiHat;
    out.rN(ind)            =H{hind}.rN;
    out.GHat(ind)          =H{hind}.GHat;
    out.IHat(ind)          =H{hind}.IHat;
    out.B0OverBBar(ind)    =H{hind}.B0OverBBar;
    out.iota(ind)          =H{hind}.iota;
    out.VPrimeHat(ind)     =H{hind}.VPrimeHat;
    out.FSABHat2(ind)      =H{hind}.FSABHat2;
    out.alpha(ind)         =H{hind}.alpha;
    out.Delta(ind)         =H{hind}.Delta;
    out.nu_n(ind)          =H{hind}.nu_n;

    %out.BDotCurlB{ind}     =H{hind}.BDotCurlB; %Only saved for VMEC input ????
    out.BHat{ind}          =H{hind}.BHat';
    out.dBHatdtheta{ind}   =H{hind}.dBHatdtheta';
    out.dBHatdzeta{ind}    =H{hind}.dBHatdzeta';
    out.BHat_sub_psi{ind}   =H{hind}.BHat_sub_psi';
    out.BHat_sup_theta{ind}   =H{hind}.BHat_sup_theta';
    out.BHat_sup_zeta{ind}   =H{hind}.BHat_sup_zeta';
    out.dBHat_sub_psi_dtheta{ind}   =H{hind}.dBHat_sub_psi_dtheta';
    out.dBHat_sub_psi_dzeta{ind}   =H{hind}.dBHat_sub_psi_dzeta';
    out.dIHat_dpsiHat{ind}   =H{hind}.dBHat_sub_theta_dpsiHat';
    out.dGHat_dpsiHat{ind}   =H{hind}.dBHat_sub_zeta_dpsiHat';
    out.dBHat_sup_theta_dpsiHat{ind}   =H{hind}.dBHat_sup_theta_dpsiHat';
    out.dBHat_sup_theta_dzeta{ind}   =H{hind}.dBHat_sup_theta_dzeta';
    out.dBHat_sup_zeta_dpsiHat{ind}   =H{hind}.dBHat_sup_zeta_dpsiHat';
    out.dBHat_sup_zeta_dtheta{ind}   =H{hind}.dBHat_sup_zeta_dtheta';
    out.dBHatdpsiHat{ind}   =H{hind}.dBHatdpsiHat';    
    
    if isfield(H{hind},'uHat')
      out.uHat{ind}         =H{hind}.uHat';
    end
    
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
      out.dnHatdrN(ind,:)    =H{hind}.dnHatdrN;
      out.dnHatdrHat(ind,:)  =H{hind}.dnHatdrHat;
      out.dTHatdpsiN(ind,:)  =H{hind}.dTHatdpsiN;      
      out.dTHatdrN(ind,:)    =H{hind}.dTHatdrN;      
      out.dTHatdrHat(ind,:)  =H{hind}.dTHatdrHat;      
      out.dPhiHatdpsiN(ind)  =H{hind}.dPhiHatdpsiN;
      out.dPhiHatdrN(ind)    =H{hind}.dPhiHatdrN;
      out.dPhiHatdrHat(ind)  =H{hind}.dPhiHatdrHat;
      out.EParallelHat(ind)  =H{hind}.EParallelHat;
      out.withAdiabatic(ind) =(H{hind}.withAdiabatic==H{hind}.integerToRepresentTrue);
      if out.withAdiabatic(ind)
        out.adiabaticZ(ind)=H{hind}.adiabaticZ;
        out.adiabaticMHat(ind)=H{hind}.adiabaticMHat;
        out.adiabaticNHat(ind)=H{hind}.adiabaticNHat;
        out.adiabaticTHat(ind)=H{hind}.adiabaticTHat;
      end
      if isfield(H{hind},'withNBIspec')
        out.withNBIspec(ind)   =(H{hind}.withNBIspec==H{hind}.integerToRepresentTrue);
        if out.withNBIspec(ind)
          out.NBIspecZ(ind)=H{hind}.NBIspecZ;
          out.NBIspecNHat(ind)=H{hind}.NBIspecNHat;
        end
      else
        out.withNBIspec(ind)=0;
      end
      out.includePhi1(ind)   =(H{hind}.includePhi1==H{hind}.integerToRepresentTrue);
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
      if out.finished(ind) || (load_unfinished_too)% && out.includePhi1(ind))
        if isfield(H{hind},'NTV')
          if out.includePhi1(ind)==0
            out.NTV(ind,:)                  =H{hind}.NTV';
            out.particleFlux_vm_rHat(ind,:) =H{hind}.particleFlux_vm_rHat';
            out.particleFlux_vm_rN(ind,:)   =H{hind}.particleFlux_vm_rN';
            out.particleFlux_vm_psiN(ind,:) =H{hind}.particleFlux_vm_psiN';
            out.particleFlux_vm0_psiN(ind,:)=H{hind}.particleFlux_vm0_psiN';
            out.heatFlux_vm_rHat(ind,:)     =H{hind}.heatFlux_vm_rHat';
            out.heatFlux_vm_rN(ind,:)       =H{hind}.heatFlux_vm_rN';
            out.heatFlux_vm_psiN(ind,:)     =H{hind}.heatFlux_vm_psiN';
            out.heatFlux_vm0_psiN(ind,:)    =H{hind}.heatFlux_vm0_psiN';
            out.momentumFlux_vm_rHat(ind,:) =H{hind}.momentumFlux_vm_rHat';
            out.momentumFlux_vm_rN(ind,:)   =H{hind}.momentumFlux_vm_rN';
            out.momentumFlux_vm_psiN(ind,:) =H{hind}.momentumFlux_vm_psiN';
            out.FSABFlow(ind,:)             =H{hind}.FSABFlow';
            out.FSABjHat(ind,1)             =H{hind}.FSABjHat';
            if isfield(H{hind},'classicalParticleFlux_rHat')
              out.classicalParticleFlux_rHat(ind,:)=H{hind}.classicalParticleFlux_rHat';
              out.classicalParticleFlux_rN(ind,:)  =H{hind}.classicalParticleFlux_rN';
              out.classicalParticleFlux_psiN(ind,:)=H{hind}.classicalParticleFlux_psiN';
            end
            if Nsp>1
              out.flow{ind}                    =permute(H{hind}.flow,[3,2,1]);
              out.densityPerturbation{ind}     =permute(H{hind}.densityPerturbation,[3,2,1]);
              out.pressurePerturbation{ind}    =permute(H{hind}.pressurePerturbation,[3,2,1]);
              out.pressureAnisotropy{ind}      =permute(H{hind}.pressureAnisotropy,[3,2,1]);
              out.NTVBeforeSurfaceIntegral{ind}=permute(H{hind}.NTVBeforeSurfaceIntegral,[3,2,1]);
            else
              out.flow{ind}                    =H{hind}.flow';
              out.densityPerturbation{ind}     =H{hind}.densityPerturbation';
              out.pressurePerturbation{ind}    =H{hind}.pressurePerturbation';
              out.pressureAnisotropy{ind}      =H{hind}.pressureAnisotropy';
              out.NTVBeforeSurfaceIntegral{ind}=H{hind}.NTVBeforeSurfaceIntegral';
            end
          else %includePhi1=1
            out.NTV(ind,:)                  =H{hind}.NTV(:,end)';
            out.particleFlux_vd_rHat(ind,:) =H{hind}.particleFlux_vd_rHat(:,end)';
            out.particleFlux_vd_rN(ind,:)   =H{hind}.particleFlux_vd_rN(:,end)';
            out.particleFlux_vd_psiN(ind,:) =H{hind}.particleFlux_vd_psiN(:,end)';
            out.heatFlux_vd_rHat(ind,:)     =H{hind}.heatFlux_vd_rHat(:,end)';
            out.heatFlux_vd_rN(ind,:)       =H{hind}.heatFlux_vd_rN(:,end)';
            out.heatFlux_vd_psiN(ind,:)     =H{hind}.heatFlux_vd_psiN(:,end)';
            out.momentumFlux_vd_rHat(ind,:) =H{hind}.momentumFlux_vd_rHat(:,end)';
            out.momentumFlux_vd_rN(ind,:)   =H{hind}.momentumFlux_vd_rN(:,end)';
            out.momentumFlux_vd_psiN(ind,:) =H{hind}.momentumFlux_vd_psiN(:,end)';
            out.FSABFlow(ind,:)             =H{hind}.FSABFlow(:,end)';
            out.FSABjHat(ind,1)             =H{hind}.FSABjHat(end)';
            if isfield(H{hind},'classicalParticleFlux_rHat')
              out.classicalParticleFlux_rHat(ind,:)=H{hind}.classicalParticleFlux_rHat(:,end)';
              out.classicalParticleFlux_rN(ind,:)  =H{hind}.classicalParticleFlux_rN(:,end)';
              out.classicalParticleFlux_psiN(ind,:)=H{hind}.classicalParticleFlux_psiN(:,end)';
            else
              out.classicalParticleFlux_rHat(ind,:)=NaN*zeros(Nsp,1);
              out.classicalParticleFlux_rN(ind,:)  =NaN*zeros(Nsp,1);
              out.classicalParticleFlux_psiN(ind,:)=NaN*zeros(Nsp,1);
            end
            out.Phi1Hat{ind}                =squeeze(H{hind}.Phi1Hat(:,:,end));
            if Nsp>1
              out.flow{ind}                    =permute(squeeze(H{hind}.flow(:,:,:,end)),[3,2,1]);
              out.densityPerturbation{ind}     =permute(squeeze(H{hind}.densityPerturbation(:,:,:,end)),[3,2,1]);
              out.pressurePerturbation{ind}    =permute(squeeze(H{hind}.pressurePerturbation(:,:,:,end)),[3,2,1]);
              out.pressureAnisotropy{ind}      =permute(squeeze(H{hind}.pressureAnisotropy(:,:,:,end)),[3,2,1]);
              out.NTVBeforeSurfaceIntegral{ind}=permute(squeeze(H{hind}.NTVBeforeSurfaceIntegral(:,:,:,end)),[3,2,1]);
            else
              out.flow{ind}                    =squeeze(H{hind}.flow(:,:,end))';
              out.densityPerturbation{ind}     =squeeze(H{hind}.densityPerturbation(:,:,end))';
              out.pressurePerturbation{ind}    =squeeze(H{hind}.pressurePerturbation(:,:,end))';
              out.pressureAnisotropy{ind}      =squeeze(H{hind}.pressureAnisotropy(:,:,end))';
              out.NTVBeforeSurfaceIntegral{ind}=squeeze(H{hind}.NTVBeforeSurfaceIntegral(:,:,end))';
              out.dPhi1Hatdtheta{ind}          =squeeze(H{hind}.dPhi1Hatdtheta(:,:,end))';
              out.dPhi1Hatdzeta{ind}           =squeeze(H{hind}.dPhi1Hatdzeta(:,:,end))';
            end
          end
        else
          out.finished(ind)=-1;
          warning(['Run number ',num2str(ind),...
                   ' in the directory ',out.run(ind).dir ,...
                   ' says ''finished'' but output results are missing']);
          out.NTV(ind,:)                 =NaN*zeros(Nsp,1);
          if out.includePhi1(ind)==1
            out.Phi1Hat{ind}=NaN;
            out.particleFlux_vd_rHat(ind,:)=NaN*zeros(Nsp,1);
            out.particleFlux_vd_rN(ind,:)  =NaN*zeros(Nsp,1);
            out.particleFlux_vd_psiN(ind,:)=NaN*zeros(Nsp,1);
            out.heatFlux_vd_rN(ind,:)      =NaN*zeros(Nsp,1);
            out.heatFlux_vd_rHat(ind,:)    =NaN*zeros(Nsp,1);
            out.heatFlux_vd_psiN(ind,:)    =NaN*zeros(Nsp,1);
            out.momentumFlux_vd_rHat(ind,:)=NaN*zeros(Nsp,1);
            out.momentumFlux_vd_rN(ind,:)  =NaN*zeros(Nsp,1);
            out.momentumFlux_vd_psiN(ind,:)=NaN*zeros(Nsp,1);            
          else
            out.particleFlux_vm_rHat(ind,:)=NaN*zeros(Nsp,1);
            out.particleFlux_vm_rN(ind,:)  =NaN*zeros(Nsp,1);
            out.particleFlux_vm_psiN(ind,:)=NaN*zeros(Nsp,1);
            out.particleFlux_vm0_psiN(ind,:)=NaN*zeros(Nsp,1);
            out.heatFlux_vm_rN(ind,:)      =NaN*zeros(Nsp,1);
            out.heatFlux_vm_rHat(ind,:)    =NaN*zeros(Nsp,1);
            out.heatFlux_vm_psiN(ind,:)    =NaN*zeros(Nsp,1);
            out.momentumFlux_vm_rHat(ind,:)=NaN*zeros(Nsp,1);
            out.momentumFlux_vm_rN(ind,:)  =NaN*zeros(Nsp,1);
            out.momentumFlux_vm_psiN(ind,:)=NaN*zeros(Nsp,1);
          end
          out.FSABFlow(ind,:)            =NaN*zeros(Nsp,1);
          out.FSABjHat(ind,1)            =NaN;
          out.flow{ind}                  =NaN*zeros(Nsp,out.Ntheta(ind),out.Nzeta(ind));
          out.densityPerturbation{ind}   =NaN*zeros(Nsp,out.Ntheta(ind),out.Nzeta(ind));
          out.pressurePerturbation{ind}  =NaN*zeros(Nsp,out.Ntheta(ind),out.Nzeta(ind));
          out.pressureAnisotropy{ind}    =NaN*zeros(Nsp,out.Ntheta(ind),out.Nzeta(ind));
          out.NTVBeforeSurfaceIntegral{ind}=NaN*zeros(Nsp,out.Ntheta(ind),out.Nzeta(ind));
          %if isfield(H{hind},'classicalParticleFlux_rHat')
          out.classicalParticleFlux_rHat(ind,:)=NaN*zeros(Nsp,1);
          out.classicalParticleFlux_rN(ind,:)  =NaN*zeros(Nsp,1);
          out.classicalParticleFlux_psiN(ind,:)=NaN*zeros(Nsp,1);
          %end
        end
      else
        if out.includePhi1(ind)==1
          out.Phi1Hat{ind}=NaN;
          out.particleFlux_vd_rHat(ind,:) =NaN*ones(out.Nspecies(ind),1);
          out.particleFlux_vd_psiN(ind,:) =NaN*ones(out.Nspecies(ind),1);
          out.particleFlux_vd_rN(ind,:)   =NaN*ones(out.Nspecies(ind),1);
          out.heatFlux_vd_psiN(ind,:)     =NaN*ones(out.Nspecies(ind),1);
          out.heatFlux_vd_rHat(ind,:)     =NaN*ones(out.Nspecies(ind),1);
          out.heatFlux_vd_rN(ind,:)     =NaN*ones(out.Nspecies(ind),1);
          out.momentumFlux_vd_psiN(ind,:) =NaN*ones(out.Nspecies(ind),1);
          out.momentumFlux_vd_rHat(ind,:) =NaN*ones(out.Nspecies(ind),1);
          out.momentumFlux_vd_psiN(ind,:) =NaN*ones(out.Nspecies(ind),1);
        else  
          out.particleFlux_vm_rHat(ind,:) =NaN*ones(out.Nspecies(ind),1);
          out.particleFlux_vm_rN(ind,:)   =NaN*ones(out.Nspecies(ind),1);
          out.particleFlux_vm_psiN(ind,:) =NaN*ones(out.Nspecies(ind),1);
          out.particleFlux_vm0_psiN(ind,:) =NaN*ones(out.Nspecies(ind),1);
          out.heatFlux_vm_psiN(ind,:)     =NaN*ones(out.Nspecies(ind),1);
          out.heatFlux_vm_rHat(ind,:)     =NaN*ones(out.Nspecies(ind),1);
          out.heatFlux_vm_rN(ind,:)       =NaN*ones(out.Nspecies(ind),1);
          out.momentumFlux_vm_psiN(ind,:) =NaN*ones(out.Nspecies(ind),1);
          out.momentumFlux_vm_rHat(ind,:) =NaN*ones(out.Nspecies(ind),1);
          out.momentumFlux_vm_psiN(ind,:) =NaN*ones(out.Nspecies(ind),1);
        end
        out.NTV(ind,:)                  =NaN*ones(out.Nspecies(ind),1);
        out.FSABFlow(ind,:)             =NaN*ones(out.Nspecies(ind),1);
        out.FSABjHat(ind,1)             =NaN;
        out.flow{ind}                   =NaN*zeros(Nsp,out.Ntheta(ind),out.Nzeta(ind));
        out.densityPerturbation{ind}    =NaN*zeros(Nsp,out.Ntheta(ind),out.Nzeta(ind));
        out.pressurePerturbation{ind}   =NaN*zeros(Nsp,out.Ntheta(ind),out.Nzeta(ind));
        out.pressureAnisotropy{ind}     =NaN*zeros(Nsp,out.Ntheta(ind),out.Nzeta(ind));
        out.NTVBeforeSurfaceIntegral{ind}=NaN*zeros(Nsp,out.Ntheta(ind),out.Nzeta(ind));
        %if isfield(H{hind},'classicalParticleFlux_rHat')
        out.classicalParticleFlux_rHat(ind,:)=NaN*zeros(out.Nspecies(ind),1);
        out.classicalParticleFlux_rN(ind,:)  =NaN*zeros(out.Nspecies(ind),1);
        out.classicalParticleFlux_psiN(ind,:)=NaN*zeros(out.Nspecies(ind),1);
        %end
      end
    end    
    
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
      if out.includePhi1(ind)==1
        out.NTVfromFlux(ind,:)=out.psiAHat(ind)*out.iota(ind)*...
            out.Zs(ind,:).*out.particleFlux_vd_psiN(ind,:)/(out.Delta(ind)/2);
      else
        out.NTVfromFlux(ind,:)=out.psiAHat(ind)*out.iota(ind)*...
            out.Zs(ind,:).*out.particleFlux_vm_psiN(ind,:)/(out.Delta(ind)/2);        
      end
    else
      out.NTVfromFlux(ind,1:out.Nspecies(ind))=NaN;
    end
  end
end
    
out.NumElements=ind;
    
    
    
    
    
    
    
    
    
    
    
    