function sfincsAlbert()

% SFINCS:
% The Stellarator Fokker-Planck Iterative Neoclassical Conservative Solver.
% Single species version.
% Original version written in 2013 by Matt Landreman
% Massachusetts Institute of Technology
% Plasma Science & Fusion Center
    
% Dimensional quantities in this program are normalized to "reference" values:
% \bar{B} = reference magnetic field, typically 1 Tesla.
% \bar{R} = reference length, typically 1 meter.
% \bar{n} = reference density, typically 10^19 m^{-3}, 10^20 m^{-3}, or something similar.
% \bar{m} = reference mass, typically either the mass of hydrogen or deuterium.
% \bar{T} = reference temperature in energy units, typically 1 eV or 1 keV.
% \bar{v} = \sqrt{2 * \bar{T} / \bar{m}} = reference speed
% \bar{Phi} = reference electrostatic potential, typically 1 V or 1 kV.

% You can choose any reference parameters you like, not just the values
% suggested here. The code "knows" about the reference values only through
% the 3 combinations Delta, omega, and nuN or nuPrime, input below.

% Radial gradients of density, temperature, and electrostatic potential are
% specified as derivatives with respect to psi_N, where psi_N is the 
% toroidal flux normalized to the value at the last closed flux surface. 
% (psi_N=0 is the magnetic axis, and psi_N=1 is the last closed flux 
% surface.)
    
% --------------------------------------------------
% Program flow control parameters:
% --------------------------------------------------

programMode = 2;
% 1 = single run.
% 2 = Do a convergence scan and save the results.
% 3 = Load a previous convergence scan and plot the results. (Doesn't do any new solves.)
% 4 = Do a nuPrime scan and save the results.
% 5 = Load a previous nuPrime scan and plot the results. (Doesn't do any new solves.)

% The setting below matters for programMode=3 or programMode=5 only:
dataFileToPlot = 'm20130318_02_SFINCS_2013-03-18_14-47_convergenceScan_convergenceScan.mat';

RHSMode = 1;
% 1 = Use a single right-hand side.
% 2 = Use multiple right-hand sides to compute the transport matrix.

% The variable below is set to true only for rare testing
% circumstances. Typically it should be false.
%testQuasisymmetryIsomorphism = true;
testQuasisymmetryIsomorphism = false;

%saveStuff = true;
saveStuff = false;

% The string below is appended to the filename of the data file
% (but before the .mat extension.)
filenameNote = 'EUTERPE_benchmark';

% --------------------------------------------------
% Geometry parameters:
% --------------------------------------------------

geometryScheme = 1;
% 1 = Two-helicity model
% 2 = Three-helicity approximation of the LHD standard configuration
% 3 = Four-helicity approximation of the LHD inward-shifted configuration
% 4 = Three-helicity approximation of the W7-X Standard configuration
% 10= Read the boozer coordinate data from the file specified as "fort996boozer_file" below
% 11= Read the boozer coordinate data from the file specified as "JGboozer_file" below

% Additional parameters used only when geometryScheme=1:
% B = BBar * B0OverBBar * [1 + epsilon_t * cos(theta) + epsilon_h * cos(helicity_l * theta - helicity_n * zeta)]
B0OverBBar = 1.55;
epsilon_t = -0.178539;
%epsilon_t = 0.178539;
epsilon_h = 0.0;
helicity_l = 0;
helicity_n = 0;
% iota is the rotational transform = 1 / (safety factor q)
iota = 0.397543;
%iota = 0.397543;
% G is c/2 * the poloidal current outside the flux
% surface. Equivalently, G is the coefficient of grad zeta in the
% covariant representation of vector B. GHat is G normalized by \bar{B}\bar{R}.
%GHat = 6.2;
GHat = 6.2;
% I is c/2 * the toroidal current inside the flux
% surface. Equivalently, I is the coefficient of grad theta in the
% covariant representation of vector B. IHat is I normalized by \bar{B}\bar{R}.
IHat = 0.0785669;
%IHat = 0.0785669;

% geometryScheme=10 parameters
fort996boozer_file='TJII-midradius_example_s_0493_fort.996';
% Note that PsiA is not stored in the fort.996 file, so we use the
% PsiAHat setting below

% geometryScheme=11 parameters
JGboozer_file='w7x-sc1.bc';
normradius_wish=0.5;   %The calculation will be performed for the radius
                       %closest to this one in the JGboozer_file
min_Bmn_to_load=1e-2;  %Filter out any Bmn components smaller than this

% --------------------------------------------------
% Physics parameters:
% --------------------------------------------------
% Roughly speaking, Delta is rho_* at the reference parameters.
% More precisely, 
% Delta = c * \bar{m} * \bar{v} / (e * \bar{B} * \bar{R}) in Gaussian units,
% Delta =     \bar{m} * \bar{v} / (e * \bar{B} * \bar{R}) in SI units,
% where
% c = speed of light
% e = proton charge

Delta_e = 1.0664e-4; %electrons, reference values: \bar{T}=1 keV, \bar{n}=10^20 m^-3,
                     %\bar{Phi}=1 kV, \bar{B}=1 T, \bar{R}=1 m
Delta_p = 4.5694e-3; %protons, reference values: \bar{T}=1 keV, \bar{n}=10^20 m^-3,
                     %\bar{Phi}=1 kV, \bar{B}=1 T, \bar{R}=1 m

omega_e = 5.3318e-5; %electrons, reference values: \bar{T}=1 keV, \bar{n}=10^20 m^-3,
                     %\bar{Phi}=1 kV, \bar{B}=1 T, \bar{R}=1 m
omega_p = 2.2847e-3; %protons, reference values: \bar{T}=1 keV, \bar{n}=10^20 m^-3,
                     %\bar{Phi}=1 kV, \bar{B}=1 T, \bar{R}=1 m

nu_nbar_pp=8.4774e-3; 
nu_nbar_ee=7.3013e-3;

species='p';
if species=='p'
  Delta   = Delta_p;
  omega   = omega_p;
  vbar    = sqrt(2*1e3*1.6022e-19/1.6726e-27);
  nu_nbar = nu_nbar_pp;
elseif species=='e'
  Delta = Delta_e;
  omega = omega_e;
  vbar    = sqrt(2*1e3*1.6022e-19/9.1094e-31);
  nu_nbar = nu_nbar_ee;  
  disp('Warning: only electron-electron collisions included!')
end

% psiAHat = psi_a / (\bar{B} * \bar{R}^2) (in both Gaussian and SI units)
% where 2*pi*psi_a is the toroidal flux at the last closed flux surface
% (the surface where psi_N = 1.)
% The value of psiAHat here is over-written for geometryScheme = 2, 3, 4, or 11.
psiAHat = 0.323266;
THat = 2.0;
nHat = 0.0596712;

% Radial electric field.
dPhiHatdpsi = 0;

% The following two quantities matter for RHSMode=1 but not for RHSMode=2:
dTHatdpsi = 0;
dnHatdpsi = -0.202428;
EHat = 0;

% There are 2 different ways to specify the collisionality: nuN and nuPrime.
% If RHSMode == 1, nuN is used and nuPrime is ignored.
% If RHSMode == 2, nuPrime is used and nuN is ignored.
% 
% nuN = nu_ii * \bar{R} / \bar{v}
% and
% nuPrime = nu_ii * (G + iota * I) / (v_i * B_0)
%         = BBarOverB0 / sqrt(THat) * (GHat + iota * IHat) * nuN
%
% where
% v_i = sqrt(2 * T_i / m_i) and
%
%                  4 * sqrt{2*pi} * n_i * Z^4 * e^4 * ln(Lambda)
% nu_ii = -----------------------------------------------------------   (SI units)
%             3 * (4 * pi * epsilon_0)^2 * sqrt(m_i} * T_i^(3/2)
%
% or, equivalently,
%
%                  4 * sqrt{2*pi} * n_i * Z^4 * e^4 * ln(Lambda)
% nu_ii = -----------------------------------------------------------   (Gaussian units)
%                       3 * sqrt(m_i} * T_i^(3/2)
%
% The definition of nuPrime is motivated by the fact that the
% transport matrix elements depend on the density and temperature
% only through nuPrime, not individually. Hence, nuPrime is used as
% the measure of collisionality when RHSMode=2. However, 
% the code originally used the different collisionality
% definition nuN. Hence, for historical reasons, nuN is used instead of nuPrime when RHSMode=1.
%
% Notice that collisionality is defined differently in the multi-species code!

nuN =  nu_nbar * nHat/THat^(3/2);
%nuhat_Helander = 1e-5; nuN = 4/3/sqrt(pi)/vbar*nuhat_Helander; %nuhat as in eq(3.48) Helander&Sigmar 

% If testQuasisymmetryIsomorphism is true, the value of nuN is changed so the physical collisionality
% stays constant as the helicity is changed.

nuPrime = NaN;%1;  %This value is over-written for geometryScheme = 4, 11

collisionOperator = 1;
% 0 = Full linearized Fokker-Planck operator
% 1 = Pitch angle scattering, with no momentum conservation
% 2 = Pitch angle scattering, with a model momentum-conserving field term
%Added by AM
% 3 = Pitch angle scattering, with no momentum conservation, and
% monoenergetic nuD

%Added by AM
%Monoenergetic deflection frequency, only used for collisionOperator = 3
nuDconst = 0.1;

% Unless you know what you are doing, keep constraintScheme = -1.
constraintScheme = -1;
% -1 = Automatic: if collisionOperator==0 then use constraintScheme=1, otherwise use constraintScheme=2.
%  0 = No constraints
%  1 = 2 constraints: <n_1> = 0 and <p_1> = 0.
%  2 = Nx constraints: <f>=0 at each x.

%includeXDotTerm = true;
includeXDotTerm = false;

%includeElectricFieldTermInXiDot = true;
includeElectricFieldTermInXiDot = false;

% If useDKESExBDrift=true, the ExB drift term in the df/dtheta and df/dzeta terms is taken
% to be E x B / <B^2> instead of E x B / B^2.
useDKESExBDrift = true;
%useDKESExBDrift = false;

% --------------------------------------------------
% Numerical resolution parameters:
% --------------------------------------------------

% For each of the quantities below, the 'Converged' value is used except
% when that quantity is being varied in a convergence scan, in which case
% each value in the array that follows (e.g. Nthetas, NLs, etc.) is used.

% Number of grid points in the poloidal direction.
% Memory and time requirements DO depend strongly on this parameter.
NthetaConverged = 30;
Nthetas = floor(linspace(25,35,3));

% Number of grid points in the toroidal direction
% (per identical segment of the stellarator.)
% Memory and time requirements DO depend strongly on this parameter.
NzetaConverged = 5;%25
Nzetas = floor(linspace(3,5,1));

% Number of Legendre polynomials used to represent the distribution
% function.
% Memory and time requirements DO depend strongly on this parameter.
% The value of this parameter required for convergence depends strongly on
% the collisionality. At high collisionality, this parameter can be as low
% as ~ 5. At low collisionality, this parameter may need to be many 10s or
% even > 100 for convergence.
NxiConverged = 80;

Nxis = floor(linspace(40,80,3));

% Number of Legendre polynomials used to represent the Rosenbluth
% potentials: (Typically 2 or 4 is plenty.)
% Memory and time requirements do NOT depend strongly on this parameter.
NLConverged = 4;
NLs = 2:2:6;

% Number of grid points in energy used to represent the distribution
% function.
% Memory and time requirements DO depend strongly on this parameter.
% This parameter almost always needs to be at least 5.
% Usually a value in the range 5-8 is plenty for convergence.
NxConverged = 10;
Nxs=8:2:12;

% Number of grid points in energy used to represent the Rosenbluth
% potentials.
% Memory and time requirements do NOT depend strongly on this parameter.
NxPotentialsPerVthConverged = 41;
%NxPotentialsPerVths = [40, 81];
NxPotentialsPerVths = floor(linspace(20,80,3));

% Tolerance used to define convergence of the Krylov solver.
% This parameter does not affect memory requirements but it does affect the
% time required for solution.
log10tolConverged = 5.5;
log10tols = 4.5:1:6.5;


% --------------------------------------------------
% Other numerical parameters:
% --------------------------------------------------

% For most production runs, you do want to use the Krylov solver rather
% than a direct solver, in which case this parameter should be "true".
% However, for low-resolution problems, or if the Krylov solvers are
% failing to converge, you might want to set this parameter to "false".
tryIterativeSolver = true;
%tryIterativeSolver = false;

%orderOfSolversToTry = [5, 2, 4, 3];
orderOfSolversToTry = [4, 5, 3];
% 1 = BiCGStab
% 2 = BiCGStab(l)
% 3 = CGS
% 4 = GMRES
% 5 = TFQMR

% Below are some setting for the Krylov solvers.
% Maximum number of iterations in the Krylov solver:
maxit = 200;
% You typically do not need to adjust "restart".
restart = maxit; % Used only for GMRES.

thetaGridMode = 2;
% 0 = uniform periodic spectral
% 1 = 2nd order uniform finite-difference
% 2 = 4th order uniform finite-difference
% 3 = 6th order uniform finite-difference
% This parameter should almost always be 2.

forceThetaParity = 1;
% 0 = either even or odd Ntheta is fine.
% 1 = force Ntheta to be odd.
% 2 = force Ntheta to be even.
% This parameter should almost always be 1.

% --------------------------------------------------
% Settings for the preconditioner:
% --------------------------------------------------

preconditioner_x = 1;
% 0 = keep full x coupling.
% 1 = keep only diagonal in x.
% 2 = keep upper-triangular part of x.
% 3 = Keep tridiagonal terms in x.
% 4 = Keep diagonal and superdiagonal in x.
% This parameter should almost always be 1.

preconditioner_x_min_L = 0;
% The simplified x coupling is used in the preconditioner only when the
% Legendre index L is >= this value; otherwise the full x coupling is used
% in the preconditioner.  Set to 0 to precondition at every L.
% Usually, good values for this parameter are 0, 1, or 2.

preconditioner_xi = 0;
% 0 = Keep full xi coupling
% 1 = keep only tridiagonal terms in xi.
% Either 0 or 1 may be appropriate for this parameter.

preconditioner_xi_max_L = inf;
% All L coupling is dropped for L >= this value.
% Recommended value: Inf

preconditioner_theta_min_L = 1;
% The full d/dtheta matrix is used for L < this value.
% Set this to 0 if you don't want to use the full d/dtheta matrix in the
% preconditioner for any L.
% Recommended values: 0, 1, 2, or Inf

preconditioner_theta_max_L = Inf; %4
% All theta coupling is dropped for L >= this value.
% Recommended value: Inf

%preconditioner_theta_remove_cyclic = true;
preconditioner_theta_remove_cyclic = false;
% If true, the (1,end) and (end,1) elements of the d/dtheta matrix are
% removed in the preconditioner.
% Recommended value: false

preconditioner_zeta_min_L = 1;
% The full d/dzeta matrix is used for L < this value.
% Set this to 0 if you don't want to use the full d/dzeta matrix in the
% preconditioner for any L.
% Recommended values: 0, 1, 2, or Inf

preconditioner_zeta_max_L = Inf;
% All theta coupling is dropped for L >= this value.
% Recommended value: Inf

%preconditioner_zeta_remove_cyclic = true;
preconditioner_zeta_remove_cyclic = false;
% If true, the (1,end) and (end,1) elements of the d/dzeta matrix are
% removed in the preconditioner.
% Recommended value: false

% --------------------------------------------------
% Plotting options:
% --------------------------------------------------

% The following offset is added to all figure numbers. It can sometimes
% be convenient to change this number if you want to save figures rather
% than over-write them when re-running the code.
figureOffset=20;

plotSpeedGrid = true;
%plotSpeedGrid = false;

% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% End of the input parameters.
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------

if constraintScheme < 0
    if collisionOperator == 0
        constraintScheme = 1;
    else
        constraintScheme = 2;
    end
end

if testQuasisymmetryIsomorphism
    if geometryScheme ~= 1
        error('To test the quasisymmetry isomorphism, you should set geometryScheme=1.')
    end
    
    if epsilon_t ~= 0
        error('To test the quasisymmetry isomorphism, you should set epsilon_t=0.')
    end
    
    if RHSMode ~= 1
        error('testQuasisymmetryIsomorphism and RHSMode=2 are not presently compatible.')
    end
    
    nuStarS = nuN;
    nuN = nuStarS*abs(helicity_n/iota - helicity_l);
end

FSADensityPerturbation = 0;
FSAFlow = 0;
FSAPressurePerturbation = 0;
particleFlux = 0;
momentumFlux = 0;
heatFlux = 0;
transportMatrix = zeros(3);
dPsidr = NaN;

speedGridFigureHandle = 0;
KrylovFigureHandle = 0;

iteration=0;

if plotSpeedGrid
    figure(figureOffset+7)
    clf
end

fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
fprintf('SFINCS: The Stellarator Fokker-Planck Iterative Neoclassical Conservative Solver.\n')

switch programMode
    case 1
        
        fprintf('Beginning a single run.\n')
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        Ntheta=NthetaConverged;
        Nzeta=NzetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        solveDKE();

    case 2
        fprintf('Beginning convergence scans.\n')
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        startTime=clock;
        
        switch RHSMode
            case 1
                quantitiesToRecord = {'FSAFlow','particleFlux','momentumFlux','heatFlux'};
            case 2
                quantitiesToRecord = {'L11','L12=L21','L13=L31','L12=L21','L22','L23=L32','L13=L31','L23=L32','L33'};
            otherwise
                error('Invalid RHSMode')
        end
            
        linespecs = {'.-b','.-r','.-g','.:c','.-m','.-r','.:k','.:b','.-m'};
        
        parametersToVary = {'N\theta','N\zeta','NL','N\xi','Nx','NxPotentialsPerVth','-log_{10}tol'};
        abscissae = {Nthetas, Nzetas, NLs, Nxis, Nxs, NxPotentialsPerVths, log10tols};
        convergeds = {NthetaConverged, NzetaConverged, NLConverged, NxiConverged, NxConverged, NxPotentialsPerVthConverged, log10tolConverged};
        
        numQuantities = numel(quantitiesToRecord);
        numParameters = numel(parametersToVary);
        quantities = cell(numParameters,1);
        quantities{1} = zeros(numel(Nthetas), numQuantities);
        quantities{2} = zeros(numel(Nzetas), numQuantities);
        quantities{3} = zeros(numel(NLs), numQuantities);
        quantities{4} = zeros(numel(Nxis), numQuantities);
        quantities{5} = zeros(numel(Nxs), numQuantities);
        quantities{6} = zeros(numel(NxPotentialsPerVths), numQuantities);
        quantities{7} = zeros(numel(log10tols), numQuantities);
        parameterScanNum = 1;
        
        % Vary Ntheta, keeping other numerical parameters fixed.
        Nzeta = NzetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(Nthetas)
            Ntheta=Nthetas(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1)=FSAFlow;
                    quantities{parameterScanNum}(iii,2)=particleFlux;
                    quantities{parameterScanNum}(iii,3)=momentumFlux;
                    quantities{parameterScanNum}(iii,4)=heatFlux;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        % Vary Nzeta, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(Nzetas)
            Nzeta=Nzetas(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1)=FSAFlow;
                    quantities{parameterScanNum}(iii,2)=particleFlux;
                    quantities{parameterScanNum}(iii,3)=momentumFlux;
                    quantities{parameterScanNum}(iii,4)=heatFlux;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        % Vary NL, keeping other numerical parameters fixed.
        Nzeta = NzetaConverged;
        Ntheta=NthetaConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(NLs)
            NL=NLs(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1)=FSAFlow;
                    quantities{parameterScanNum}(iii,2)=particleFlux;
                    quantities{parameterScanNum}(iii,3)=momentumFlux;
                    quantities{parameterScanNum}(iii,4)=heatFlux;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        % Vary Nxi, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        Nzeta = NzetaConverged;
        NL=NLConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(Nxis)
            Nxi=Nxis(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1)=FSAFlow;
                    quantities{parameterScanNum}(iii,2)=particleFlux;
                    quantities{parameterScanNum}(iii,3)=momentumFlux;
                    quantities{parameterScanNum}(iii,4)=heatFlux;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        
        % Vary Nx, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        Nzeta = NzetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(Nxs)
            Nx = Nxs(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1)=FSAFlow;
                    quantities{parameterScanNum}(iii,2)=particleFlux;
                    quantities{parameterScanNum}(iii,3)=momentumFlux;
                    quantities{parameterScanNum}(iii,4)=heatFlux;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        % Vary NxPotentialsPerVth, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        Nzeta = NzetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(NxPotentialsPerVths)
            NxPotentialsPerVth = NxPotentialsPerVths(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1)=FSAFlow;
                    quantities{parameterScanNum}(iii,2)=particleFlux;
                    quantities{parameterScanNum}(iii,3)=momentumFlux;
                    quantities{parameterScanNum}(iii,4)=heatFlux;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        % Vary tol, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        Nzeta = NzetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        for iii = 1:numel(log10tols)
            tol = 10^(-log10tols(iii));
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1)=FSAFlow;
                    quantities{parameterScanNum}(iii,2)=particleFlux;
                    quantities{parameterScanNum}(iii,3)=momentumFlux;
                    quantities{parameterScanNum}(iii,4)=heatFlux;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        maxs=ones(numQuantities,1)*(-1e10);
        mins=ones(numQuantities,1)*(1e10);
        for iParameter = 1:numParameters
            maxs = max([maxs, quantities{iParameter}'],[],2);
            mins = min([mins, quantities{iParameter}'],[],2);
        end
        
        temp=dbstack;
        nameOfThisProgram=sprintf('%s',temp.file);
        filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_convergenceScan_',filenameNote];
        outputFilename=[filenameBase,'.mat'];
        if saveStuff
            save(outputFilename)
        end

        plotConvergenceScan()

    case 3
        % Plot a previous convergence scan:
        
        load(dataFileToPlot)
        plotConvergenceScan()

    case 4
        % Do a nuPrime scan:
        
        if RHSMode ~= 2
            error('programMode=4 requires RHSMode=2')
        end
        
        fprintf('Beginning nuPrime scans.\n')
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        %numNus = 13;
        %nuPrimes = logspace(-1,2,numNus);
        
        numNus = 17;
        nuPrimes = logspace(-2,2,numNus);

        
        referenceNuPrimes = [0.01, 0.1, 0.3,  1, 10, 100];
        referenceNthetas =  [  15,  15,  15, 15, 15,  15];
        referenceNzetas =   [  15,  13,  13, 13, 13,  13];
        referenceNxis =     [  48,  37,  34, 13, 13,  13];
        referenceNxs =      [   5,   5,   6,  8,  8,   8];
        
        %{
        referenceNuPrimes = [0.01, 0.1];
        referenceNthetas =  [  15,  15];
        referenceNzetas =   [  15,  13];
        referenceNxis =     [  48,  37];
        referenceNxs =      [   5,   5];
        %}
        %{
        referenceNuPrimes = [0.1, 0.3,  1, 10, 100];
        referenceNthetas =  [ 19,  19, 15, 15,  15];
        referenceNzetas =   [ 13,  13, 13, 13,  13];
        referenceNxis =     [ 37,  34, 13, 13,  13];
        referenceNxs =      [  7, 7.5,  8,  8,   8];
        %}
        %{
        referenceNuPrimes = [0.1,  1, 10, 100];
        referenceNthetas =  [ 13,  9,  7,   7];
        referenceNzetas =   [  9,  9,  7,   7];
        referenceNxis =     [ 11,  8,  6,   6];
        referenceNxs =      [  4,  4,  5,   5];
        %}
        
        numReferenceNus = numel(referenceNuPrimes);
        if numel(referenceNthetas) ~= numReferenceNus || numel(referenceNzetas) ~= numReferenceNus || numel(referenceNxis) ~= numReferenceNus || numel(referenceNxs) ~= numReferenceNus
            error('Number of reference points for resolution is not consistent')
        end
        
        NthetaMultipliers = [1, 2, 1, 1, 1];
        NzetaMultipliers =  [1, 1, 2, 1, 1];
        NxiMultipliers    = [1, 1, 1, 2, 1];
        NxMultipliers     = [1, 1, 1, 1, 1.7];
    
        NConvergence = numel(NthetaMultipliers);
        if numel(NzetaMultipliers) ~= NConvergence || numel(NxiMultipliers) ~= NConvergence || numel(NxMultipliers) ~= NConvergence
            error('Sizes of multiplier arrays are not consistent')
        end
        
        numRuns = 3*NConvergence*numNus;
        scanResults = zeros(3, NConvergence, numNus, 3, 3);
        
        runNum=0;
        
        Nthetas = interp1(log10(referenceNuPrimes), referenceNthetas, log10(nuPrimes),'cubic');
        Nzetas = interp1(log10(referenceNuPrimes), referenceNzetas, log10(nuPrimes),'cubic');
        Nxis = interp1(log10(referenceNuPrimes), referenceNxis, log10(nuPrimes),'cubic');
        Nxs = interp1(log10(referenceNuPrimes), referenceNxs, log10(nuPrimes),'cubic');
        
        NL=NLConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);

        for iConvergence = 1:NConvergence
            for iCollisionOperator = 1:3
                for iNu = 1:numNus
                    runNum = runNum + 1;
                    fprintf('##################################################################\n')
                    fprintf('Run %d of %d.\n', runNum, numRuns)
                    fprintf('##################################################################\n')
                    
                    nuPrime = nuPrimes(iNu);
                    if nuPrime > 1
                        preconditioner_x_min_L = 1;
                    else
                        preconditioner_x_min_L = 0;
                    end
                    
                    collisionOperator = iCollisionOperator-1;
                    Ntheta=floor(Nthetas(iNu)*NthetaMultipliers(iConvergence));
                    Nzeta=floor(Nzetas(iNu)*NzetaMultipliers(iConvergence));
                    Nxi=floor(Nxis(iNu)*NxiMultipliers(iConvergence));
                    Nx=floor(Nxs(iNu)*NxMultipliers(iConvergence));
                    
                    solveDKE()
                    
                    scanResults(iCollisionOperator, iConvergence, iNu, :, :) = transportMatrix;
                end
            end
        end
        
        temp=dbstack;
        nameOfThisProgram=sprintf('%s',temp.file);
        filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_nuPrimeScan_',filenameNote];
        outputFilename=[filenameBase,'.mat'];
        if saveStuff
            save(outputFilename)
        end

        plotNuPrimeScan()
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Done with nuPrime scans.\n')
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case 5
        % Plot results of a previous nuPrime scan:
        
        load(dataFileToPlot)
        plotNuPrimeScan()
        
    otherwise
        error('Invalid setting for programMode.')
end

    function plotConvergenceScan()
        % ------------------------------------------------------
        % Plot results of the convergence scan:
        % ------------------------------------------------------
        
        quantityToRowMap = [1, 2, 3, 2, 4, 5, 3, 5, 6];
        switch RHSMode
            case 1
                numRows = numQuantities;
            case 2
                numRows = 6;
        end
        numCols = numParameters;
        
        figure(1+figureOffset)
        clf
        for iQuantity = 1:numQuantities
            if maxs(iQuantity) <= mins(iQuantity)
                maxs(iQuantity) = mins(iQuantity)+1;
            end
            switch RHSMode
                case 1
                    iRow = iQuantity;
                case 2
                    iRow = quantityToRowMap(iQuantity);
            end
            for iParameter = 1:numParameters
                subplot(numRows, numCols, iParameter  + (iRow - 1)*numParameters)
                plot(1./abscissae{iParameter}, quantities{iParameter}(:,iQuantity)', linespecs{iQuantity})
                hold on
                plot(1./[convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
                ylim([mins(iQuantity), maxs(iQuantity)])
                xlabel(['1/',parametersToVary{iParameter}])
                ylabel(quantitiesToRecord{iQuantity})
            end
        end
        switch RHSMode
            case 1
                stringForTop = sprintf('SFINCS convergence scan: nuN=%g. Base case: Ntheta=%d, Nzeta=%d, NL=%d, Nxi=%d, Nx=%d, NxPotentialsPerVth=%g, -log10tol=%g.', ...
                    nuN, NthetaConverged, NzetaConverged, NLConverged, NxiConverged, NxConverged, NxPotentialsPerVthConverged, log10tolConverged);
            case 2
                stringForTop = sprintf('SFINCS convergence scan: nuPrime=%g. Base case: Ntheta=%d, Nzeta=%d, NL=%d, Nxi=%d, Nx=%d, NxPotentialsPerVth=%g, -log10tol=%g.', ...
                    nuPrime, NthetaConverged, NzetaConverged, NLConverged, NxiConverged, NxConverged, NxPotentialsPerVthConverged, log10tolConverged);
        end
        annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
            'Interpreter','none','VerticalAlignment','bottom',...
            'FontSize',12,'LineStyle','none','String',stringForTop);
        
        figure(8+figureOffset)
        clf
        for iQuantity = 1:numQuantities
            if maxs(iQuantity) <= mins(iQuantity)
                maxs(iQuantity) = mins(iQuantity)+1;
            end
            switch RHSMode
                case 1
                    iRow = iQuantity;
                case 2
                    iRow = quantityToRowMap(iQuantity);
            end
            for iParameter = 1:numParameters
                subplot(numRows, numCols, iParameter  + (iRow - 1)*numParameters)
                plot(abscissae{iParameter}, quantities{iParameter}(:,iQuantity)', linespecs{iQuantity})
                hold on
                plot([convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
                ylim([mins(iQuantity), maxs(iQuantity)])
                xlabel(parametersToVary{iParameter})
                ylabel(quantitiesToRecord{iQuantity})
            end
        end
        
        annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
            'Interpreter','none','VerticalAlignment','bottom',...
            'FontSize',12,'LineStyle','none','String',stringForTop);
        
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time for convergence scans: %g seconds.\n',etime(clock,startTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    end

    function plotNuPrimeScan()
        legendText = {'base case','2x N\theta','2x N\zeta','2x N\xi','2x Nx'};
        colors = [  1,0,0;
                    0.7,0.5,0;
                    0,0.7,0;
                    0,0,1;
                    1,0,1];
        
        for iCollision = 1:3
            figure(iCollision+figureOffset)
            clf
            numRows = 2;
            numCols = 3;
            plotNum = 1;
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 1, 1))), '.-','Color',colors(iConvergence,:))
                hold on
            end
            xlabel('nuPrime')
            ylabel('-L11')
            legend(legendText)
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 1, 2))), '.-','Color',colors(iConvergence,:))
                hold on
                loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 2, 1))), '.:','Color',colors(iConvergence,:))
            end
            xlabel('nuPrime')
            ylabel('L12=L21')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                semilogx(nuPrimes, squeeze(scanResults(iCollision,iConvergence, :, 1, 3)), '.-','Color',colors(iConvergence,:))
                hold on
                semilogx(nuPrimes, squeeze(scanResults(iCollision,iConvergence, :, 3, 1)), '.:','Color',colors(iConvergence,:))
            end
            xlabel('nuPrime')
            ylabel('L13=L31')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 2, 2))), '.-','Color',colors(iConvergence,:))
                hold on
            end
            xlabel('nuPrime')
            ylabel('-L22')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                semilogx(nuPrimes, squeeze(scanResults(iCollision,iConvergence, :, 2, 3)), '.-','Color',colors(iConvergence,:))
                hold on
                semilogx(nuPrimes, squeeze(scanResults(iCollision,iConvergence, :, 3, 2)), '.:','Color',colors(iConvergence,:))
            end
            xlabel('nuPrime')
            ylabel('L23=L32')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 3, 3))), '.-','Color',colors(iConvergence,:))
                hold on
            end
            xlabel('nuPrime')
            ylabel('L33')
        end
        
        legendText = {'Fokker-Planck','pitch-angle scattering','momentum-conserving model'};
        colors = [  1,0,0;
                    0,0.7,0;
                    0,0,1];

        figure(4+figureOffset)
        clf
        numRows = 2;
        numCols = 3;
        plotNum = 1;
        iConvergence = 1;
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 1, 1))), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('-L11')
        legend(legendText)
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            data = 0.5*(scanResults(iCollision,iConvergence, :, 1, 2) + scanResults(iCollision,iConvergence, :, 2, 1));
            loglog(nuPrimes, abs(squeeze(data)), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('L12=L21')
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            data = 0.5*(scanResults(iCollision,iConvergence, :, 1, 3) + scanResults(iCollision,iConvergence, :, 3, 1));
            %loglog(nuPrimes, abs(squeeze(data)), '.-','Color',colors(iCollision,:))
            semilogx(nuPrimes, squeeze(data), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('L13=L31')
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 2, 2))), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('-L22')
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            data = 0.5*(scanResults(iCollision,iConvergence, :, 3, 2) + scanResults(iCollision,iConvergence, :, 2, 3));
            %loglog(nuPrimes, abs(squeeze(data)), '.-','Color',colors(iCollision,:))
            semilogx(nuPrimes, squeeze(data), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('L23=L32')
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 3, 3))), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('L33')
    end

% --------------------------------------------------------
% --------------------------------------------------------
% Done with the routines for convergence scans.
% Next comes the core function of the code.
% --------------------------------------------------------
% --------------------------------------------------------

    function solveDKE()
        
        startTimeForThisRun=tic;
        
        sqrtpi=sqrt(pi);
        iteration = iteration+1;
        
        
        % --------------------------------------------------------
        % --------------------------------------------------------
        % First, set up the grids, differentiation matrices, and
        % integration weights for each coordinate.
        % --------------------------------------------------------
        % --------------------------------------------------------
        
        
        switch forceThetaParity
            case 0
                % Do nothing
            case 1
                % For Ntheta to be odd
                if mod(Ntheta,2)==0
                    Ntheta=Ntheta+1;
                end
            case 2
                % For Ntheta to be even
                if mod(Ntheta,2)==1
                    Ntheta=Ntheta+1;
                end
            otherwise
                error('Invalid forceThetaParity')
        end
        
        if mod(Nzeta,2)==0
            Nzeta=Nzeta+1;
        end
        
        if iteration>1
            fprintf('********************************************************************\n')
        end
        fprintf('Ntheta = %d,  Nzeta = %d,  NL = %d,  Nxi = %d,  Nx = %d, NxPtentialsPerVth = %g, tol = %g\n',Ntheta, Nzeta,NL,Nxi,Nx,NxPotentialsPerVth,tol)
        
        tic
        
        % Generate abscissae, quadrature weights, and derivative matrix for theta grid.
        if Ntheta == 1
            theta = 0;
            thetaWeights = 2*pi;
            ddtheta = 0;
            ddtheta_preconditioner = 0;
        else
            switch thetaGridMode
                case 0
                    % Spectral uniform
                    scheme = 20;
                case 1
                    % Uniform periodic 2nd order FD
                    scheme = 0;
                case 2
                    % Uniform periodic 4th order FD
                    scheme = 10;
                case 3
                    % Uniform periodic 6th order FD
                    scheme = 70;
                otherwise
                    error('Error! Invalid thetaGridMode')
            end
            [theta, thetaWeights, ddtheta, ~] = m20121125_04_DifferentiationMatricesForUniformGrid(Ntheta, 0, 2*pi, scheme);
            
            scheme = 0;
            [~, ~, ddtheta_preconditioner, ~] = m20121125_04_DifferentiationMatricesForUniformGrid(Ntheta, 0, 2*pi, scheme);
            if preconditioner_theta_remove_cyclic
                ddtheta_preconditioner(1,end) = 0;
                ddtheta_preconditioner(end,1) = 0;
            end
           
        end
        
        
        % Generate abscissae, quadrature weights, and derivative matrix for zeta grid.
        setNPeriods()
        zetaMax = 2*pi/NPeriods;
        
        if Nzeta==1
            zeta=0;
            zetaWeights=2*pi;
            ddzeta=0;
            ddzeta_preconditioner=0;
        else
            switch thetaGridMode
                case 0
                    % Spectral uniform
                    scheme = 20;
                case 1
                    % Uniform periodic 2nd order FD
                    scheme = 0;
                case 2
                    % Uniform periodic 4th order FD
                    scheme = 10;
                case 3
                    % Uniform periodic 6th order FD
                    scheme = 70;
                otherwise
                    error('Error! Invalid thetaGridMode')
            end
            [zeta, zetaWeights, ddzeta, ~] = m20121125_04_DifferentiationMatricesForUniformGrid(Nzeta, 0, zetaMax, scheme);
            zetaWeights = zetaWeights * NPeriods;
            
            scheme = 0;
            [~, ~, ddzeta_preconditioner, ~] = m20121125_04_DifferentiationMatricesForUniformGrid(Nzeta, 0, zetaMax, scheme);
            if preconditioner_zeta_remove_cyclic
                ddzeta_preconditioner(1,end) = 0;
                ddzeta_preconditioner(end,1) = 0;                
            end
        end
        
        % Evaluate the magnetic field and its derivatives on the
        % (theta,zeta) grid:
        computeBHat()
        
        % Compute a few quantities related to the magnetic field
        VPrimeHat = thetaWeights' * (1./BHat.^2) * zetaWeights;
        FSABHat2 = 4*pi*pi/VPrimeHat;
        
        % Generate abscissae, quadrature weights, and derivative matrices for
        % the energy (x) grid used to represent the distribution function.
        k=0;
        scale=1;
        pointAtZero = false;
        [x, ddx, d2dx2, xWeights] = m20130312_02_SpectralNodesWeightsAndDifferentiationMatricesForV(Nx, k, scale, pointAtZero);
        
        % Make the energy grid and differentiation matrices for the
        % Rosenbluth potentials:
        function y=weight(x)
            x2=x.*x;
            y=exp(-x2);
        end
        xMax=max([5, max(x)]);
        xMin=0;
        NxPotentials = ceil(xMax * NxPotentialsPerVth);
        % Uniform grid with 5-point stencil for derivatives:
        scheme = 12;
        [xPotentials, ~, ddxPotentials, d2dx2Potentials] = m20121125_04_DifferentiationMatricesForUniformGrid(NxPotentials, xMin, xMax, scheme);
        
        % Make the matrices for interpolating between the two energy grids:
        regridPolynomialToUniform = m20120703_03_polynomialInterpolationMatrix(x,xPotentials,weight(x),weight(xPotentials));
        regridUniformToPolynomial = m20121127_02_makeHighOrderInterpolationMatrix(xPotentials,x,0,'f');
        
        
        if plotSpeedGrid
            if iteration == 1
                speedGridFigureHandle = figure(figureOffset+7);
            else
                set(0, 'CurrentFigure', speedGridFigureHandle);
            end
            plot(xPotentials,zeros(size(xPotentials))+iteration,'.r')
            hold on
            plot(x, zeros(size(x))+iteration,'o')
            title('Speed grid for distribution function (blue) and Rosenbluth potentials(red)')
            xlabel('x')
            ylabel('iteration')
        end
        
        % Set the size of the main linear system:
        matrixSize = Nx * Nxi * Ntheta * Nzeta;
        switch constraintScheme
            case 0
                % Nothing to do here.
            case 1
                matrixSize = matrixSize + 2;
            case 2
                matrixSize = matrixSize + Nx;
            otherwise
                error('Invalid setting for constraintScheme')
        end
        
        % To build the matrix as efficiently as possible, a reasonably
        % accurate estimate of the number of nonzeros (nnz) is needed beforehand:
        estimated_nnz = 1*(Nx*3*Nxi*nnz(ddtheta)*Nzeta + Nx*Nxi*3*Ntheta*nnz(ddzeta) + Nx*5*Nxi*Ntheta*Nzeta + Nx*Nx*4*Nxi*Ntheta*Nzeta + Nx*Nxi*Ntheta*Nzeta*4);
        estimated_nnz_original = estimated_nnz;
        fprintf('matrixSize: %d.\n',matrixSize)
        
        if iteration==1
            figure(figureOffset+4);
            clf
            numRows=3;
            numCols=3;
            plotNum=1;
            numContours=15;
            
            subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
            contourf(zeta,theta,BHat,numContours,'EdgeColor','none')
            colorbar
            xlabel('\zeta')
            ylabel('\theta')
            title('BHat')
            
            subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
            contourf(zeta,theta,dBHatdtheta,numContours,'EdgeColor','none')
            colorbar
            xlabel('\zeta')
            ylabel('\theta')
            title('dBHatdtheta')
            
            subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
            contourf(zeta,theta,dBHatdzeta,numContours,'EdgeColor','none')
            colorbar
            xlabel('\zeta')
            ylabel('\theta')
            title('dBHatdzeta')
            
            drawnow
        end
        
        if RHSMode == 2
            nuN = nuPrime * sqrt(THat) * B0OverBBar / (GHat+iota*IHat);
        end
        
        % Begin timer for matrix construction:
        tic
        
        % Order of the rows of the matrix and of the RHS:
        % --------------------------------
        % for i = 1:Nx
        %   for L = 0:(Nxi-1)
        %     for j = 1:Ntheta
        %       for k = 1:Nzeta
        %         Enforce the drift-kinetic equation.
        % Enforce any constraints.
        
        
        % Order of the vector of unknowns & of columns in the matrix:
        % --------------------------------
        % for i = 1:Nx
        %   for L = 0:(Nxi-1)
        %     for j = 1:Ntheta
        %       for k = 1:Nzeta
        %         Value of the distribution function.
        % Value of any sources.
        
        
        % ------------------------------------------------------
        % ------------------------------------------------------
        % Build the right-hand side of the main linear system.
        % ------------------------------------------------------
        % ------------------------------------------------------
        
        switch RHSMode
            case 1
                RHSSize = 1;
            case 2
                RHSSize = 3;
            otherwise
                error('Invalid RHSMode')
        end
        
        x2=x.*x;
        expx2=exp(-x2);
        sqrtTHat = sqrt(THat);
        rhs=zeros(matrixSize,RHSSize);
        
        spatialPartOfRHS_gradients = (GHat*dBHatdtheta - IHat*dBHatdzeta) ./ (2*(BHat.^3)*sqrtTHat);
        spatialPartOfRHS_EPar = 2*omega*psiAHat*(GHat+iota*IHat)./(Delta*Delta*THat*THat*FSABHat2*BHat);

        for col=1:RHSSize
            switch RHSMode
                case 1
                    dnHatdpsiToUse = dnHatdpsi;
                    dTHatdpsiToUse = dTHatdpsi;
                    dPhiHatdpsiToUse = dPhiHatdpsi;
                    EHatToUse = EHat;
                case 2
                    dPhiHatdpsiToUse = 0;
                    switch col
                        case 1
                            dnHatdpsiToUse = 1;
                            dTHatdpsiToUse = 0;
                            EHatToUse = 0;
                        case 2
                            % The next 2 lines ensure (1/n)*dn/dpsi + (3/2)*dT/dpsi = 0 while dT/dpsi is nonzero.
                            dnHatdpsiToUse = (3/2)*nHat/THat;
                            dTHatdpsiToUse = 1;
                            EHatToUse = 0;
                        case 3
                            dnHatdpsiToUse = 0;
                            dTHatdpsiToUse = 0;
                            EHatToUse = 1;
                    end
                otherwise
                    error('Invalid RHSMode')
            end
            
            xPartOfRHS_gradients = x2.*expx2.*(dnHatdpsiToUse/nHat + 2*omega*dPhiHatdpsiToUse/(Delta*THat) + (x2-3/2)*dTHatdpsiToUse/THat);
            xPartOfRHS_EPar = EHatToUse*x.*expx2;
                        
            for ix=1:Nx
                for itheta=1:Ntheta
                    L=0;
                    indices = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + (1:Nzeta);
                    rhs(indices, col) = (4/3) * spatialPartOfRHS_gradients(itheta,:)' * xPartOfRHS_gradients(ix);
                    
                    L=1;
                    indices = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + (1:Nzeta);
                    rhs(indices, col) = spatialPartOfRHS_EPar(itheta,:)' * xPartOfRHS_EPar(ix);
                    
                    L=2;
                    indices = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + (1:Nzeta);
                    rhs(indices, col) = (2/3) * spatialPartOfRHS_gradients(itheta,:)' * xPartOfRHS_gradients(ix);
                end
            end
        end
        
        
        assignin('base','rhsm',rhs)
        
        % ------------------------------------------------------
        % ------------------------------------------------------
        % Build the matrix for the main linear system.
        % ------------------------------------------------------
        % ------------------------------------------------------
        
        sparseCreatorIndex=1;
        sparseCreator_i=0;
        sparseCreator_j=0;
        sparseCreator_s=0;
        resetSparseCreator()
        
        if tryIterativeSolver
            matricesToMake=0:1;
        else
            matricesToMake=1;
        end
        
        for whichMatrixToMake = matricesToMake
            % 0 = preconditioner
            % 1 = main matrix
            
            if whichMatrixToMake==1
                ddthetaToUse = ddtheta;
                ddzetaToUse = ddzeta;
                maxLForThetaDot = Nxi-1;
                maxLForZetaDot = Nxi-1;
                maxLForXiDot = Nxi-1;
            else
                ddthetaToUse = ddtheta_preconditioner;
                ddzetaToUse = ddzeta_preconditioner;
                maxLForThetaDot = min([preconditioner_theta_max_L, Nxi-1]);
                maxLForZetaDot = min([preconditioner_zeta_max_L, Nxi-1]);
                maxLForXiDot = min([preconditioner_xi_max_L, Nxi-1]);
            end
            
            matrixStartTime = tic;
            
            % -----------------------------------------
            % Add d/dtheta terms:
            % -----------------------------------------
            
            for izeta=1:Nzeta
                if useDKESExBDrift
                    thetaPartOfExBTerm_lowL = omega*GHat*dPhiHatdpsi/(psiAHat*FSABHat2) * ddtheta;
                    thetaPartOfExBTerm_highL = omega*GHat*dPhiHatdpsi/(psiAHat*FSABHat2) * ddthetaToUse;
                else
                    thetaPartOfExBTerm_lowL = omega*GHat*dPhiHatdpsi/psiAHat * diag(1./BHat(:,izeta).^2)*ddtheta;
                    thetaPartOfExBTerm_highL = omega*GHat*dPhiHatdpsi/psiAHat * diag(1./BHat(:,izeta).^2)*ddthetaToUse;
                end
                thetaPartOfStreamingTerm_lowL = iota*sqrtTHat*diag(1./BHat(:,izeta))*ddtheta;
                thetaPartOfStreamingTerm_highL = iota*sqrtTHat*diag(1./BHat(:,izeta))*ddthetaToUse;
                for L=0:maxLForThetaDot
                    if L < preconditioner_theta_min_L
                        thetaPartOfStreamingTerm = thetaPartOfStreamingTerm_lowL;
                        thetaPartOfExBTerm = thetaPartOfExBTerm_lowL;
                    else
                        thetaPartOfStreamingTerm = thetaPartOfStreamingTerm_highL;
                        thetaPartOfExBTerm = thetaPartOfExBTerm_highL;
                    end
                    
                    for ix=1:Nx
                        rowIndices = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + ((1:Ntheta)-1)*Nzeta + izeta;
                        
                        % Diagonal term
                        addSparseBlock(rowIndices, rowIndices, thetaPartOfExBTerm)
                        
                        % Super-diagonal term
                        if (L<maxLForThetaDot)
                            colIndices = rowIndices + Ntheta*Nzeta;
                            addSparseBlock(rowIndices, colIndices, x(ix)*(L+1)/(2*L+3)*thetaPartOfStreamingTerm)
                        end
                        
                        % Sub-diagonal term
                        if (L>0)
                            colIndices = rowIndices - Ntheta*Nzeta;
                            addSparseBlock(rowIndices, colIndices, x(ix)*L/(2*L-1)*thetaPartOfStreamingTerm)
                        end
                        
                    end
                end
            end
            
            % -----------------------------------------
            % Add d/dzeta terms:
            % -----------------------------------------
            
            for itheta=1:Ntheta
                if useDKESExBDrift
                    zetaPartOfExBTerm_lowL = -omega*IHat*dPhiHatdpsi/(psiAHat*FSABHat2) *ddzeta;
                    zetaPartOfExBTerm_highL = -omega*IHat*dPhiHatdpsi/(psiAHat*FSABHat2) *ddzetaToUse;
                else
                    zetaPartOfExBTerm_lowL = -omega*IHat*dPhiHatdpsi/psiAHat * diag(1./BHat(itheta,:).^2)*ddzeta;
                    zetaPartOfExBTerm_highL = -omega*IHat*dPhiHatdpsi/psiAHat * diag(1./BHat(itheta,:).^2)*ddzetaToUse;
                end
                zetaPartOfStreamingTerm_lowL = sqrtTHat*diag(1./BHat(itheta,:))*ddzeta;
                zetaPartOfStreamingTerm_highL = sqrtTHat*diag(1./BHat(itheta,:))*ddzetaToUse;
                for L=0:maxLForZetaDot
                    if L < preconditioner_zeta_min_L
                        zetaPartOfExBTerm = zetaPartOfExBTerm_lowL;
                        zetaPartOfStreamingTerm = zetaPartOfStreamingTerm_lowL;
                    else
                        zetaPartOfExBTerm = zetaPartOfExBTerm_highL;
                        zetaPartOfStreamingTerm = zetaPartOfStreamingTerm_highL;
                    end
                    
                    for ix=1:Nx
                        rowIndices = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta +(1:Nzeta);
                        
                        % Diagonal term
                        addSparseBlock(rowIndices, rowIndices, zetaPartOfExBTerm)
                        
                        % Super-diagonal term
                        if (L<maxLForZetaDot)
                            colIndices = rowIndices + Ntheta*Nzeta;
                            addSparseBlock(rowIndices, colIndices, x(ix)*(L+1)/(2*L+3)*zetaPartOfStreamingTerm)
                        end
                        
                        % Sub-diagonal term
                        if (L>0)
                            colIndices = rowIndices - Ntheta*Nzeta;
                            addSparseBlock(rowIndices, colIndices, x(ix)*L/(2*L-1)*zetaPartOfStreamingTerm)
                        end
                        
                    end
                end
            end
            
            
            % -----------------------------------------
            % Add d/dxi terms:
            % -----------------------------------------
            
            for itheta=1:Ntheta
                spatialPartOfOldMirrorTerm = -sqrtTHat*(iota*dBHatdtheta(itheta,:)+dBHatdzeta(itheta,:))./(2*BHat(itheta,:).^2);
                spatialPartOfNewMirrorTerm = omega*dPhiHatdpsi*(GHat*dBHatdtheta(itheta,:) - IHat*dBHatdzeta(itheta,:))./(2*psiAHat*BHat(itheta,:).^3);
                for ix=1:Nx
                    for L=0:maxLForXiDot
                        rowIndices = (ix-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + (1:Nzeta);
                        
                        % Super-diagonal term
                        if (L<maxLForXiDot)
                            colIndices = rowIndices + Ntheta*Nzeta;
                            addToSparse(rowIndices, colIndices, x(ix)*(L+1)*(L+2)/(2*L+3)*spatialPartOfOldMirrorTerm)
                        end
                        
                        % Sub-diagonal term
                        if (L>0)
                            colIndices = rowIndices - Ntheta*Nzeta;
                            addToSparse(rowIndices, colIndices, x(ix)*(-L)*(L-1)/(2*L-1)*spatialPartOfOldMirrorTerm)
                        end
                        
                        if includeElectricFieldTermInXiDot
                            % Diagonal term
                            addToSparse(rowIndices, rowIndices, L*(L+1)/((2*L-1)*(2*L+3))*spatialPartOfNewMirrorTerm)
                            
                            if (whichMatrixToMake==1 || preconditioner_xi==0)
                                % Super-super-diagonal term:
                                if (L < maxLForXiDot-1)
                                    colIndices = rowIndices + 2*Ntheta*Nzeta;
                                    addToSparse(rowIndices, colIndices, (L+1)*(L+2)*(L+3)/((2*L+5)*(2*L+3))*spatialPartOfNewMirrorTerm)
                                end
                                
                                % Sub-sub-diagonal term:
                                if (L > 1)
                                    colIndices = rowIndices - 2*Ntheta*Nzeta;
                                    addToSparse(rowIndices, colIndices, -L*(L-1)*(L-2)/((2*L-3)*(2*L-1))*spatialPartOfNewMirrorTerm)
                                end
                            end
                        end
                    end
                end
            end
            
            
            % -----------------------------------------
            % Add the collisionless d/dx term:
            % -----------------------------------------
            
            if includeXDotTerm
                xPartOfXDot = diag(x)*ddx;
                if (whichMatrixToMake==1)
                    xPartOfXDotForLargeL = xPartOfXDot;
                else
                    % We're making the preconditioner, so simplify matrix
                    % if needed:
                    switch preconditioner_x
                        case 0
                            xPartOfXDotForLargeL = xPartOfXDot;
                        case 1
                            xPartOfXDotForLargeL = diag(diag(xPartOfXDot));
                        case 2
                            xPartOfXDotForLargeL = triu(xPartOfXDot);
                        case 3
                            mask = eye(Nx) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1);
                            xPartOfXDotForLargeL = xPartOfXDot .* mask;
                        case 4
                            mask = eye(Nx) + diag(ones(Nx-1,1),1);
                            xPartOfXDotForLargeL = xPartOfXDot .* mask;
                        otherwise
                            error('Invalid setting for preconditioner_x')
                    end
                end
                for L=0:(Nxi-1)
                    if L >= preconditioner_x_min_L
                        xPartOfXDotToUse = xPartOfXDotForLargeL;
                    else
                        xPartOfXDotToUse = xPartOfXDot;
                    end
                    for itheta=1:Ntheta
                        for izeta=1:Nzeta
                            spatialPart = omega*dPhiHatdpsi*(GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))/(2*psiAHat*BHat(itheta,izeta)^3);
                            
                            rowIndices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                            
                            % Diagonal term
                            addSparseBlock(rowIndices, rowIndices, 2*(3*L*L+3*L-2)/((2*L+3)*(2*L-1))*spatialPart*xPartOfXDotToUse)
                            
                            if (whichMatrixToMake==1 || preconditioner_xi==0)
                                % Super-super-diagonal in L
                                if (L<Nxi-2)
                                    colIndices = rowIndices + 2*Ntheta*Nzeta;
                                    addSparseBlock(rowIndices, colIndices, (L+1)*(L+2)/((2*L+5)*(2*L+3))*spatialPart*xPartOfXDotToUse)
                                end
                                
                                % Sub-sub-diagonal in L
                                if (L>1)
                                    colIndices = rowIndices - 2*Ntheta*Nzeta;
                                    addSparseBlock(rowIndices, colIndices, L*(L-1)/((2*L-3)*(2*L-1))*spatialPart*xPartOfXDotToUse)
                                end
                                
                            end
                        end
                    end
                end
            end
            
            % -----------------------------------------
            % Add collision operator:
            % -----------------------------------------
            
            erfs=erf(x);
            x2 = x.*x;
            x3 = x2.*x;
            expx2 = exp(-x.*x);
            % Psi is the Chandrasekhar function:
            Psi = (erfs - 2/sqrtpi*x .* expx2) ./ (2*x.*x);
            % Energy-dependent deflection frequency:
            nuD = 3*sqrtpi/4*(erfs - Psi) ./ x3;
            
            switch collisionOperator
                case 0
                    % Full linearized Fokker-Planck operator
                    
                    xWith0s = [0; xPotentials(2:(end-1)); 0];
                    M21 = 4*pi*diag(xWith0s.^2) * regridPolynomialToUniform;
                    M32 = -2*diag(xWith0s.^2);
                    LaplacianTimesX2WithoutL = diag(xPotentials.^2)*d2dx2Potentials + 2*diag(xPotentials)*ddxPotentials;
                    
                    PsiPrime = (-erfs + 2/sqrtpi*x.*(1+x.*x) .* expx2) ./ x3;
                    xPartOfCECD = 3*sqrtpi/4*(diag(Psi./x)*d2dx2   +  diag((PsiPrime.*x  + Psi + 2*Psi.*x2)./x2)*ddx + diag(2*PsiPrime + 4*Psi./x)) + 3*diag(expx2);
                    M12IncludingX0 = nuN * 3/(2*pi)*diag(expx2)*regridUniformToPolynomial;
                    M13IncludingX0 = -nuN * 3/(2*pi) * diag(x2.*expx2) * regridUniformToPolynomial* d2dx2Potentials;
                    
                    for L=0:(Nxi-1)
                        M11 = -nuN * (-0.5*diag(nuD)*L*(L+1) + xPartOfCECD);
                        
                        if L <= (NL-1)
                            % Add Rosenbluth potential stuff
                            
                            M13 = M13IncludingX0;
                            M12 = M12IncludingX0;
                            
                            M22 = LaplacianTimesX2WithoutL-L*(L+1)*eye(NxPotentials);
                            % Add Dirichlet or Neumann boundary condition for
                            % potentials at x=0:
                            if L==0
                                M22(1,:)=ddxPotentials(1,:);
                            else
                                M22(1,:) = 0;
                                M22(1,1) = 1;
                                M12(:,1) = 0;
                                M13(:,1) = 0;
                            end
                            M33 = M22;
                            
                            % Add Robin boundary condition for potentials at x=xMax:
                            M22(NxPotentials,:) = xMax*ddxPotentials(NxPotentials,:);
                            M22(NxPotentials,NxPotentials) = M22(NxPotentials,NxPotentials) + L+1;
                            
                            M33(NxPotentials,:) = xMax*xMax*d2dx2Potentials(NxPotentials,:) + (2*L+1)*xMax*ddxPotentials(NxPotentials,:);
                            M33(NxPotentials,NxPotentials) = M33(NxPotentials,NxPotentials) + (L*L-1);
                            if L~=0
                                M22(NxPotentials,1)=0;
                                M33(NxPotentials,1)=0;
                            end
                            
                            KWithoutThetaPart = M11 -  (M12 - M13 * (M33 \ M32)) * (M22 \ M21);
                            
                        else
                            KWithoutThetaPart = M11;
                        end
                        
                        if (whichMatrixToMake==0 && L >= preconditioner_x_min_L)
                            % We're making the preconditioner, so simplify
                            % the matrix if needed.
                            switch preconditioner_x
                                case 0
                                    % Do nothing.
                                case 1
                                    KWithoutThetaPart = diag(diag(KWithoutThetaPart));
                                case 2
                                    KWithoutThetaPart = triu(KWithoutThetaPart);
                                case 3
                                    mask = eye(Nx) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1);
                                    KWithoutThetaPart = KWithoutThetaPart .* mask;
                                case 4
                                    mask = eye(Nx) + diag(ones(Nx-1,1),1);
                                    KWithoutThetaPart = KWithoutThetaPart .* mask;
                                otherwise
                                    error('Invalid setting for preconditioner_x')
                            end
                        end
                        
                        for itheta=1:Ntheta
                            for izeta=1:Nzeta
                                indices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                                addSparseBlock(indices, indices, (GHat+iota*IHat)/(BHat(itheta,izeta)^2)*KWithoutThetaPart)
                            end
                        end
                    end
                    
                case {1,2,3}
                    % Pitch angle scattering operator
                    if collisionOperator==3
                      nuD = 3*sqrtpi/4;%nuDconst;%3*sqrtpi/4*(erf(1) - (erf(1) - 2/sqrtpi .* exp(-1)) ./ (2));
                    end
                    
                    for itheta=1:Ntheta
                        for izeta=1:Nzeta
                            spatialPart = -nuN*(GHat+iota*IHat)/(BHat(itheta,izeta)^2);
                            for L=1:(Nxi-1)
                                indices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                                addToSparse(indices, indices, -0.5*L*(L+1)*spatialPart*nuD)
                            end
                        end
                    end
                    
                    if collisionOperator==2 %|| collisionOperator==3
                        % Add model field term
                        
                        L=1;
                        fieldTerm = (nuD.*x.*expx2) * ((xWeights.*x.*x2.*nuD)')/0.354162849836926;
                        
                        if (whichMatrixToMake==0 && L >= preconditioner_x_min_L)
                            switch preconditioner_x
                                case 0
                                    % Nothing to do here.
                                case 1
                                    fieldTerm = diag(diag(fieldTerm));
                                case 2
                                    fieldTerm = triu(fieldTerm);
                                case 3
                                    mask = eye(Nx) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1);
                                    fieldTerm = fieldTerm .* mask;
                                case 4
                                    mask = eye(Nx) + diag(ones(Nx-1,1),1);
                                    fieldTerm = fieldTerm .* mask;
                                otherwise
                                    error('Invalid setting for preconditioner_x')
                            end
                        end
                        
                        for itheta=1:Ntheta
                            for izeta=1:Nzeta
                                indices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                                addSparseBlock(indices, indices, -nuN*(GHat+iota*IHat)/(BHat(itheta,izeta)^2)*fieldTerm)
                            end
                        end
                        
                    end
                otherwise
                    error('Invalid setting for collisionOperator')
            end
            
            
            % --------------------------------------------------
            % Add constraints.
            % --------------------------------------------------
            
            switch constraintScheme
                case 0
                    % Do nothing.
                    
                case 1
                    
                    L=0;
                    for itheta=1:Ntheta
                        for izeta=1:Nzeta
                            colIndices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                            
                            rowIndex = matrixSize-1;
                            addSparseBlock(rowIndex, colIndices, (x2.*xWeights)' / (BHat(itheta,izeta)^2))
                            
                            rowIndex = matrixSize;
                            addSparseBlock(rowIndex, colIndices, (x2.*x2.*xWeights)' / (BHat(itheta,izeta)^2))
                        end
                    end
                    
                case 2
                    L=0;
                    spatialPart = 1./(BHat.*BHat);
                    spatialPart = reshape(spatialPart',[Ntheta*Nzeta,1])';
                    for ix=1:Nx
                        rowIndex = Nx*Nxi*Ntheta*Nzeta + ix;
                        colIndices = (ix-1)*Nxi*Ntheta*Nzeta + (1:(Ntheta*Nzeta));
                        addSparseBlock(rowIndex, colIndices, spatialPart)
                    end
                    
                otherwise
                    error('Invalid constraintScheme')
            end
            
            % --------------------------------------------------
            % Add sources.
            % --------------------------------------------------
            
            switch constraintScheme
                case 0
                    % Do nothing
                    
                case 1
                    xPartOfSource1 = (x2-5/2).*expx2;
                    xPartOfSource2 = (x2-3/2).*expx2;
                    
                    for ix=1:Nx
                        rowIndices = (ix-1)*Nxi*Ntheta*Nzeta + (1:(Ntheta*Nzeta));
                        
                        colIndex = matrixSize-1;
                        addSparseBlock(rowIndices, colIndex, xPartOfSource1(ix)*ones(Ntheta*Nzeta,1))
                        
                        colIndex = matrixSize;
                        addSparseBlock(rowIndices, colIndex, xPartOfSource2(ix)*ones(Ntheta*Nzeta,1))
                    end
                    
                case 2
                    
                    for ix=1:Nx
                        rowIndices = (ix-1)*Nxi*Ntheta*Nzeta + (1:(Ntheta*Nzeta));
                        colIndex = Nx*Nxi*Ntheta*Nzeta + ix;
                        addSparseBlock(rowIndices, colIndex, ones(Ntheta*Nzeta,1))
                    end
                otherwise
                    error('Invalid constraintScheme')
            end
            
            % --------------------------------------------------
            % End of adding entries to the matrix.
            % --------------------------------------------------
            
            switch whichMatrixToMake
                case 0
                    fprintf('Time to contruct preconditioner: %g seconds.\n',toc(matrixStartTime))
                    tic
                    preconditionerMatrix = createSparse();
                    fprintf('Time to sparsify preconditioner: %g seconds.\n',toc)
                case 1
                    fprintf('Time to contruct main matrix: %g seconds.\n',toc(matrixStartTime))
                    tic
                    matrix = createSparse();
                    fprintf('Time to sparsify main matrix: %g seconds.\n',toc)
                otherwise
                    error('Program should not get here')
            end
        end
        
        
        % ------------------------------------------------------
        % Finalize matrix
        % ------------------------------------------------------
        
        if ~ tryIterativeSolver
            preconditionerMatrix = matrix;
        end
        fprintf('Actual nnz of matrix: %d,  of preconditioner: %d,  predicted: %d\n',nnz(matrix),nnz(preconditionerMatrix),estimated_nnz_original)
        fprintf('Fraction of nonzero entries in matrix: %g,  in preconditioner: %g\n', nnz(matrix)/numel(matrix), nnz(preconditionerMatrix)/numel(preconditionerMatrix))
        fprintf('nnz(preconditioner)/nnz(matrix): %g\n', nnz(preconditionerMatrix)/nnz(matrix))
        
        assignin('base','mm',matrix)
        if tryIterativeSolver
            assignin('base','pm',preconditionerMatrix)
        end
        
        
        if tryIterativeSolver
            fprintf('LU-factorizing preconditioner...')
            tic
            [preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q] = lu(preconditionerMatrix);
            fprintf('done.  Took %g seconds.\n',toc)
        end
        
        
        function solnVector=preconditioner(rhsVector)
            solnVector = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * rhsVector)));
        end
        
        % ------------------------------------------------------
        % Solve the main linear system
        % ------------------------------------------------------
        
        if ~tryIterativeSolver
            fprintf('Applying sparse direct solver...\n')
            tic
            soln = matrix \ rhs;
            fprintf('Done. Time to solve system: %g seconds.\n',toc)
        else
            % Use an iterative Krylov-space solver.
            
            soln = zeros(size(rhs));
            numRHSs = size(rhs,2);
            
            for col=1:numRHSs
                if numRHSs > 1
                    fprintf('--- Solving linear system for RHS column %d of %d. ---\n',col, numRHSs)
                end
                attempt=0;
                keepTrying = true;
                x0 = zeros(matrixSize,1);
                while keepTrying
                    attempt = attempt+1;
                    tic
                    switch orderOfSolversToTry(attempt)
                        case 1
                            fprintf('Attempting iterative solve using BiCGStab...\n')
                            solverName = 'BiCGStab';
                            [soln0,fl0,rr0,it0,rv0]=bicgstab(matrix,rhs(:,col),tol,maxit,@preconditioner, [], x0);
                        case 2
                            solverName = 'BiCGStab(l)';
                            fprintf('Attempting iterative solve using BiCGStab(l)...\n')
                            [soln0,fl0,rr0,it0,rv0]=bicgstabl(matrix,rhs(:,col),tol,maxit,@preconditioner, [], x0);
                        case 3
                            fprintf('Attempting iterative solve using CGS...\n')
                            solverName = 'CGS';
                            [soln0,fl0,rr0,it0,rv0]=cgs(matrix,rhs(:,col),tol,maxit,@preconditioner, [], x0);
                        case 4
                            fprintf('Attempting iterative solve using GMRES...\n')
                            solverName = 'GMRES';
                            [soln0,fl0,rr0,it0,rv0]=gmres(matrix,rhs(:,col),restart,tol,maxit/restart,@preconditioner, [], x0);
                        case 5
                            fprintf('Attempting iterative solve using TFQMR...\n')
                            solverName = 'TFQMR';
                            [soln0,fl0,rr0,it0,rv0]=tfqmr(matrix,rhs(:,col),tol,maxit,@preconditioner, [], x0);
                        otherwise
                            error('Invalid setting for orderOfSolversToTry.')
                    end
                    switch fl0
                        case 0
                            fprintf('Converged!\n')
                        case 1
                            fprintf('Did not converge :(\n')
                        case 2
                            fprintf('Preconditioner was ill-conditioned\n')
                        case 3
                            fprintf('Stagnated :(\n')
                    end
                    fprintf('Time to apply solver: %g seconds.\n',toc)
                    if iteration==1
                        KrylovFigureHandle = figure(3 + figureOffset);
                    else
                        set(0, 'CurrentFigure', KrylovFigureHandle);
                    end
                    clf
                    semilogy(rv0/rv0(1),'-o');
                    xlabel('Iteration number');
                    ylabel('Relative residual');
                    title(['Convergence of Krylov solver ',solverName]);
                    drawnow
                    fprintf('Minimum residual: %g.\n',min(rv0)/rv0(1))
                    if fl0==0
                        keepTrying=false;
                        soln(:,col) = soln0;
                    else
                        if attempt >= numel(orderOfSolversToTry)
                            keepTrying=false;
                        else
                            x0 = soln0;
                            fprintf('Iterative solver failed, so trying again with backup solver.\n')
                        end
                    end
                end
                
                % If last iterative solver failed, use direct solver.
                if fl0 ~= 0
                    fprintf('Switching to direct solution since iterative solver(s) failed...\n')
                    tic
                    soln = matrix \ rhs;
                    fprintf('Done. Time to solve system: %g seconds.\n',toc)
                    break
                end
            end
            
            
        end
        
        fprintf('Total elapsed time: %g sec.\n',toc(startTimeForThisRun))

        computeOutputs()
        
        % ------------------------------------------------------
        % Calculate radial heat and particle fluxes
        % ------------------------------------------------------
        function computeOutputs()
            
            for col = 1:RHSSize
                
                if RHSSize > 1
                    fprintf('--- Analyzing solution vector %d of %d. ---\n',col, RHSSize)
                end
                
                switch constraintScheme
                    case 0
                        % Do nothing
                    case 1
                        sources = soln(end-1:end, col);
                        fprintf('Sources: %g,  %g\n',sources(1),sources(2))
                    case 2
                        sources = soln((Nx*Nxi*Ntheta*Nzeta+1):end, col);
                        fprintf('min source: %g,   max source: %g\n',min(sources),max(sources))
                    otherwise
                        error('Invalid constraintScheme')
                end
                
                
                densityPerturbation = zeros(Ntheta,Nzeta);
                flow = zeros(Ntheta,Nzeta);
                pressurePerturbation = zeros(Ntheta,Nzeta);
                
                particleFluxBeforeSurfaceIntegral = zeros(Ntheta,Nzeta);
                momentumFluxBeforeSurfaceIntegral = zeros(Ntheta,Nzeta);
                heatFluxBeforeSurfaceIntegral = zeros(Ntheta,Nzeta);
                
                densityPerturbationIntegralWeight = x.^2;
                flowIntegralWeight = x.^3;
                pressurePerturbationIntegralWeight = x.^4;
                
                particleFluxIntegralWeight = x.^4;
                momentumFluxIntegralWeight = x.^5;
                heatFluxIntegralWeight = x.^6;
                
                for itheta=1:Ntheta
                    for izeta = 1:Nzeta
                        L=0;
                        indices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                        fSlice = soln(indices, col);
                        densityPerturbation(itheta,izeta) = xWeights' * (densityPerturbationIntegralWeight .* fSlice);
                        pressurePerturbation(itheta,izeta) = xWeights' * (pressurePerturbationIntegralWeight .* fSlice);
                        particleFluxBeforeSurfaceIntegral(itheta,izeta) = (8/3)*xWeights' * (particleFluxIntegralWeight .* fSlice);
                        heatFluxBeforeSurfaceIntegral(itheta,izeta) = (8/3)*xWeights' * (heatFluxIntegralWeight .* fSlice);
                        
                        L=1;
                        indices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                        fSlice = soln(indices, col);
                        flow(itheta,izeta) = xWeights' * (flowIntegralWeight .* fSlice);
                        momentumFluxBeforeSurfaceIntegral(itheta,izeta) = (16/15)*xWeights' * (momentumFluxIntegralWeight .* fSlice);
                        
                        L=2;
                        indices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                        fSlice = soln(indices, col);
                        particleFluxBeforeSurfaceIntegral(itheta,izeta) = particleFluxBeforeSurfaceIntegral(itheta,izeta) ...
                            + (4/15)*xWeights' * (particleFluxIntegralWeight .* fSlice);
                        heatFluxBeforeSurfaceIntegral(itheta,izeta) = heatFluxBeforeSurfaceIntegral(itheta,izeta) ...
                            + (4/15)*xWeights' * (heatFluxIntegralWeight .* fSlice);
                        
                        L=3;
                        indices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                        fSlice = soln(indices, col);
                        momentumFluxBeforeSurfaceIntegral(itheta,izeta) = momentumFluxBeforeSurfaceIntegral(itheta,izeta) ...
                            + (4/35)*xWeights' * (momentumFluxIntegralWeight .* fSlice);
                        
                    end
                end
                
                densityPerturbation = 4*Delta*THat*sqrtTHat/(sqrtpi*psiAHat)*densityPerturbation;
                flow = 4*THat*THat/(3*sqrtpi) * flow;
                pressurePerturbation = (8/3)*Delta*THat*sqrtTHat/(sqrtpi*psiAHat)*pressurePerturbation;
                
                particleFluxBeforeSurfaceIntegral = -(THat^(5/2))*(GHat*dBHatdtheta-IHat*dBHatdzeta)./(sqrtpi*BHat.^3) ...
                    .* particleFluxBeforeSurfaceIntegral;
                
                momentumFluxBeforeSurfaceIntegral = -(THat^3)*(GHat*dBHatdtheta-IHat*dBHatdzeta)./(sqrtpi*BHat.^3) ...
                    .* momentumFluxBeforeSurfaceIntegral;
                
                heatFluxBeforeSurfaceIntegral = -(THat^(7/2))*(GHat*dBHatdtheta-IHat*dBHatdzeta)./(2*sqrtpi*BHat.^3) ...
                    .* heatFluxBeforeSurfaceIntegral;
                
                FSADensityPerturbation = (1/VPrimeHat) * thetaWeights' * (densityPerturbation./(BHat.^2)) * zetaWeights;
                FSAFlow = (1/VPrimeHat) * thetaWeights' * (flow./BHat) * zetaWeights;
                FSADensityPerturbation = (1/VPrimeHat) * thetaWeights' * (pressurePerturbation./(BHat.^2)) * zetaWeights;
                
                particleFlux = thetaWeights' * particleFluxBeforeSurfaceIntegral * zetaWeights;
                momentumFlux = thetaWeights' * momentumFluxBeforeSurfaceIntegral * zetaWeights;
                heatFlux = thetaWeights' * heatFluxBeforeSurfaceIntegral * zetaWeights;
                
                fprintf('FSADensityPerturbation:  %g\n',FSADensityPerturbation)
                fprintf('FSAFlow:                 %g\n',FSAFlow)
                fprintf('FSAPressurePerturbation: %g\n',FSAPressurePerturbation)
                fprintf('particleFlux:            %g\n',particleFlux)
                fprintf('momentumFlux:            %g\n',momentumFlux)
                fprintf('heatFlux:                %g\n',heatFlux)
                
                if RHSMode == 2
                    VPrimeHatWithG = VPrimeHat*(GHat+iota*IHat);
                    switch col
                        case 1
                            transportMatrix(1,1) = 4*(GHat+iota*IHat)*particleFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat^(3/2))*GHat);
                            transportMatrix(2,1) = 8*(GHat+iota*IHat)*heatFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat^(5/2))*GHat);
                            transportMatrix(3,1) = 2*nHat*FSAFlow/(GHat*THat);
                        case 2
                            transportMatrix(1,2) = 4*(GHat+iota*IHat)*particleFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*GHat);
                            transportMatrix(2,2) = 8*(GHat+iota*IHat)*heatFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*THat*GHat);
                            transportMatrix(3,2) = 2*FSAFlow/(GHat);
                        case 3
                            transportMatrix(1,3) = particleFlux*Delta*Delta*FSABHat2/(VPrimeHatWithG*GHat*psiAHat*omega);
                            transportMatrix(2,3) = 2*Delta*Delta*heatFlux*FSABHat2/(GHat*VPrimeHatWithG*psiAHat*THat*omega);
                            transportMatrix(3,3) = FSAFlow*Delta*Delta*sqrtTHat*FSABHat2/((GHat+iota*IHat)*2*psiAHat*omega*B0OverBBar);
                    end
                end
            end
            
            if RHSMode == 2
                format longg
                transportMatrix
                
                if geometryScheme==11 || geometryScheme==4
                  TSI=THat*1e3*1.6022e-19; %Assuming Tbar = 1 keV
                  if species=='p'
                    vT=sqrt(2*TSI/1.6726e-27);
                    q=1.6022e-19;
                  elseif species=='e'
                    vT=sqrt(2*TSI/9.1094e-31);
                    q=-1.6022e-19;
                  end
                  
                  transportMatrixSI=zeros(3);
                  transportMatrixSI(1:2,1:2)=transportMatrix(1:2,1:2) / ...
                      (-B0OverBBar*(GHat+iota*IHat)/GHat^2/TSI^2*vT* ...
                       q^2*dPsidr^2);
                  transportMatrixSI(1:2,3)=transportMatrix(1:2,3) / ...
                      (B0OverBBar*q/TSI/GHat*dPsidr);
                  transportMatrixSI(3,1:2)=transportMatrix(3,1:2) / ...
                      (-B0OverBBar*q/TSI/GHat*dPsidr);
                  transportMatrixSI(3,3)=transportMatrix(3,3) / ...
                      (B0OverBBar/vT/(GHat+iota*IHat));

                  disp(['In SI units for species ',species])
                  transportMatrixSI
                end
            end
            
            if testQuasisymmetryIsomorphism
                modifiedHeatFluxThatShouldBeConstant = heatFlux*abs(helicity_n/iota-helicity_l)/((helicity_n*IHat+helicity_l*GHat)^2);
                modifiedFSAFlowThatShouldBeConstant = FSAFlow*(helicity_n/iota-helicity_l)/(helicity_n*IHat+helicity_l*GHat);
                fprintf('   > Testing quasisymmetry isomorphism.\n')
                fprintf('   > Below are the two modified quantities that should be independent of helicity:\n')
                fprintf('   > Modified heat flux: %g\n',modifiedHeatFluxThatShouldBeConstant)
                fprintf('   > Modified FSA flow: %g\n',modifiedFSAFlowThatShouldBeConstant)
            end
            
            if programMode == 1
                figure(4+figureOffset)
                
                subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
                contourf(zeta,theta,densityPerturbation,numContours,'EdgeColor','none')
                colorbar
                xlabel('\zeta')
                ylabel('\theta')
                title('densityPerturbation')
                
                subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
                contourf(zeta,theta,flow,numContours,'EdgeColor','none')
                colorbar
                xlabel('\zeta')
                ylabel('\theta')
                title('flow')
                
                subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
                contourf(zeta,theta,pressurePerturbation,numContours,'EdgeColor','none')
                colorbar
                xlabel('\zeta')
                ylabel('\theta')
                title('pressurePerturbation')
            end
        end
        
        
        
        
        % --------------------------------------------------------
        % Below are some utilities for building sparse matrices.
        % --------------------------------------------------------
        
        function resetSparseCreator()
            sparseCreatorIndex=1;
            sparseCreator_i=zeros(estimated_nnz,1);
            sparseCreator_j=zeros(estimated_nnz,1);
            sparseCreator_s=zeros(estimated_nnz,1);
        end
        
        function addToSparse(i,j,s)
            n=numel(i);
            if n ~= numel(j)
                error('Error A');
            end
            if n ~= numel(s)
                error('Error B');
            end
            if any(i<1)
                error('Error Q: i<1');
            end
            if any(j<1)
                error('Error Q: j<1');
            end
            sparseCreator_i(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = i;
            sparseCreator_j(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = j;
            sparseCreator_s(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = s;
            sparseCreatorIndex = sparseCreatorIndex+n;
            if sparseCreatorIndex > estimated_nnz
                fprintf('Error! estimated_nnz is too small.\n')
            end
        end
        
        function addSparseBlock(rowIndices, colIndices, block)
            s=size(block);
            if (s(1) ~= numel(rowIndices)) || (s(2) ~= numel(colIndices))
                s
                size(rowIndices)
                size(colIndices)
                error('Error in addSparseBlock!')
            end
            [rows, cols, values] = find(block);
            addToSparse(rowIndices(rows),colIndices(cols),values)
        end
        
        function sparseMatrix = createSparse()
            fprintf('estimated nnz: %d   Actual value required: %d\n',estimated_nnz_original, sparseCreatorIndex)
            sparseMatrix = sparse(sparseCreator_i(1:(sparseCreatorIndex-1)), sparseCreator_j(1:(sparseCreatorIndex-1)), sparseCreator_s(1:(sparseCreatorIndex-1)), matrixSize, matrixSize);
            resetSparseCreator()
        end
        
        % ------------------------------------------------------
        % ------------------------------------------------------
        % Below are routines to set the magnetic geometry.
        % ------------------------------------------------------
        % ------------------------------------------------------
        
        function setNPeriods()
            switch geometryScheme
                case 1
                    NPeriods = max([1, helicity_n]);
                case {2,3}
                    NPeriods = 10;
                case 4
                    NPeriods = 5;
                case 10
                    fid = fopen(fort996boozer_file);
                    if fid<0
                        error('Unable to open file %s\n',fort996boozer_file)
                    end
                    try
                      NPeriods = fscanf(fid,'%d',1);
                      fclose(fid);
                    catch me
                      error('%s\n\nFile\n\t%s\ndoes not seem to be a valid vmec fort.996 output file.\n',...
                            me.message, fort996boozer_file)
                    end
                case 11
                    fid = fopen(JGboozer_file);
                    if fid<0
                        error('Unable to open file %s\n',JGboozer_file)
                    end
                    try
                        tmp_str=fgetl(fid);       %Skip comment line
                        while strcmp(tmp_str(1:2),'CC');
                            tmp_str=fgetl(fid);     %Skip comment line
                        end
                        header=fscanf(fid,'%d %d %d %d %f %f %f\n',7);
                        NPeriods = header(4);
                        fclose(fid);
                    catch me
                      error('%s\n\nFile\n\t%s\ndoes not seem to be a valid vmec .bc output file.\n',...
                            me.message, JGboozer_file)
                    end
                otherwise
                    error('Invalid setting for geometryScheme')
            end
        end
        
        function computeBHat()
            % Eventually, this subroutine should be expanded to allow more options
            
            [zeta2D, theta2D] = meshgrid(zeta,theta);

            switch geometryScheme
               case 1
                  % 2-helicity model:
                  BHarmonics_l = [1, helicity_l];
                  if helicity_n==0
                    BHarmonics_n = [0, 0];
                  else
                    BHarmonics_n = [0, 1];
                  end
                  BHarmonics_amplitudes = [epsilon_t, epsilon_h];
                  
               case 2
                  % LHD standard configuration.
                  % Values taken from Table 1 of
                  % Beidler et al, Nuclear Fusion 51, 076001 (2011).
                  iota = 0.4542;
                  BHarmonics_l = [1, 2, 1];
                  BHarmonics_n = [0, 1, 1];
                  BHarmonics_amplitudes = [-0.07053, 0.05067, -0.01476];
                  
                  B0OverBBar = 1; % (Tesla)
                  R0 = 3.7481; % (meters)
                  a = 0.5585; % (meters)
                  GHat = B0OverBBar * R0;
                  %IHat = GHat*3; % Change this to 0 eventually.
                  IHat = 0;
                  psiAHat = B0OverBBar*a^2/2;
                  
               case 3
                  % LHD inward-shifted configuration.
                  % Values taken from Table 1 of
                  % Beidler et al, Nuclear Fusion 51, 076001 (2011).
                  iota = 0.4692;
                  BHarmonics_l = [1, 2, 1, 0];
                  BHarmonics_n = [0, 1, 1, 1];
                  BHarmonics_amplitudes = [-0.05927, 0.05267, -0.04956, 0.01045];
                  
                  B0OverBBar = 1; % (Tesla)
                  R0 = 3.6024; % (meters)
                  a = 0.5400; % (meters)
                  GHat = B0OverBBar * R0;
                  IHat = 0;
                  psiAHat = B0OverBBar*a^2/2;
                  
               case 4
                  % W7-X Standard configuration
                  % Values taken from Table 1 of
                  % Beidler et al, Nuclear Fusion 51, 076001 (2011).
                  iota=0.8700;
                  BHarmonics_l = [0, 1, 1];
                  BHarmonics_n = [1, 1, 0];
                  BHarmonics_amplitudes = [0.04645, -0.04351, -0.01902];
                  
                  B0OverBBar = 3.089; % (Tesla)
                  %R0 = 5.5267; % (meters)
                  a = 0.5109; % (meters)
                  %psiAHat = -B0OverBBar*a^2/2;
                  GHat = -17.885;%B0OverBBar * R0;
                  IHat = 0;
                  psiAHat = -0.384935;
                  radius=0.2555; %m, radius of the flux surface
                  dPsidr=2*psiAHat/a*(radius/a);
                  %nuPrime=nuN*abs(GHat+iota*IHat)/B0OverBBar/sqrt(THat);
                  nuPrime=nuN*(GHat+iota*IHat)/B0OverBBar/sqrt(THat);

               case 10
                  fid = fopen(fort996boozer_file);
                  if fid<0
                      error('Unable to open file %s\n',fort996boozer_file)
                  end
                  
                  % File description:
                  % 1st line: 2 integers:     nfp,ns
                  % 2nd line: 4 real numbers: aspect,rmax,rmin,betaxis
                  % 3rd line: 3 integers:     mboz, nboz, mnboz
                  % 4th line: 7 real numbers: iota,pres,beta,phip,phi,bvco,buco
                  %
                  % Then, you have 'mnboz' lines.
                  % If 'mn' is a dummy integer variable that goes from 1 to mnboz,
                  % for each value of mn you read
                  %
                  % m(mn),n(mn),bmn(mn),rmnc(mn),zmns(mn)pmns(m,n),gmn(mn)
                  try
                    header=fscanf(fid,'%d %d\n %f %f %f %f\n %d %d %d %f %f %f %f %f %f %f',16);
                    mnboz=header(9);
                    modes =fscanf(fid,'%d %d %g %g %g %g %g',[7,mnboz]);
                    fclose(fid);
                    
                    % scalar values
                    %Nper = header(1); %number of field periods
                    iota = header(10);
                    Ihat = header(16);  % Covariant theta comp. of B, known as I in sfincs (in meter * Tesla)
                    Ghat = header(15);  % Covariant phi comp. of B, known as G in sfincs (in meter * Tesla)
                    % Note that the flux at the separatrix is not stored in the
                    % file, so we set PsiAHat in the Physics parameters
                    % section in the beginning of the program
                    
                    % mode amplitudes
                    if modes(1,1)==0 && modes(2,1)==0 
                      B0OverBBar=modes(3,1); %The B00 component in Tesla
                    else
                      error('The first fort996boozer_file entry is not the B00 component')
                    end
                    BHarmonics_l = modes(1,2:end);
                    BHarmonics_n = modes(2,2:end);
                    BHarmonics_amplitudes = modes(3,2:end)/B0OverBBar; % Store the values normalised to the B00 component. 
                    
                  catch me
                    error('%s\n\nFile\n\t%s\ndoes not seem to be a valid vmec fort.996 output file.\n',...
                        me.message, fort996boozer_file)
                  end
               case 11
              
                  fid = fopen(JGboozer_file);
                  if fid<0
                      error('Unable to open file %s\n',JGboozer_file)
                  end
                  
                  try
                      tmp_str=fgetl(fid);
                      while strcmp(tmp_str(1:2),'CC');
                          tmp_str=fgetl(fid); %Skip comment line
                      end
                      header=fscanf(fid,'%d %d %d %d %f %f %f\n',7);
                      fgetl(fid);  %Skip variable name line
                      
                      NPeriods = header(4);
                      psiAHat  = header(5)/2/pi; %Convert the flux from Tm^2 to Tm^2/rad
                      a        = header(6);      %minor radius %m
                      
                      max_no_of_modes=500;
                      modesm_new=NaN*zeros(1,max_no_of_modes);
                      modesn_new=NaN*zeros(1,max_no_of_modes);
                      modesb_new=NaN*zeros(1,max_no_of_modes);
                      normradius_new=-inf;
                      no_of_modes_new=NaN;
                      iota_new=NaN;
                      G_new=NaN;
                      I_new=NaN;
                      end_of_file=0;
                      
                      while (normradius_new<normradius_wish) && not(end_of_file)
                          normradius_old=normradius_new;
                          no_of_modes_old=no_of_modes_new;
                          modesm_old=modesm_new;
                          modesn_old=modesn_new;
                          modesb_old=modesb_new;
                          iota_old=iota_new;
                          G_old=G_new;
                          I_old=I_new;
                          
                          fgetl(fid);
                          surfheader=fscanf(fid,'%f %f %f %f %f %f\n',6);
                          
                          normradius_new=sqrt(surfheader(1)); %r/a=sqrt(psi/psi_a)
                          iota_new=surfheader(2);
                          G_new=surfheader(3)*NPeriods/2/pi*(4*pi*1e-7); %Tesla*meter
                          I_new=surfheader(4)/2/pi*(4*pi*1e-7);          %Tesla*meter
                          
                          fgetl(fid); %Skip units line
                          proceed=1;
                          modeind=0;
                          while proceed
                              tmp_str=fgetl(fid);
                              if length(tmp_str)==1
                                  if tmp_str==-1 %End of file has been reached
                                      proceed=0;
                                      end_of_file=1;
                                  end
                              elseif tmp_str(8)=='s' %Next flux surface has been reached
                                  proceed=0;
                              else
                                  tmp=sscanf(tmp_str,'%d %d %f %f %f %f',6);
                                  if abs(tmp(6))>min_Bmn_to_load
                                      modeind=modeind+1;
                                      %if modeind > max_no_of_modes %Unnecessary to check this in matlab
                                      %  error(' modeind > max_no_of_modes !')
                                      %end
                                      modesm_new(modeind)=tmp(1);
                                      modesn_new(modeind)=tmp(2);
                                      modesb_new(modeind)=tmp(6);
                                  end
                              end
                          end
                          no_of_modes_new=modeind;
                          modesm_new(no_of_modes_new+1:end)=NaN;
                          modesn_new(no_of_modes_new+1:end)=NaN;
                          modesb_new(no_of_modes_new+1:end)=NaN;
                      end
                      fclose(fid);
                  catch me
                      error('%s\n\nFile\n\t%s\ndoes not seem to be a valid .bc geometry file.\n',...
                          me.message, JGboozer_file)
                  end

                  [~,minind]=min([(normradius_old-normradius_wish)^2,...
                                  (normradius_new-normradius_wish)^2]);
                  if minind==1
                    BHarmonics_l = modesm_old(1:no_of_modes_old);
                    BHarmonics_n = modesn_old(1:no_of_modes_old);
                    BHarmonics_amplitudes = modesb_old(1:no_of_modes_old);
                    iota=iota_old;
                    GHat=G_old;
                    IHat=I_old;
                    normradius=normradius_old;
                  else %minind=2
                    BHarmonics_l = modesm_new(1:no_of_modes_new);
                    BHarmonics_n = modesn_new(1:no_of_modes_new);
                    BHarmonics_amplitudes = modesb_new(1:no_of_modes_new);
                    iota=iota_new;
                    GHat=G_new;
                    IHat=I_new;
                    normradius=normradius_new;
                  end
                  disp(['The calculation is performed for radius ' ...
                        ,num2str(normradius*a),' m , r/a=',num2str(normradius)])
                  
                  m0inds=find(BHarmonics_l==0);
                  n0m0inds=find(BHarmonics_n(m0inds)==0);
                  if isempty(n0m0inds)
                    error(' B00 component is missing!')
                  end
                  nm00ind=m0inds(n0m0inds);
                  B0OverBBar=BHarmonics_amplitudes(nm00ind); %Assumes \bar{B}=1T
                  BHarmonics_amplitudes=[BHarmonics_amplitudes(1:nm00ind-1),...
                                         BHarmonics_amplitudes(nm00ind+1:end)]...
                                        /B0OverBBar;
                  BHarmonics_l = [BHarmonics_l(1:nm00ind-1),...
                                  BHarmonics_l(nm00ind+1:end)];
                  BHarmonics_n = [BHarmonics_n(1:nm00ind-1),...
                                  BHarmonics_n(nm00ind+1:end)];
                  
                  dPsidr=2*psiAHat/a*normradius;
                  %nuPrime=nuN*abs(GHat+iota*IHat)/B0OverBBar/sqrt(THat);
                  nuPrime=nuN*(GHat+iota*IHat)/B0OverBBar/sqrt(THat);
                  
                otherwise
                  error('Invalid setting for geometryScheme')
            end
            
            NHarmonics = numel(BHarmonics_amplitudes);
            BHat = B0OverBBar * ones(Ntheta,Nzeta);
            dBHatdtheta = zeros(Ntheta,Nzeta);
            dBHatdzeta = zeros(Ntheta,Nzeta);
            for i=1:NHarmonics
              BHat = BHat + B0OverBBar * BHarmonics_amplitudes(i) *...
                     cos(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods * zeta2D);
              dBHatdtheta = dBHatdtheta - B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_l(i) *...
                  sin(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods * zeta2D);
              dBHatdzeta = dBHatdzeta + B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_n(i) * NPeriods *...
                  sin(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods * zeta2D);            
            end
        end  
    end
end
        



