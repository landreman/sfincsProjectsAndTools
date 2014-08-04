function Jparallel()

% --------------------------------------------------
% Geometry parameters:
% --------------------------------------------------

geometryScheme = 11;
% 1 = Two-helicity model
% 2 = Three-helicity approximation of the LHD standard configuration
% 3 = Four-helicity approximation of the LHD inward-shifted configuration
% 4 = Three-helicity approximation of the W7-X Standard configuration
% 10= Read the boozer coordinate data from the file specified as "fort996boozer_file" below
% 11= Read the boozer coordinate data from the file specified as "JGboozer_file" below

% Additional parameters used only when geometryScheme=1:
% B = BBar * B0OverBBar * [1 + epsilon_t * cos(theta) + epsilon_h * cos(helicity_l * theta - helicity_n * zeta)]
B0OverBBar = 1.55;
epsilon_t = 0.1;
%epsilon_t = 0.178539;
epsilon_h = 0.0;
helicity_l = 0;
helicity_n = 0;
% iota is the rotational transform = 1 / (safety factor q)
iota = 0.397543;
%iota = 1.0;
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
THat = 1.0;
nHat = 1.0;

% Radial electric field.
dPhiHatdpsi = 0;

% The following two quantities matter for RHSMode=1 but not for RHSMode=2:
dTHatdpsi = -0.0;
dnHatdpsi = -0.2;
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
nuN =  100/5.79;% AM
%nuhat_Helander = 1e-5; nuN = 4/3/sqrt(pi)/vbar*nuhat_Helander; %nuhat as in eq(3.48) Helander&Sigmar 

% If testQuasisymmetryIsomorphism is true, the value of nuN is changed so the physical collisionality
% stays constant as the helicity is changed.

nuPrime = NaN;%1;  %This value is over-written for geometryScheme = 4, 11
nuPrime = 100;

% --------------------------------------------------
% Numerical resolution parameters:
% --------------------------------------------------

% For each of the quantities below, the 'Converged' value is used except
% when that quantity is being varied in a convergence scan, in which case
% each value in the array that follows (e.g. Nthetas, NLs, etc.) is used.

% Number of grid points in the poloidal direction.
% Memory and time requirements DO depend strongly on this parameter.
NthetaConverged = 200;
%Nthetas = floor(linspace(25,35,3));

% Number of grid points in the toroidal direction
% (per identical segment of the stellarator.)
% Memory and time requirements DO depend strongly on this parameter.
NzetaConverged = 200;%25
%Nzetas = floor(linspace(3,5,1));

% Number of Legendre polynomials used to represent the distribution
% function.
% Memory and time requirements DO depend strongly on this parameter.
% The value of this parameter required for convergence depends strongly on
% the collisionality. At high collisionality, this parameter can be as low
% as ~ 5. At low collisionality, this parameter may need to be many 10s or
% even > 100 for convergence.
%NxiConverged = 80;

%Nxis = floor(linspace(40,80,3));

% Tolerance used to define convergence of the Krylov solver.
% This parameter does not affect memory requirements but it does affect the
% time required for solution.
log10tolConverged = 5.5;
%log10tols = 4.5:1:6.5;

tol = 10^(-log10tolConverged);

% Number of Cosine Harmonics to keep in 1/B^2 = h = SUM(m,n) h_{m,n} * cos(m * theta - n * zeta)]
% h_{0,0} always included
% OBS!! N_Harmonics_m SHOULD NOT BE LARGER THAN NthetaConverged AND
% N_Harmonics_n SHOULD NOT BE LARGER THAN NzetaConverged
N_Harmonics_m = 100;
N_Harmonics_n = 40;



% --------------------------------------------------
% Other numerical parameters:
% --------------------------------------------------

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



% --------------------------------------------------------
% --------------------------------------------------------
% First, set up the grids, differentiation matrices, and
% integration weights for each coordinate.
% --------------------------------------------------------
% --------------------------------------------------------

Ntheta=NthetaConverged;
Nzeta=NzetaConverged;
        

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
% iteration=0;
% if iteration>1
%     fprintf('********************************************************************\n')
% end
fprintf('Ntheta = %d,  Nzeta = %d, tol = %g\n',Ntheta, Nzeta,tol)

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
    %             if preconditioner_theta_remove_cyclic
    %                 ddtheta_preconditioner(1,end) = 0;
    %                 ddtheta_preconditioner(end,1) = 0;
    %             end
    
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
    %             if preconditioner_zeta_remove_cyclic
    %                 ddzeta_preconditioner(1,end) = 0;
    %                 ddzeta_preconditioner(end,1) = 0;
    %             end
end

% Evaluate the magnetic field and its derivatives on the
% (theta,zeta) grid:
computeBHat()

% Compute a few quantities related to the magnetic field
VPrimeHat = thetaWeights' * (1./BHat.^2) * zetaWeights;
FSABHat2 = 4*pi*pi/VPrimeHat;

theta
zeta
theta2D
zeta2D
BHat

thetaWeights

zetaWeights

% ------------------------------------------------------------
% Calculate parallel current u from cosine harmonics in 1/B^2.
% ------------------------------------------------------------

hHat = 1./(BHat.^2);

N_Harmonics_m = ceil(abs(N_Harmonics_m));
N_Harmonics_n = ceil(abs(N_Harmonics_n));

if N_Harmonics_m < 1 && N_Harmonics_n < 1
    error('Invalid setting for number of harmonics in 1/B^2')
end

if N_Harmonics_m > Ntheta || N_Harmonics_n > Nzeta
    fprintf('WARNING! Number of Harmonics used in Fourier expansion larger than grid size.\n')
end

% NHarmonics_u = N_Harmonics_m*N_Harmonics_n;

% hHarmonics_m = zeros(NHarmonics_u);
% hHarmonics_n = zeros(NHarmonics_u);

% uHat is u normalized to Rbar/Bbar, i.e. uHat = u * Bbar / Rbar
onesMatrix = ones(Ntheta,Nzeta);
uHat = zeros(Ntheta,Nzeta);
hHatApprox = 1/(4*pi^2)*(thetaWeights' * hHat * zetaWeights) * onesMatrix;

% Derivatives of uHat
duHatdtheta = zeros(Ntheta,Nzeta);
duHatdzeta = zeros(Ntheta,Nzeta);

NPeriods

for m=0:N_Harmonics_m-1
    for n=0:N_Harmonics_n-1
        if m==0 && n==0 %Don't need this component
            continue
        end
        % DO INTEGRATION
        % uHarmonics_amplitude = integral2(fun,0,2*pi,0,2*pi) / (2*pi^2);
        % u_{n,m} = (J*m + I*n)/(n - iota*m) * h_{n,m}
        hHatHarmonics_amplitude = 2/(Ntheta*Nzeta)  * sum(sum(cos(m * theta2D - n * NPeriods * zeta2D).*hHat ));
        %hHatHarmonics_amplitude = 1/(2*pi^2) * thetaWeights' * (cos(m * theta2D - n * NPeriods * zeta2D).*hHat) * zetaWeights;
        % uHatHarmonics_amplitude = (GHat*m + IHat*n)/(n - iota*m) * 2/(Ntheta*Nzeta)  * sum(sum(cos(m * theta2D - n * NPeriods * zeta2D).*hHat));
        
        uHatHarmonics_amplitude = iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude;
        
%         % TEST
%         if abs(n - iota*m)<10^(-2)
%             n - iota*m
%             uHatHarmonics_amplitude
%         end
        
        % ADD TO CURRENT u
        hHatApprox = hHatApprox + hHatHarmonics_amplitude * cos(m * theta2D - n * NPeriods * zeta2D);
        uHat = uHat + uHatHarmonics_amplitude * cos(m * theta2D - n * NPeriods * zeta2D);
        duHatdtheta = duHatdtheta - uHatHarmonics_amplitude * m * sin(m * theta2D - n * NPeriods * zeta2D);
        duHatdzeta = duHatdzeta + uHatHarmonics_amplitude * n * NPeriods * sin(m * theta2D - n * NPeriods * zeta2D);
    end
end

close all;

figure
surf(theta2D, zeta2D, BHat); % view(30, 30);
shading interp;
light;
lighting phong;
title('Bhat', 'FontName', 'Courier', 'FontSize', 14);

figure
surf(theta2D, zeta2D, hHat); % view(30, 30);
shading interp;
light;
lighting phong;
title('1/Bhat^2', 'FontName', 'Courier', 'FontSize', 14);

figure
surf(theta2D, zeta2D, hHatApprox); % view(30, 30);
shading interp;
light;
lighting phong;
title('1/Bhat^2 approximation', 'FontName', 'Courier', 'FontSize', 14);

figure
surf(theta2D, zeta2D, uHat); % view(30, 30);
shading interp;
light;
lighting phong;
title('Parallel current', 'FontName', 'Courier', 'FontSize', 14);

figure
surf(theta2D, zeta2D, duHatdtheta); % view(30, 30);
shading interp;
light;
lighting phong;
title('du/dtheta', 'FontName', 'Courier', 'FontSize', 14);

figure
surf(theta2D, zeta2D, duHatdzeta); % view(30, 30);
shading interp;
light;
lighting phong;
title('du/dzeta', 'FontName', 'Courier', 'FontSize', 14);


% Simakov analytical axisymmetry + Flux function
uHatAxi = - 1*(GHat * hHat - 1*GHat * 1/(4*pi^2)*(thetaWeights' * hHat * zetaWeights));%/iota;

 figure
% hold on
 surf(theta2D, zeta2D, uHatAxi); % view(30, 30);
% 
% hold off

test=0;
% COMPUTE QUANTITES FOR CALCULATING FLUX-SURFACE AVERAGE HEAT FLUX
% VPrimeHat = thetaWeights' * (1./BHat.^2) * zetaWeights;
% onesMatrix = ones(Ntheta,Nzeta);
BhatsquareAvg = thetaWeights' * onesMatrix * zetaWeights / VPrimeHat
% uHatBhatsquareAvg = thetaWeights' * (uHat./BHat.^2) * zetaWeights / VPrimeHat
uHatBhatsquareAvg = thetaWeights' * (uHat +  test*onesMatrix) * zetaWeights / VPrimeHat
% uHatBhatsquareAvg = thetaWeights' * (uHat./BHat.^2) * zetaWeights / VPrimeHat;
% uHatsquareBhatsquareAvg = thetaWeights' * (uHat.^2./BHat.^2) * zetaWeights / VPrimeHat
uHatsquareBhatsquareAvg = thetaWeights' * ((uHat +  test*onesMatrix).^2) * zetaWeights / VPrimeHat

% FLUX-SURFACE AVERAGE DERIVATIVES
GradParallelLnBhat = (iota*dBHatdtheta + dBHatdzeta)/(iota*IHat + GHat);
GradParallelBhat = BHat.*GradParallelLnBhat;
GradParalleluHatBhatsquare = BHat/(iota*IHat + GHat) .* (iota*BHat.^2.*duHatdtheta + 2*iota*BHat.*(uHat +  test*onesMatrix).*dBHatdtheta + BHat.^2.*duHatdzeta + 2*BHat.*(uHat +  test*onesMatrix).*dBHatdzeta);

%GradParallelBhatAvg = thetaWeights' * ((iota*dBHatdtheta + dBHatdzeta) ./ (BHat * (iota*IHat + GHat))) * zetaWeights / VPrimeHat;
%GradParallelLnBhatAvg = thetaWeights' * ((iota*dBHatdtheta + dBHatdzeta) ./ (BHat.^2 * (iota*IHat + GHat))) * zetaWeights / VPrimeHat;

%GradParalleluHatBhatsquareAvg = thetaWeights' * (1./(BHat * (iota*IHat + GHat)) .* (iota*BHat.^2.*duHatdtheta + 2*iota*BHat.*uHat.*dBHatdtheta + BHat.^2.*duHatdzeta + 2*BHat.*uHat.*dBHatdzeta)) * zetaWeights / VPrimeHat;

G1Hat = 1*(thetaWeights' * (GradParallelLnBhat.*GradParalleluHatBhatsquare./BHat.^2) * zetaWeights)^2 / ((thetaWeights' * ((GradParallelBhat).^2 ./BHat.^2) * zetaWeights) * VPrimeHat)     -    thetaWeights' * ((GradParalleluHatBhatsquare./BHat).^2 ./BHat.^2) * zetaWeights / VPrimeHat
% The Heat flux
% Modify this
%HeatFluxAvg = sqrt(2)/5 * (iota*IHat + GHat) * nHat * VPrimeHat /(iota^2 * sqrt(THat)) * dTHatdpsi * (uHatBhatsquareAvg^2/BhatsquareAvg - uHatsquareBhatsquareAvg) * nu_nbar

B0OverBBar
iota
GHat
IHat

nuPrime

Factor = sqrt(2)/5 * 8/(iota^2 * GHat^2) * (B0OverBBar)^2 * (uHatBhatsquareAvg^2/BhatsquareAvg - uHatsquareBhatsquareAvg)

minusL22 = - Factor * nuPrime

VPrimeHat

FactorL11 = 0.96*3/(2*sqrt(2)) * (iota*IHat + GHat)^2/(GHat^2*iota^2) * G1Hat

L11 = FactorL11 / nuPrime


L31 = uHatBhatsquareAvg / (iota*GHat) - 1*BhatsquareAvg/(2*iota*GHat) * (thetaWeights' * (GradParallelLnBhat.*GradParalleluHatBhatsquare./BHat.^2) * zetaWeights) / (thetaWeights' * ((GradParallelBhat).^2 ./BHat.^2) * zetaWeights) 



% hHarmonics_amplitudes

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
