function m20140424_01_launchSfincsW7XFullProfileTest()

%profilesFilename='C:\Users\landreman\Box Sync\MATLAB\dataForSfincs\profiles.txt';
%JGboozer_file = 'C:\Users\landreman\Box Sync\MATLAB\dataForSfincs\w7x-sc1-ecb2-scaled-truncated.bc';
profilesFilename='/global/scratch2/sd/landrema/W7X_files_from_Yuriy/profiles.txt';
JGboozer_file = '/global/scratch2/sd/landrema/W7X_files_from_Yuriy/w7x-sc1-ecb2-scaled-truncated.bc';

NErs_at_each_psi = 5;
Er_range = 4; % in kV/m
Er_bunchedness = 0.7; % Must be a value in (0,1).

%desired_r_over_a_values = 0.2:0.2:1;
%desired_r_over_a_values = 0.1:0.2:1;
desired_r_over_a_values = 0.1:0.1:1;

Npsi = numel(desired_r_over_a_values);

% ---------------------------------------------
% Read in profile info from Yuriy:

data=importdata(profilesFilename);

r_over_a = data.data(:,1);
ne = data.data(:,3);
ni = data.data(:,4);
Te = data.data(:,5);
Ti = data.data(:,6);
Er = data.data(:,7);
particleFlux = data.data(:,8);
Qe = data.data(:,10);
Qi = data.data(:,11);
je = data.data(:,12);
ji = data.data(:,13);

psiN = r_over_a .^ 2;

r_over_a_fine = linspace(0,1,200);

figure(1)
clf
plot(r_over_a,Er,'.-')
hold on

for ipsi = 1:Npsi
    
    normradius_wish = desired_r_over_a_values(ipsi);
    
    dashWidth=0.1;
    
    % *************************************************************************
    % First, find the radii that are available in the .bc file, and select the
    % closest available radius:
    % *************************************************************************
    %normradius = normradius_wish;
    
    
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
        
        normradius_new=-inf;
        end_of_file=0;
        
        while (normradius_new<normradius_wish) && not(end_of_file)
            normradius_old=normradius_new;
            
            fgetl(fid);
            surfheader=fscanf(fid,'%f %f %f %f %f %f\n',6);
            
            normradius_new=sqrt(surfheader(1)); %r/a=sqrt(psi/psi_a)
            
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
                elseif not(isempty(find(tmp_str=='s'))) %Next flux surface has been reached
                    proceed=0;
                else
                    tmp=sscanf(tmp_str,'%d %d %f %f %f %f',6);
                end
            end
        end
        fclose(fid);
    catch me
        error('%s\n\nFile\n\t%s\ndoes not seem to be a valid .bc geometry file.\n',...
            me.message, JGboozer_file)
    end
    
    [~,minind]=min([(normradius_old-normradius_wish)^2,...
        (normradius_new-normradius_wish)^2]);
    if minind==1
        normradius=normradius_old;
    else %minind=2
        normradius=normradius_new;
    end
    %disp(['The calculation is performed for radius ' ...
    %    ,num2str(normradius*a),' m , r/a=',num2str(normradius)])
    
    fprintf('r/a requested: %g, used: %g\n',normradius_wish, normradius)
    
    
    % *************************************************************************
    % Done scanning the .bc file for available radii.
    % *************************************************************************
    
    radiusDirectoryName = ['r_over_a_',num2str(normradius)];
    mkdir(radiusDirectoryName)
    cd(radiusDirectoryName)
    
    interpolationAlgorithm = 'cubic';
    ne_to_use = interp1(r_over_a, ne, normradius, interpolationAlgorithm);
    Te_to_use = interp1(r_over_a, Te, normradius, interpolationAlgorithm);
    Ti_to_use = interp1(r_over_a, Ti, normradius, interpolationAlgorithm);
    
    ne_fine = interp1(r_over_a, ne, r_over_a_fine, interpolationAlgorithm);
    Te_fine = interp1(r_over_a, Te, r_over_a_fine, interpolationAlgorithm);
    Ti_fine = interp1(r_over_a, Ti, r_over_a_fine, interpolationAlgorithm);
    
    delta = -1e-3;

    data = ne;
    dnedrho_to_use = (interp1(r_over_a, data, normradius+delta, interpolationAlgorithm) - interp1(r_over_a, data, normradius, interpolationAlgorithm))/delta;
    
    data = Te;
    dTedrho_to_use = (interp1(r_over_a, data, normradius+delta, interpolationAlgorithm) - interp1(r_over_a, data, normradius, interpolationAlgorithm))/delta;

    data = Ti;
    dTidrho_to_use = (interp1(r_over_a, data, normradius+delta, interpolationAlgorithm) - interp1(r_over_a, data, normradius, interpolationAlgorithm))/delta;


    %{    
    [f,gof,out] = fit(r_over_a, Te, 'smoothingspline');
    Te_fine = f(r_over_a_fine);
    Te_to_use = f(normradius);
    dTedrho_to_use = differentiate(f,normradius);
    
    [f,gof,out] = fit(r_over_a, Ti, 'smoothingspline');
    Ti_fine = f(r_over_a_fine);
    Ti_to_use = f(normradius);
    dTidrho_to_use = differentiate(f,normradius);
    
    [f,gof,out] = fit(r_over_a, ne, 'smoothingspline');
    ne_fine = f(r_over_a_fine);
    ne_to_use = f(normradius);
    dnedrho_to_use = differentiate(f,normradius);
    %}
    
    %{
    [f,gof,out] = fit(psiN, ne, 'smoothingspline');
    ne_test = f(normradius^2);
    dnedpsiN_test = differentiate(f,normradius^2);
    %}
    
    dnedpsiN_to_use = dnedrho_to_use / (2*normradius);
    dTedpsiN_to_use = dTedrho_to_use / (2*normradius);
    dTidpsiN_to_use = dTidrho_to_use / (2*normradius);

    lnLambda = 25.3-1.15*log10(ne_to_use*1e14)+2.3*log10(Te_to_use*1000);
    epsilon0 = 8.8542e-12;
    mBar = 1.6726e-27;
    e = 1.6022e-19;
    nu = 4*sqrt(2*pi)*(1e20)*(e^4)*lnLambda/(3*(4*pi*epsilon0)^2*sqrt(mBar)*(1000*e)^(3/2));
    vBar = sqrt(2*1000*e/mBar);
    nu_n_to_use = nu/vBar;
    
    %{
    figure(1)
    clf
    numRows = 2;
    numCols = 2;
    plotNum = 1;
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(r_over_a,ne,'x')
    hold on
    plot(r_over_a,ni,'+r')
    plot(r_over_a_fine,ne_fine,'-')
    plot(normradius, ne_to_use, '.')
    plot([normradius-dashWidth, normradius+dashWidth], [ne_to_use - dashWidth*dnedrho_to_use, ne_to_use + dashWidth*dnedrho_to_use], ':')
    xlabel('r/a')
    ylabel('n [10^{20} m^{-3}]')
    plot(normradius, ne_to_use, '.')
    legend('n_e from Yuriy','n_i from Yuriy','Fit','Value at radius to use','Fit to gradient')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(r_over_a,Te,'x')
    hold on
    plot(r_over_a,Ti,'+r')
    plot(r_over_a_fine,Te_fine,'-')
    plot(r_over_a_fine,Ti_fine,'-r')
    xlabel('r/a')
    ylabel('T [keV]')
    plot(normradius, Te_to_use, '.')
    plot(normradius, Ti_to_use, '.r')
    plot([normradius-dashWidth, normradius+dashWidth], [Te_to_use - dashWidth*dTedrho_to_use, Te_to_use + dashWidth*dTedrho_to_use], ':')
    plot([normradius-dashWidth, normradius+dashWidth], [Ti_to_use - dashWidth*dTidrho_to_use, Ti_to_use + dashWidth*dTidrho_to_use], ':r')
    legend('T_e from Yuriy','T_i from Yuriy','Fit to T_e','Fit to T_i','Value ot T_e at radius to use','Value of T_i at radius to use','Fit to gradient of T_e','Fit to gradient of T_i')
    
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(psiN,ne,'x')
    hold on
    plot(normradius^2, ne_to_use, '+')
    plot(normradius^2, ne_test, '.')
    plot([normradius^2-dashWidth, normradius^2+dashWidth], [ne_to_use - dashWidth*dnedpsiN_to_use, ne_to_use + dashWidth*dnedpsiN_to_use], ':')
    plot([normradius^2-dashWidth, normradius^2+dashWidth], [ne_test - dashWidth*dnedpsiN_test, ne_test + dashWidth*dnedpsiN_test], ':r')
    xlabel('\psi_N')
    ylabel('n [10^{20} m^{-3}]')
    legend('n_e from Yuriy','ne to use','ne test','Fit to d/d rho','Fit to d/d psi_N')
    %}
    
    %{
    fprintf('\n')
    fprintf('Numbers for species namelist in sfincs:\n')
    fprintf('nHats = %g %g\n',ne_to_use,ne_to_use)
    fprintf('dnHatdpsiNs = %g %g\n',dnedpsiN_to_use,dnedpsiN_to_use)
    fprintf('THats = %g %g\n',Ti_to_use,Te_to_use)
    fprintf('dTHatdpsiNs = %g %g\n',dTidpsiN_to_use,dTedpsiN_to_use)
    fprintf('\n')
    
    fprintf('Er from DKES:                 %g\n',interp1(r_over_a, Er, normradius, 'cubic'));
    fprintf('Particle flux from DKES:      %g\n',interp1(r_over_a, particleFlux, normradius, 'cubic'));
    fprintf('Electron heat flux from DKES: %g\n',interp1(r_over_a, Qe, normradius, 'cubic'));
    fprintf('Ion heat flux from DKES:      %g\n',interp1(r_over_a, Qi, normradius, 'cubic'));
    fprintf('Electron j_bs from DKES:      %g\n',interp1(r_over_a, je, normradius, 'cubic'));
    fprintf('Ion j_bs from DKES:           %g\n',interp1(r_over_a, ji, normradius, 'cubic'));
    %}
    
    Er_DKES = interp1(r_over_a, Er, normradius, 'cubic');
    Er_to_try = tan(pi/2*linspace(-1,1,NErs_at_each_psi)*Er_bunchedness) / tan(pi/2*Er_bunchedness) * Er_range/2 + Er_DKES;
    plot(ones(size(Er_to_try))*normradius, Er_to_try, '.r')
    
    a = 0.52625; % meters
    dPhiHatdpsiNs = - a/(2 * normradius) * Er_to_try;
    
    for iPhysics = 1:4
        switch iPhysics
            case 1
                physicsDirName='FP_full';
            case 2
                physicsDirName='FP_DKES';
            case 3
                physicsDirName='PAS_full';
            case 4
                physicsDirName='PAS_DKES';
            otherwise
                error('Should not get here')
        end
        mkdir(physicsDirName)
        cd(physicsDirName)
        
        for iEr = 1:NErs_at_each_psi
            ErDirName = ['Er_',num2str(Er_to_try(iEr))];
            mkdir(ErDirName)
            cd(ErDirName)
            
            % Copy pbs job file, editing job name:
            
            fid_out = fopen('job.sfincsScan','w');
            if fid_out == -1
                error('Could not open output job file')
            end
            fid_in = fopen('../../../job.sfincsScan','r');
            if fid_in == -1
                error('Could not open input job.sfincsScan file')
            end
            
            line = fgetl(fid_in);
            %fprintf(fid_out,line);
            fprintf(fid_out,'%s\n',line);

            jobName = [radiusDirectoryName,'_',physicsDirName,'_',ErDirName];
            fprintf(fid_out,'#PBS -N %s\n',jobName);
            
            while ~ feof(fid_in)
                line = fgetl(fid_in);
                fprintf(fid_out,'%s\n',line);
            end
            
            fclose(fid_out);
            fclose(fid_in);
            
            % Done copying PBS job file.
            
            % Next, copy input.namelist file, editing parameters as needed.

            fid_out = fopen('input.namelist','w');
            if fid_out == -1
                error('Could not open output namelist file')
            end
            fid_in = fopen('../../../input.namelist','r');
            if fid_in == -1
                error('Could not open input namelist file')
            end
            
            while ~ feof(fid_in)
                line = strtrim(fgetl(fid_in));
                
                if lineContains(line,'normradius_wish')
                    line = ['normradius_wish = ',num2str(normradius)];
                end
                
                if lineContains(line,'dPhiHatdpsiN')
                    line = ['dPhiHatdpsiN = ',num2str(dPhiHatdpsiNs(iEr))];
                end
                
                if lineContains(line,'nHats')
                    line = sprintf('nHats = %g %g', ne_to_use, ne_to_use);
                end
                
                if lineContains(line,'dnHatdpsiNs')
                    line = sprintf('dnHatdpsiNs = %g %g', dnedpsiN_to_use, dnedpsiN_to_use);
                end
                
                if lineContains(line,'THats')
                    line = sprintf('THats = %g %g', Ti_to_use, Te_to_use);
                end
                
                if lineContains(line,'dTHatdpsiNs')
                    line = sprintf('dTHatdpsiNs = %g %g', dTidpsiN_to_use, dTedpsiN_to_use);
                end
                
                if lineContains(line,'nu_n')
                    line = sprintf('nu_n = %g', nu_n_to_use);
                end
                
                if lineContains(line,'includeXDotTerm')
                    if iPhysics==2 || iPhysics==4
                        line = 'includeXDotTerm = .f.';
                    else
                        line = 'includeXDotTerm = .t.';
                    end
                end
                
                if lineContains(line,'includeElectricFieldTermInXiDot')
                    if iPhysics==2 || iPhysics==4
                        line = 'includeElectricFieldTermInXiDot = .f.';
                    else
                        line = 'includeElectricFieldTermInXiDot = .t.';
                    end
                end
                
                if lineContains(line,'useDKESExBDrift')
                    if iPhysics==2 || iPhysics==4
                        line = 'useDKESExBDrift = .t.';
                    else
                        line = 'useDKESExBDrift = .f.';
                    end
                end
                
                if lineContains(line,'collisionOperator')
                    if iPhysics==1 || iPhysics==2
                        line = 'collisionOperator = 0';
                    else
                        line = 'collisionOperator = 1';
                    end
                end
                
                fprintf(fid_out,'%s\n',line);
            end
            
            fclose(fid_out);
            fclose(fid_in);

            % Submit job!
            !qsub job.sfincsScan
            
            cd('..')
        end
        
        cd('..')
    end
    
    cd('..')
end
xlabel('r/a')
ylabel('E_r [kV/m]')


    function tf = lineContains(line,pattern)
        line = strtrim(line);
        
        if numel(line)<1
            tf = false;
            return
        end
        
        if line(1)=='!'
            tf=false;
            return;
        end
        
        if numel(line) < numel(pattern)
            tf = false;
            return
        end
        
        if strncmp(line,pattern,numel(pattern))
            % Get the rest of the line after the pattern.
            % This test is required so the line 'Nxi=2' does not contain
            % 'Nx'.
            temp = strtrim(line((numel(pattern)+1):end));
            if temp(1) == '='
                tf = true;
            else
                tf = false;
            end
        else
            tf = false;
        end
    end
end