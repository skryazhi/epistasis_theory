%%% This code test whether the results of Theorem 1 hold when the
%%% the effects of mutations are not infinitesimally small

%%% Update: 3 January 2021 (SK)
%%% (1) Added a probability for a reaction to have zero rate constant in genKineticParamsLin.m
%%% (2) Fixed error in fixed point calculation for linear function fit

%%% Update: 22 Dec 2020 (SK). Added new plot combining results for
%%% topological classes Mio and Mi

%%% Update: 22 Oct 2020 (SK)

% As a result of these calculations structure eps is generated with the
% following fields:
%
% eps.xvec is an neps by 1 vector of epsx values
%
% eps.yvec is an neps by ParamSetCnt by nMut array with epsy values for
% each parameter and mutational perturbation combination
%
% eps.ft is a ParamSetCnt by nMut cell array with the cfit objects
% describing the hyperbolic function that maps epsx onto epsy. If the fit
% failed, nan is returned
%
% eps.yfit is a neps by ParamSetCnt by nMut array that contains the values
% of the fit at the sampled values of epsx
%
% eps.fp is an ParamSetCnt by nMut array which contains the location of
% the fixed point of the epsy as function of epsx if it exists, otherwise
% contains nan
%
% eps.slope_fp is an ParamSetCnt by nMut array which contains the slopes
% of the epsy as a function of epsx at the fixed point values, if the fixed
% point exists, otherwise contains nan
%
% eps.TF is a logical 4 by ParamSetCnt by nMut array which evaluates how
% the function epsy behaves at 0 and 1. row 1 = is eps x = 0 admissible (it
% may be inadmissible if the resulting x^{AB} becomes negative)?; row 2 = is
% epsy ( 0 ) <= 0? row 3 = is eps x = 1 admissible?; row 4 = is eps y ( 1 ) >= 1?


%% Parameter definitions

% 1. Specify class of the single-marked module (Mio or Mi)
ModuleClass = 'Mi'; 

% 2. Specify what to do
IFGENPARAMS = false; % generate new parameters or use previously generated ones?
IFCALCULATE = false; % carry out a new calculations or plot only?
IFPLOT = false;     % plot during calculations? Useful to check curve fits

% 3. Specify how many random parameter combinations to explore
ParamSetCnt = 1000;

% Mutational perturbation sizes
mutMat = [
    -0.01 -0.01;
    -0.1 -0.1;
    -0.5 -0.5;
    0.01  0.01;
    0.1   0.1;
    0.5   0.5;
    -0.01 0.01;
    -0.1  0.1;
    -0.5  0.5];

% Vector of epsilon x values to test:
eps.xvec = (-1:0.2:2)';
TF01 = eps.xvec >= 0 & eps.xvec <= 1;

% Distribution of kinetic parameters
DistrType = 'Exp0';
DistrParam  = [0.25 , 1]; % 25% probability of having no edge, parameter of the exponential is 1

nMut = size(mutMat,1 );
neps = length(eps.xvec);

%%% Adjacency matrix:
if strcmp(ModuleClass, 'Mio')
    AdjMat = ones(3,3) - eye(3,3);
    MutRxn = [1 3];
elseif strcmp(ModuleClass, 'Mi')    
    AdjMat = ones(4,4) - eye(4,4);
    MutRxn = [3 4];
else
    error('Please specify a single-marked module class Mio or Mi');
end


%% Generating random kinetic parameters
if IFGENPARAMS
    VmMat0 = nan( size(AdjMat,1), size(AdjMat,2), ParamSetCnt);
    KeqMat0 = nan( size(AdjMat,1), size(AdjMat,2), ParamSetCnt);
    
    for iParam = 1:ParamSetCnt
        
        ISMUTRXNPRESENT = false;
        
        while ~ISMUTRXNPRESENT 
            [VmMatTMP, KeqMatTMP] = genKineticParamsLin( AdjMat, DistrType, DistrParam );
            if VmMatTMP( MutRxn(1), MutRxn(2) ) > 0
                ISMUTRXNPRESENT = true;
            end
        end
        
        VmMat0(:,:,iParam) = VmMatTMP;
        KeqMat0(:,:,iParam) = KeqMatTMP;
    end
    
    filename = sprintf('params_%s.mat', ModuleClass);
    save(filename, 'VmMat0', 'KeqMat0');
end
clear iParam ISMUTRXNPRESENT VmMatTMP KeqMatTMP;


%% Calculating how eps y depends on eps x
if IFCALCULATE
    
    % Loading previously stored kinetic parameters unless new parameters
    % were just generated
    if ~IFGENPARAMS
        filename = sprintf('params_%s.mat', ModuleClass);
        load(filename);
    end
    
    close all;

    eps.yvec = nan(neps, ParamSetCnt, nMut);
    eps.slope_fp = nan(ParamSetCnt, nMut);
    eps.fp = nan(ParamSetCnt, nMut);
    eps.ft = cell(ParamSetCnt, nMut);
    eps.yfit = nan(neps, ParamSetCnt, nMut);
    eps.TF = false(4,ParamSetCnt, nMut);
    
    for iMut = 1:size(mutMat,1)        
        dxA = mutMat(iMut,1);
        dxB = mutMat(iMut,2);
        
        for iParam = 1:ParamSetCnt
            
            fprintf('iMut = %d, iParam = %d: ', iMut, iParam);
            
            KeqMat = KeqMat0(:,:,iParam);
            VmMat = VmMat0(:,:,iParam);
            
            % wildtype effective rate
            yWT = get_effective_rate( VmMat , KeqMat );            
            
            % effective rate of A mutant and its relative effect
            VmMat = VmMat0(:,:,iParam);
            VmMat(MutRxn(1),MutRxn(2)) = VmMat(MutRxn(1),MutRxn(2)) * (1 + dxA );
            VmMat(MutRxn(2),MutRxn(1)) = VmMat(MutRxn(1),MutRxn(2)) / KeqMat(MutRxn(1),MutRxn(2));
            dyA = get_effective_rate( VmMat , KeqMat ) / yWT - 1;

            % effective rate of B mutant and its relative effect
            VmMat = VmMat0(:,:,iParam);
            VmMat(MutRxn(1),MutRxn(2)) = VmMat(MutRxn(1),MutRxn(2)) * (1 + dxB );
            VmMat(MutRxn(2),MutRxn(1)) = VmMat(MutRxn(1),MutRxn(2)) / KeqMat(MutRxn(1),MutRxn(2));            
            dyB = get_effective_rate( VmMat , KeqMat ) / yWT - 1;
            
            % effective rate of AB mutant, its relative effect and the resulting epsy
            for ieps = 1:neps
                VmMat = VmMat0(:,:,iParam);
                VmMat(MutRxn(1),MutRxn(2)) = VmMat(MutRxn(1),MutRxn(2)) * (1 + dxA + dxB + 2 * eps.xvec(ieps) * dxA * dxB );
                VmMat(MutRxn(2),MutRxn(1)) = VmMat(MutRxn(1),MutRxn(2)) / KeqMat(MutRxn(1),MutRxn(2));
                
                if VmMat(MutRxn(1),MutRxn(2)) < 0
                    continue;
                end
                
                dyAB = get_effective_rate( VmMat , KeqMat ) / yWT - 1;
                
                eps.yvec(ieps,iParam,iMut) = (dyAB - dyA - dyB)/2/dyA/dyB;
            end
            
            % evaluating the bounds for epsx (to later check that the fixed
            % point is within them)
            if dxA * dxB > 0
                epsxmin = ( -1 - dxA - dxB)/2/dxA/dxB;
                epsxmax = Inf;
            else
                epsxmin = -Inf;
                epsxmax = ( -1  - dyA - dyB)/2/dxA/dxB;                
            end
                        
            
            % plotting the calculated epsy values
            if IFPLOT
                clf;
                hold on, box on;
                plot(eps.xvec, eps.yvec(:,iParam, iMut), 'o' );
            end
            
            % Fitting the hyperbolic function
            [eps.ft{iParam,iMut}, f, slope, fp] = ...
                fit_hyperbolic(eps.xvec, eps.yvec(:,iParam,iMut) );
            
            % Checking the behavior at epsx = 0 and epsx = 1
            epsy01 = eps.yvec(TF01,iParam,iMut);
            if ~isnan( epsy01(1) )
                eps.TF(1,iParam,iMut) = true;
                if epsy01(1) <= 0
                    eps.TF(2,iParam,iMut) = true;
                end
            end
            
            if ~isnan( epsy01(end) )
                eps.TF(3,iParam,iMut) = true;
                if epsy01(end) >= 1
                    eps.TF(4,iParam,iMut) = true;
                end
            end

            if ~isa(f, 'function_handle')
                fprintf('no fit\n');
                continue;
            end
            
            % Calculating fitted values at the sampled epsx points
            TF = ~isnan( eps.yvec(:,iParam,iMut) );
            eps.yfit(TF,iParam,iMut) = f( eps.xvec(TF) );
            
            % Recording the fixed point
            if ~isnan( fp(1) )
                eps.fp(iParam,iMut) = fp(1);
            end
            
            % Checking that the fixed point is in the admissible range
            if eps.fp(iParam,iMut) < epsxmin || eps.fp(iParam,iMut) > epsxmax
                eps.fp(iParam,iMut) = nan;
            end
            
            % Recording the slope at the fixed point
            if ~isnan(eps.fp(iParam,iMut))
                eps.slope_fp(iParam,iMut) = slope( fp(1) );
            end
            
            % Plot the fit
            if IFPLOT
                plot(eps.xvec, eps.yfit(:,iParam, iMut), '-' );
                plot(eps.xvec, eps.xvec, 'k--');
                
                pause(0.25);
            end
            fprintf('fp = %.2f\n', fp(1) );
        end
    end
    clear iMut dxA dxB iParam yWT VmMat KeqMat dyA dyB dyAB ieps TF f slope fp;
    
    filename = sprintf('params_%s.mat', ModuleClass);
    save(filename, 'VmMat0', 'KeqMat0', 'eps');
end







%% DATA SUMMARIES AND VISUALIZATION

ModuleClassList = {'Mio', 'Mi'};

% Will display only bar-eps values between these numbers
xMin = 0; 
xMax = 4;
xLen = xMax - xMin;
xTickLabels = arrayfun(@(x) sprintf('(%d, %d]', x, x+1),  xMin:(xMax-1), 'UniformOutput', false);
xTickLabels = [{sprintf('< %d', xMin)} xTickLabels {sprintf('> %d', xMax), 'No f.p.'}];

close all;

edges = xMin:1:xMax; % edges for the histogram for the fixed point location
X = (xMin-1):(xLen+1);      % X locations for the histogram plot
Y = nan(1,xLen+3); % histrogram counts for a single topological class;
Ytot = zeros( nMut, xLen+3); % total histrogram counts (across both topological classes);

for iMC = 1:2
    ModuleClass = ModuleClassList{iMC};
    
    filename = sprintf('params_%s.mat', ModuleClass);
    load(filename);
    
    %%%% Check if the hyperbolic fits are good, if haven't checked during
    %%%% the calcuation. Otherwise skip
    %
    % iMut = 3;
    %
    % for iParam = 1:100
    %     clf;
    %     hold on, box on;
    %     plot(eps.xvec, eps.yvec(:,iParam, iMut), 'o' );
    %
    %     plot(eps.xvec, eps.yfit(:,iParam, iMut), '-' );
    %     plot(eps.xvec, eps.xvec, 'k--');
    %
    %     pause(0.25);
    % end


    %% Print some summary statistics
    fprintf('=== Module class %s ===\n', ModuleClass);
    fprintf('Mut. perturb set\tno f.p.\tf.p. < %d\tf.p. > 1\tf.p. > %d\tf.p.slope < 1\n', xMin, xMax)
    fprintf('---\n');
    
    for iMut = 1:nMut
        fprintf('%d\t%d\t%d\t%d\t%d\t%d\n',...
            iMut,... % mutational perturbation set
            nnz(isnan(eps.fp(:,iMut))),... % existence of fixed point
            nnz( eps.fp(:,iMut) < xMin ), ... % fixed point < 0
            nnz( eps.fp(:,iMut) > 1 ), ... % fixed point > 1
            nnz( eps.fp(:,iMut) > xMax ), ... % fixed point > 3
            nnz( eps.slope_fp(:,iMut) < 1 ) ... % fixed point stable ?
            );
        fprintf('---\n');
    end
    clear iMut;
    
    % For cases, where the fixed point does not exist, checking whether phi(x)
    % is above or below the diagonal (for positive x)
    for iMut = 1:nMut
        if nnz(isnan(eps.fp(:,iMut))) > 0
            ixvec = find(isnan(eps.fp(:,iMut)));
            cntpos = 0;
            for ix = ixvec'
                % checkin the value of the phi function at 0:
                if feval(eps.ft{1,1}, 0) < 0 
                    cntpos = cntpos+1;
                end
            end
            fprintf('In mut. perturb. set %d, phi(x) < x in %d out of %d cases with no f.p.\n',...
                iMut, cntpos, nnz(isnan(eps.fp(:,iMut))) );
        end
    end
    clear iMut ixvec cntpos ix;
        
    
    %% Plotting histograms for individual topological classes (Figure 4 supplements)
    
    %%% Specify the following dimensions:
    fdim.spwa = 5; % subplotwidth in cm
    fdim.spha = 5; % subplotheight in cm
    
    fdim.nx = 3; % number of panels along the horizontal dimension
    fdim.ny = 3; % number of panels along the vertical dimension
    
    fdim.xma = [1 0.2]; % left right horizontal margin in cm
    fdim.yma = [1.5 0.2]; % bottom top vertical margin cm
    
    fdim.dxa = 0.2; % horizontal distance between panels in cm
    fdim.dya = 0.2; % vertical distance between panels in cm
    
    fdim.tickfs = 8; % tick font size
    fdim.labelfs = 10; % axis label font size
    
    %%% These will be computed automatically:
    fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
    fdim.fh = fdim.spha * fdim.ny + fdim.dya * (fdim.ny - 1) + sum(fdim.yma);
    
    fdim.spwr = fdim.spwa / fdim.fw;
    fdim.sphr = fdim.spha / fdim.fh;
    fdim.xmr = fdim.xma / fdim.fw;
    fdim.ymr = fdim.yma / fdim.fh;
    fdim.dxr = fdim.dxa / fdim.fw;
    fdim.dyr = fdim.dya / fdim.fh;
    
    figure;
    
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [0 0 fdim.fw fdim.fh]);
        
    fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
    fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );
    
    % fdim
    
    for iMut = 1:nMut
        ix = mod(iMut-1,3)+1;
        iy = floor((iMut-1)/3)+1;
        
        subplot('Position', [fdim.spxvec(ix) fdim.spyvec(iy) fdim.spwr fdim.sphr]),
        hold on, box on;
        set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
        
        [Y(2:xLen+1), edges] = histcounts( eps.fp(:,iMut), edges );
        Y(1) = nnz( eps.fp(:,iMut) <= xMin );
        Y(xLen+2) = nnz( eps.fp(:,iMut) > xMax );
        Y(xLen+3) = nnz( isnan(eps.fp(:,iMut)) );
        Ytot(iMut,:) = Ytot(iMut,:) + Y;
        Y = Y/size(eps.fp,1);
        
        bar(X,Y, 0.9, 'FaceColor', 0.4*[1 1 1], 'EdgeColor', 'none', 'ShowBaseLine', 'off');        
        
        set(gca, 'XLim', [xMin-1 xMax+1] + [-0.5 0.5] + (xLen+2)*0.05*[-1 +1],...
            'XTick', xMin-1:xMax+1 );
        set(gca, 'YLim', [-0.05, 1.05]);
        
        if iy == 3
            set(gca, 'XTickLabel', xTickLabels);
            if ix == 2
                xlabel('Fixed point position', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
            end
        else
            set(gca, 'XTickLabel', {});
        end
        
        xtickangle(45);
        
        if ix > 1
            set(gca, 'YTickLabel', {});
        else
            if iy == 2
                ylabel('Fraction of sampled modules', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
            end
        end
        
        lim = axis;
        
        if ix == 1
            s = sprintf('\\delta^Ay_\\mu = %.2f\n\\delta^By_\\mu = %.2f', mutMat(iMut,1), mutMat(iMut,2)  );
        else
            s = sprintf('\\delta^Ay_\\mu = %.1f\n\\delta^By_\\mu = %.1f', mutMat(iMut,1), mutMat(iMut,2)  );
        end
        
        text( lim(2) - (lim(2)-lim(1))*0.02, lim(4) - (lim(4)-lim(3))*0.02,...
            s, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
            'FontName', 'Helvetica', 'FontSize', fdim.labelfs);
        
    end    
    clear fdim cc iMut ix iy Y lim;         
    
end

%% Plotting Figure 4

Ytot = Ytot ./ repmat(sum(Ytot,2), 1, xLen+3);

%%% Specify the following dimensions:
fdim.spwa = 5; % subplotwidth in cm
fdim.spha = 5; % subplotheight in cm

fdim.nx = 3; % number of panels along the horizontal dimension
fdim.ny = 1; % number of panels along the vertical dimension

fdim.xma = [1 0.2]; % left right horizontal margin in cm
fdim.yma = [1.5 0.2]; % bottom top vertical margin cm

fdim.dxa = 0.2; % horizontal distance between panels in cm
fdim.dya = 0.2; % vertical distance between panels in cm

fdim.tickfs = 8; % tick font size
fdim.labelfs = 10; % axis label font size

%%% These will be computed automatically:
fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
fdim.fh = fdim.spha * fdim.ny + fdim.dya * (fdim.ny - 1) + sum(fdim.yma);

fdim.spwr = fdim.spwa / fdim.fw;
fdim.sphr = fdim.spha / fdim.fh;
fdim.xmr = fdim.xma / fdim.fw;
fdim.ymr = fdim.yma / fdim.fh;
fdim.dxr = fdim.dxa / fdim.fw;
fdim.dyr = fdim.dya / fdim.fh;

figure;

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 fdim.fw fdim.fh]);

fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );

% fdim

for iMut = 4:6
    ix = iMut-3;
    iy = 1;
    
    subplot('Position', [fdim.spxvec(ix) fdim.spyvec(iy) fdim.spwr fdim.sphr]),
    hold on, box on;
    set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
    
    bar(X,Ytot(iMut,:), 0.9, 'FaceColor', 0.4*[1 1 1], 'EdgeColor', 'none', 'ShowBaseLine', 'off');
    
    set(gca, 'XLim', [xMin-1 xMax+1] + [-0.5 0.5] + (xLen+2)*0.05*[-1 +1],...
            'XTick', xMin-1:xMax+1 );
    set(gca, 'YLim', [-0.05, 1.05]);
    
    set(gca, 'XTickLabel', xTickLabels);
    
    xtickangle(45);
    
    if ix > 1
        set(gca, 'YTickLabel', {});
        if ix == 2
            xlabel('Fixed point position', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
        end
    else
        ylabel('Fraction of sampled modules', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
    end
    
    lim = axis;
       
    s = sprintf('+%.0f%% effect', mutMat(iMut,1)*100);
    
    text( lim(2) - (lim(2)-lim(1))*0.02, lim(4) - (lim(4)-lim(3))*0.02,...
        s, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontName', 'Helvetica', 'FontSize', 12);
    
end

fprintf('%.1f%% cases with +50%% effect have bar eps > 1\n', sum(Ytot(6,3:xLen+2))*100 );

clear fdim cc iMut ix iy X Ytot edges lim s xLen xMin xMax xTickLabels;



