%%% This code test whether the results of Theorem 1 hold when the
%%% the effects of mutations are not infinitesimally small

%%% Update: 3 January 2021 (SK)
%%% Added a probability for a reaction to have zero rate constant in genKineticParamsLin.m

%%% Update 22 Dec 2020 (SK)

%%% Update 23 Oct 2020 (SK)

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

filename = 'AllSerParGenRxnSets.mat';

if exist(filename, 'file') == 2
    load(filename);
else
    error('Run discover_generating_topologies.m to create data file AllSerParGenRxnSets.mat that contains all strictly serial and strictly parallel topologies');
end

TopoRelation = 'parallel'; % 'parallel' or 'serial';

% 2. Specify what to do
IFGENPARAMS = false; % generate new parameters or use previously generated ones?
IFCALCULATE = false; % carry out a new calculations or plot only?
IFLOAD = true;       % load pre-calculated data?

% 3. Specify how many random parameter combinations to explore
ParamSetCnt = 10000;

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


% Distribution of kinetic parameters
DistrType = 'Exp0';
DistrParam  = [0.25 , 1]; % 25% probability of having no edge, parameter of the exponential is 1

% Precision threshold:
MutTol = 1e-5; % mutational effects with absolute value below this threshold are considered to be zero
EpsTol = 1e-2; % epsilon values close to zero and 1 are evaluated with this tolerance

nMut = size(mutMat,1 );


%% Generating random kinetic parameters
if IFGENPARAMS
    
    if strcmp(TopoRelation, 'parallel')
        A = P;
    elseif strcmp(TopoRelation, 'serial')
        A = S;
    else
        error('Specify TopoRelation as either ''parallel'' or ''serial''');
    end
    
    TopoClasses = fieldnames( A );
    
    %%% Structure Params holds the descriptions of all generating
    %%% topologies
    Params.TC = {};      % Topological class
    Params.n_int = [];   % Number of internal metabolites in this class
    Params.id = [];      % ID of the generating topology
    Params.AdjMat = {};  % Adjacency matrix
    
    
    %%% First, identify the structure of adjacency matrices
    for iTC = 1:length(TopoClasses)
        TC = TopoClasses{iTC};

        nGenRxnSet = size( A.(TC).rmrxns_nosym, 1);
        
        % If there are no relations of the specified type in this
        % topological class, move to the next class
        if nGenRxnSet == 0
            continue;
        end
                       
        AdjMatF = ones( ClassDef.(TC).n+2 ) - eye( ClassDef.(TC).n+2 );
                
        if nGenRxnSet == 1 && size( A.(TC).rmrxns_nosym{1}, 1 ) == 0 % The full reaction set is a generating reaction set
            Params.TC = [Params.TC ; {TC} ];
            Params.n_int = [Params.n_int; ClassDef.(TC).n];
            Params.id = [Params.id ; 0];
            Params.AdjMat = [Params.AdjMat ; {AdjMatF} ];
        else
            % Cycle through all generating reaction sets
            for iid = 1:nGenRxnSet
                Params.TC = [Params.TC ; {TC} ];
                Params.n_int = [Params.n_int; ClassDef.(TC).n];
                Params.id = [Params.id ; iid];
                
                AdjMat = AdjMatF;
                
                % Remove all required edges to obtain the current generating set
                for ib = 1:size( A.(TC).rmrxns_nosym{iid}, 1 )
                    im1 = A.(TC).rmrxns_nosym{iid}(ib,1);
                    im2 = A.(TC).rmrxns_nosym{iid}(ib,2);
                    AdjMat(im1,im2) = 0;
                    AdjMat(im2,im1) = 0;
                end
                
                Params.AdjMat = [Params.AdjMat ; {AdjMat} ];
            end
        end
    end
    clear TopoClasses iTC TC nGenRxnSet AdjMatF AdjMat iid ib im1 im2;
    clear ClassList P S;
    
    Params.nGT = length( Params.TC );
    
    
    %%% Second, generate parameters based on adjacency matrices
    Params.VmMat0 = cell( Params.nGT, 1);
    Params.KeqMat0 = cell( Params.nGT, 1);
    
    fprintf('--- Generating random parameters ---\n');
    for iGT = 1:Params.nGT
        
        fprintf('Generating topology %d out of %d\n', iGT, Params.nGT);
        
        AdjMat = Params.AdjMat{iGT};        
        
        Params.VmMat0{iGT} = nan( size(AdjMat,1), size(AdjMat,2), ParamSetCnt);
        Params.KeqMat0{iGT} = nan( size(AdjMat,1), size(AdjMat,2), ParamSetCnt);
        
        MutRxn = nan(2,2);
        MutRxn(1,:) = ClassDef.(Params.TC{iGT}).a;
        MutRxn(2,:) = ClassDef.(Params.TC{iGT}).b;

        
        for iParam = 1:ParamSetCnt
            
            ISMUTRXNPRESENT = false;
            
            while ~ISMUTRXNPRESENT
                [VmMat, KeqMatTMP] = genKineticParamsLin( AdjMat, DistrType, DistrParam );

                if VmMat( MutRxn(1,1), MutRxn(1,2) ) > 0 && VmMat( MutRxn(2,1), MutRxn(2,2) ) > 0
                    ISMUTRXNPRESENT = true;
                end
            end

            Params.VmMat0{iGT}(:,:,iParam) = VmMat;
            Params.KeqMat0{iGT}(:,:,iParam) = KeqMatTMP;
        end
    end
    clear iGT AdjMat iParam ISMUTRXNPRESENT VmMatTMP KeqMatTMP MutRxn AdjMat;
    
    filename = sprintf('params_%s.mat', TopoRelation);
    save(filename, 'Params', 'ClassDef');
end


%% Calculating eps y
if IFCALCULATE
    
    % Loading previously stored kinetic parameters unless new parameters
    % were just generated
    if ~IFGENPARAMS
        filename = sprintf('params_%s.mat', TopoRelation);
        load(filename);
    end
    
    close all;

    eps = nan( ParamSetCnt , nMut, Params.nGT );
    
    fprintf('--- Calculating epistasis ---\n');
    for iGT = 1:Params.nGT
        
        fprintf('Generating topology %d out of %d\n', iGT, Params.nGT);
        
        MutRxn = nan(2,2);
        MutRxn(1,:) = ClassDef.(Params.TC{iGT}).a;
        MutRxn(2,:) = ClassDef.(Params.TC{iGT}).b;
        
        for iParam = 1:ParamSetCnt                        
            
            KeqMat = Params.KeqMat0{iGT}(:,:,iParam);
            VmMat = Params.VmMat0{iGT}(:,:,iParam);
            
            % wildtype effective rate
            yWT = get_effective_rate( VmMat , KeqMat );

            for iMut = 1:nMut
                dxA = mutMat(iMut,1);
                dxB = mutMat(iMut,2);
                
                % effective rate of A mutant and its relative effect
                VmMat = Params.VmMat0{iGT}(:,:,iParam);
                VmMat(MutRxn(1,1),MutRxn(1,2)) = VmMat(MutRxn(1,1),MutRxn(1,2)) * (1 + dxA );
                VmMat(MutRxn(1,2),MutRxn(1,1)) = VmMat(MutRxn(1,1),MutRxn(1,2)) / KeqMat(MutRxn(1,1),MutRxn(1,2));
                dyA = get_effective_rate( VmMat , KeqMat ) / yWT - 1;
                
                if abs(dyA) < MutTol
                    continue;
                end
                
                % effective rate of B mutant and its relative effect
                VmMat = Params.VmMat0{iGT}(:,:,iParam);
                VmMat(MutRxn(2,1),MutRxn(2,2)) = VmMat(MutRxn(2,1),MutRxn(2,2)) * (1 + dxB );
                VmMat(MutRxn(2,2),MutRxn(2,1)) = VmMat(MutRxn(2,1),MutRxn(2,2)) / KeqMat(MutRxn(2,1),MutRxn(2,2));
                dyB = get_effective_rate( VmMat , KeqMat ) / yWT - 1;
                
                if abs(dyB) < MutTol
                    continue;
                end
                
                % effective rate of AB mutant, its relative effect and the resulting epsy
                VmMat(MutRxn(1,1),MutRxn(1,2)) = VmMat(MutRxn(1,1),MutRxn(1,2)) * (1 + dxA );
                VmMat(MutRxn(1,2),MutRxn(1,1)) = VmMat(MutRxn(1,1),MutRxn(1,2)) / KeqMat(MutRxn(1,1),MutRxn(1,2));                
                
                dyAB = get_effective_rate( VmMat , KeqMat ) / yWT - 1;

                eps(iParam,iMut,iGT) = (dyAB - dyA - dyB)/2/dyA/dyB;
            end
        end
    end     
            
    clear iGT MutRxn iParam KeqMat VmMat yWT iMut dxA dxB dxAB;
        
    for iGT = 1:Params.nGT
        fprintf('Generating topology %d:\t <= 0\t| (0,1)\t| >= 1\n', iGT);
        for iMut = 1:nMut
            fprintf('Mut. pert %d:\t', iMut);
            fprintf('\t %d', nnz( eps(:,iMut,iGT) <= EpsTol ) );
            fprintf('\t| %d', nnz( eps(:,iMut,iGT) > EpsTol & eps(:,iMut,iGT) < 1 - EpsTol ) );
            fprintf('\t| %d\n', nnz( eps(:,iMut,iGT) >= 1 - EpsTol ) );
        end
    end

    
    filename = sprintf('params_%s.mat', TopoRelation);
    save(filename, 'Params', 'ClassDef', 'eps');
end


%% PLOTTING

close all;

if IFLOAD
    filename = sprintf('params_%s.mat', TopoRelation);
    load(filename);
end

GTID = Params.TC;
for iGT = 1:Params.nGT
    if Params.id(iGT) == 0
        GTID{iGT} = sprintf('%sF', GTID{iGT});
    else
        if strcmp(TopoRelation, 'parallel')
            GTID{iGT} = sprintf('%sP%d', GTID{iGT}, Params.id(iGT));
        elseif strcmp(TopoRelation, 'serial')
            GTID{iGT} = sprintf('%sS%d', GTID{iGT}, Params.id(iGT));
        end
    end
end

%%% Plot histrogram of fixed point locations
if strcmp(TopoRelation, 'parallel')
    fprintf('There are %d (%.1g%%) cases with eps y > 0 and %d (%.1g%%) cases with eps y >= 1\n',...
        nnz(eps > EpsTol ), nnz(eps > EpsTol )/nnz( ~isnan(eps) ) * 100,...
        nnz(eps >= 1-EpsTol ), nnz(eps >= 1-EpsTol )/nnz( ~isnan(eps) )*100  );
    
elseif strcmp(TopoRelation, 'serial')
    fprintf('There are %d (%.1g%%) cases with eps y < 1 and %d (%.1g%%) cases with eps y <= 0\n',...
        nnz(eps < 1 - EpsTol ), nnz(eps < 1-EpsTol )/nnz( ~isnan(eps) ) * 100,...
        nnz(eps <= EpsTol ), nnz(eps <= EpsTol )/nnz( ~isnan(eps) )*100  );
end

fprintf('Only cases between 0 and 1 are plotted\n');

X = (1:Params.nGT)';
EpsCnt = nan( Params.nGT , nMut ); % # of cases with 0 < eps < 1
TotCnt = nan( Params.nGT , nMut ); % total number of cases where epistasis could be computed

for iMut = 1:nMut    
    for iGT = 1:Params.nGT
        EpsCnt( iGT, iMut) = nnz( eps(:,iMut,iGT) > EpsTol & eps(:,iMut,iGT) < 1 - EpsTol );        
        TotCnt(iGT, iMut) = nnz( ~isnan(eps(:,iMut,iGT)) );
    end
end


close all;

%%% Specify the following dimensions:
fdim.spwa = 5; % subplotwidth in cm
fdim.spha = 5; % subplotheight in cm

fdim.nx = 3; % number of panels along the horizontal dimension
fdim.ny = 3; % number of panels along the vertical dimension

fdim.xma = [1.2 0.2]; % left right horizontal margin in cm
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

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 fdim.fw fdim.fh]);


fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );

% fdim

Ymax = 0.4;

if Ymax == 0
    Ymax = 1;
end

for iMut = 1:nMut
    ix = mod(iMut-1,3)+1;
    iy = floor((iMut-1)/3)+1;
        
    subplot('Position', [fdim.spxvec(ix) fdim.spyvec(iy) fdim.spwr fdim.sphr]),
    hold on, box on;
    set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');    
              
    bar(X,EpsCnt(:,iMut)./TotCnt(:,iMut), 0.9, 'FaceColor', 0.4*[1 1 1], 'EdgeColor', 'none', 'ShowBaseLine', 'off');
                
    set(gca, 'XLim', [0.25, Params.nGT+0.75], 'XTick', 1:Params.nGT);
    set(gca, 'YLim', Ymax*[-0.05 1.05]);
    xtickangle(90);
      
    if iy == 3
        set(gca, 'XTickLabel', GTID);
        if ix == 2
            xlabel('Generating topology', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
        end
    else
        set(gca, 'XTickLabel', {});
    end
    
    if ix > 1
        set(gca, 'YTickLabel', {});
    else
        if iy == 2
            ylabel('Fraction of sampled modules with 0 < \epsilon y_\mu < 1', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
        end
    end
        
    lim = axis;
    
    if ix == 1
        s = sprintf('\\delta^Ax_{ij} = %.2f\n\\delta^Bx_{kl} = %.2f', mutMat(iMut,1), mutMat(iMut,2)  );
    else
        s = sprintf('\\delta^Ax_{ij} = %.1f\n\\delta^Bx_{kl} = %.1f', mutMat(iMut,1), mutMat(iMut,2)  );
    end
                
    text( lim(2) - (lim(2)-lim(1))*0.02, lim(4) - (lim(4)-lim(3))*0.02,...
        s, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontName', 'Helvetica', 'FontSize', fdim.labelfs);
    
end

clear fdim iMut ix iy X Y edges lim iGT s;






%% Plotting the summary data across all strictly serial generating topologies (Figure 4B)

if strcmp(TopoRelation, 'serial')
    X = (1:3)';
    EpsCnt = nan( Params.nGT , nMut ); % # of cases with 0 < eps < 1
    TotCnt = nan( Params.nGT , nMut ); % total number of cases where epistasis could be computed
    
    for iMut = 1:nMut
        for iGT = 1:Params.nGT
            EpsCnt( iGT, iMut) = nnz( eps(:,iMut,iGT) > EpsTol & eps(:,iMut,iGT) < 1 - EpsTol );
            TotCnt(iGT, iMut) = nnz( ~isnan(eps(:,iMut,iGT)) );
        end
    end
    Y = mean( EpsCnt ./ TotCnt, 1 );
    Y = reshape( Y' , 3, 3);
    
    cc = [
        0, 114, 178;    % blue
        204, 121, 167;   % raddish purple
        213, 94, 0;     % vermillion
        86, 180, 233;   % sky blue
        230 159, 0;     % orange
        0, 158, 115;    % bluish green
        240, 228, 66   % yellow
        ]./256;
    
    %%% Specify the following dimensions:
    fdim.spwa = 5; % subplotwidth in cm
    fdim.spha = 5; % subplotheight in cm
    
    fdim.nx = 1; % number of panels along the horizontal dimension
    fdim.ny = 1; % number of panels along the vertical dimension
    
    fdim.xma = [1.8 0.2]; % left right horizontal margin in cm
    fdim.yma = [1 0.2]; % bottom top vertical margin cm
    
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
    
    subplot('Position', [fdim.spxvec(1) fdim.spyvec(1) fdim.spwr fdim.sphr]),
    hold on, box on;
    set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
    
    for i1 = 1:3
        plot(X, Y(:,i1), 'o-', 'Color', cc(i1,:), 'LineWidth', 2, ...
            'MarkerSize', 10, 'MarkerFaceColor', cc(i1,:), 'MarkerEdgeColor', 'w');
    end
    
    legend({'-,-', '+,+', '+,-'}, 'Location', 'NorthWest');
    
    set(gca, 'XLim', [0.9, 3.1], 'XTick', 1:3, 'XTickLabel', {'1%', '10%', '50%'});
    set(gca, 'YLim', [0, 0.2] + 0.2*0.05 * [-1 1], 'YTick', 0:0.05:Ymax);
    
    xlabel('Effect of mutations', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
    ylabel(sprintf('Fraction of sampled modules\nwith \\gamma < \\epsilon y_\\mu < 1-\\gamma'), 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
    
    clear fdim iMut i1 X Y cc Ymax;
end
