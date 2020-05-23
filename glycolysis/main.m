curr_dir = '/Users/skryazhi/epistasis_theory/data/Chassagnole_etal';
runId = '2020-05-23';

Glu = 'Low'; % 'Low' 'Med' or 'Inf';

IFINIT = false; % initialize models or load previously initialized?
% Should always re-initialized if something changes in the definition of
% any model

ModuleList = {'FULL', 'UGPP', 'GPP', 'LG'};

data_dir = sprintf('%s/%s/data [Glu] = %s', curr_dir, runId, Glu); 
fig_dir = sprintf('%s/%s/figures [Glu] = %s', curr_dir, runId, Glu); 

modelfile = sprintf('%s/BIOMD0000000051.xml', curr_dir);

%% Initializing

for iMod = 1:length(ModuleList)
    ModelName = ModuleList{iMod};
    filename = sprintf('%s/model_%s_%s.mat', data_dir, ModelName, Glu);
    
    if exist(filename, 'file') == 0 || IFINIT
        fprintf('Initializing %s model\n', ModelName);
        [WT, ModelSpecs] = initialize_model(modelfile, data_dir, ModelName, Glu);
    end
end




%% Plot Figure 3B: all FCC for the FULL model

close all;

Glu = 'Low';
ModelName = 'FULL';

filename = sprintf('%s/model_%s_%s.mat', data_dir, ModelName, Glu);
load( filename );

RxnNames = cellfun(@(rxn) get_short_name( ModelSpecs.EnzNames, rxn), WT.FluxDistr.Name, 'UniformOutput', false);
FCC = WT.FluxDistr.FCC;

[FCC, ix] = sort( FCC, 'descend');
RxnNames = RxnNames( ix );
nRxn = length(RxnNames);
        
clear filename ix;


%%% Specify the following dimensions:
fdim.spwa = 5; % subplotwidth in cm
fdim.spha = 4; % subplotheight in cm

fdim.nx = 1; % number of panels along the horizontal dimension
fdim.ny = 1; % number of panels along the vertical dimension

fdim.xma = [1 0.2]; % left right horizontal margin in cm
fdim.yma = [2 0.2]; % bottom top vertical margin cm

fdim.dxa = 0.3; % horizontal distance between panels in cm
fdim.dya = 0; % vertical distance between panels in cm

fdim.tickfs = 8;
fdim.labelfs = 12;


%%% --- These will be computed automatically:
fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
fdim.fh = sum(fdim.spha) + sum(fdim.dya) + sum(fdim.yma);

fdim.spwr = fdim.spwa / fdim.fw;
fdim.sphr = fdim.spha / fdim.fh;
fdim.xmr = fdim.xma / fdim.fw;
fdim.ymr = fdim.yma / fdim.fh;
fdim.dxr = fdim.dxa / fdim.fw;
fdim.dyr = fdim.dya / fdim.fh;

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 fdim.fw fdim.fh]);

fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
tmp = cumsum(fliplr([fdim.sphr 0]));
tmp = tmp(1:end-1) + cumsum(fliplr([fdim.dyr 0]));
fdim.spyvec = fdim.ymr(1) + fliplr(tmp);

subplot('Position', [fdim.spxvec(1) fdim.spyvec(1) fdim.spwr fdim.sphr(1)]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');

cc = [
    0, 114, 178;    % blue
    213, 94, 0;     % vermillion
    204, 121, 167;   % raddish purple
    86, 180, 233;   % sky blue
    230 159, 0;     % orange
    0, 158, 115;    % bluish green
    240, 228, 66   % yellow
    ]./256;


bar((1:nRxn)', FCC, 'FaceColor', 'k', 'EdgeColor', 'none');

set(gca, 'XLim', [0.3 nRxn+0.7], 'XTick', 1:nRxn);
set(gca, 'XTickLabel', RxnNames );
xtickangle(90);
set(gca, 'YLim', [-0.05, 0.35], 'YTick', 0:0.1:0.3, 'YGrid', 'off');
ylabel('FCC', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs); %'Position', [0.53, 0.5]
xlabel('Reaction', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs, 'Position', [10, -0.18]);
clear  cc;

text(10, 0.34, 'Module FULL', 'FontName', 'Helvetica', 'FontSize', 12,...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

filename =  sprintf('%s/eps/FCC FUL [Glu]=%s.eps', fig_dir, Glu);
% saveas(gcf, filename, 'epsc');






%% Plot Figure 3C (FCC and EPI for selected rxns)

Glu = 'Low';
RxnList = {'PGI', 'PGDH', 'PFK'};
ModuleList = {'UGPP'; 'GPP'; 'FULL'};

nMod = length(ModuleList);
nRxn = length(RxnList);
nEps = nRxn * (nRxn - 1)/2;

FCC = nan( nMod, nRxn );
EPS = nan( nMod, nEps );

for iMod = 1:nMod
    ModelName = ModuleList{iMod};
    filename = sprintf('%s/model_%s_%s.mat', data_dir, ModelName, Glu);
    load( filename );
        
    for irxn = 1:length(RxnList)
        ix = find( strcmp( WT.FluxDistr.Name, get_full_name( ModelSpecs.EnzNames , RxnList{ irxn }) ) );
        FCC( iMod, irxn ) = WT.FluxDistr.FCC( ix );
    end
    
    ieps = 1;
    for irxn1 = 1:length(RxnList)-1
        ix1 = find( strcmp( WT.FluxDistr.Name, get_full_name( ModelSpecs.EnzNames , RxnList{ irxn1 }) ) );
        
        for irxn2 = irxn1+1:length(RxnList)
            ix2 = find( strcmp( WT.FluxDistr.Name, get_full_name( ModelSpecs.EnzNames , RxnList{ irxn2 }) ) );
            EPS( iMod, ieps ) = WT.FluxDistr.FIC( ix1, ix2 ) / 2 / FCC( iMod, irxn1 ) / FCC( iMod, irxn2 );
            ieps = ieps + 1;
        end 
    end    
end
clear iMod irxn ix irxn1 irx1 irxn2 ix2 ieps WT ModelSpecs;

figure;

%%% Specify the following dimensions:
fdim.spwa = 5; % subplotwidth in cm
fdim.spha = [3.5 2 2 2]; % subplotheight in cm

fdim.nx = 1; % number of panels along the horizontal dimension
fdim.ny = nEps+1; % number of panels along the vertical dimension

fdim.xma = [1.2 0.2]; % left right horizontal margin in cm
fdim.yma = [1.1 0.2]; % bottom top vertical margin cm

fdim.dxa = 0.3; % horizontal distance between panels in cm
fdim.dya = [1.5 0.2 0.2]; % vertical distance between panels in cm

fdim.tickfs = 8;
fdim.labelfs = 12;

%%% These will be computed automatically:
fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
fdim.fh = sum(fdim.spha) + sum(fdim.dya) + sum(fdim.yma);

fdim.spwr = fdim.spwa / fdim.fw;
fdim.sphr = fdim.spha / fdim.fh;
fdim.xmr = fdim.xma / fdim.fw;
fdim.ymr = fdim.yma / fdim.fh;
fdim.dxr = fdim.dxa / fdim.fw;
fdim.dyr = fdim.dya / fdim.fh;

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 fdim.fw fdim.fh]);

fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
tmp = cumsum(fliplr([fdim.sphr 0]));
tmp = tmp(1:end-1) + cumsum(fliplr([fdim.dyr 0]));
fdim.spyvec = fdim.ymr(1) + fliplr(tmp);

fdim.minx = 1 - 0.2;
fdim.maxx = nMod + 1 - fdim.minx;

subplot('Position', [fdim.spxvec(1) fdim.spyvec(1) fdim.spwr fdim.sphr(1)]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');

%%% Plot FCC:
% cc = repmat(linspace(0.3, 0.7, nRxn)', 1, 3);
cc = [
    0, 114, 178;    % blue
    213, 94, 0;     % vermillion
    204, 121, 167;   % raddish purple
    86, 180, 233;   % sky blue
    230 159, 0;     % orange
    0, 158, 115;    % bluish green
    240, 228, 66   % yellow
    ]./256;

for iRxn = 1:nRxn
    plot( (1:nMod)' , FCC(:,iRxn), '-', 'LineWidth', 2, 'Color', cc(iRxn, :));
    plot( (1:nMod)' , FCC(:,iRxn), 'o', 'LineWidth', 1, ...
        'MarkerFaceColor', cc(iRxn, :), 'MarkerEdgeColor', 'w', 'MarkerSize', 10);
end

set(gca, 'XLim', [fdim.minx fdim.maxx], 'XTick', 1:nMod);
set(gca, 'XTickLabel', ModuleList );
set(gca, 'YLim', [-0.1, 1.1], 'YTick', [0 1], 'YGrid', 'on');
ylabel('FCC', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs,...
    'Position', [0.53, 0.5]);
xlabel('Module', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs,...
    'Position', [2, -0.27]);
text(1.08, 0.17, RxnList{1}, 'Color', cc(1,:), 'VerticalAlignment', 'bottom');
text(1.08, 0.08, RxnList{2}, 'Color', cc(2,:), 'VerticalAlignment', 'top');
text(1.08, 0.83, RxnList{3}, 'Color', cc(3,:), 'VerticalAlignment', 'top');
clear iRxn cc;




%%% Plot eps
for iEps = 1:nEps
    subplot('Position', [fdim.spxvec(1) fdim.spyvec(1+iEps) fdim.spwr fdim.sphr(iEps+1)]),
    hold on, box on;
    set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
        
    plot( (1:nMod)' , EPS(:,iEps), '-', 'LineWidth', 1, 'Color', 'k');
    plot( (1:nMod)' , EPS(:,iEps), 'o', 'LineWidth', 1, ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerSize', 10);

    plot([fdim.minx, fdim.maxx], [0 0], 'Color', 0.8*[1 1 1]);

    set(gca, 'XLim', [fdim.minx fdim.maxx], 'XTick', 1:nMod, 'XTickLabel', {});
    
    set(gca, 'YGrid', 'off');
    
    
    if iEps == 1
        set(gca, 'YLim', [-35, 5], 'YTick', [-30 0]);
        text(0.85, -34.5, sprintf('%s-%s', RxnList{1}, RxnList{2}),...
            'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    elseif iEps == 2
        plot([fdim.minx, fdim.maxx], [1 1], 'Color', 0.8*[1 1 1]);
        set(gca, 'YLim', [-0.5, 1.2], 'YTick', [0 1]);
        ylabel('Epistasis for flux', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs, ...
            'Position', [0.53, 0.35]);
        text(0.85, -0.47, sprintf('%s-%s', RxnList{1}, RxnList{3}),...
            'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    elseif iEps == 3
        plot([fdim.minx, fdim.maxx], [1 1], 'Color', 0.8*[1 1 1]);
        set(gca, 'YLim', [-0.5, 1.2], 'YTick', [0 1]);
        set(gca, 'XTickLabel', ModuleList );
        xlabel('Module', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs,...
            'Position', [2, -0.85]);
        text(0.85, -0.47, sprintf('%s-%s', RxnList{2}, RxnList{3}),...
            'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    end    
end


filename =  sprintf('%s/eps/FCC EPI [Glu]=%s.eps', fig_dir, Glu);
% saveas(gcf, filename, 'epsc');




%% Plot the figure for the supplementary material (all Rxns)

figure;

Glu = 'Low'; % 'Inf';
AbsTol = 1e-4;
ModuleList = {'LG'; 'UGPP'; 'GPP'; 'FULL'};
nMod = length(ModuleList);

%%% Specify RxnList explicitly:
RxnList = {...
    'PGI','Glucose-6-phosphate isomerase';
    'PGDH','6-Phosphogluconate dehydrogenase'
    'PFK','Phosphofructokinase';
    'ALDO','Aldolase';
    'TIS','Triosephosphate isomerase';    
    'GAPDH','Glyceraldehyde-3-phosphate dehydrogenase';
    'PGK','Phosphoglycerate kinase';
    'PGluMu','Phosphoglycerate mutase';
    'ENO','Enolase';
    'PK','Pyruvate kinase';
    'PTS','Phosphotransferase system'};

%%% Or load it from the FULL model:
% ModelName = 'FULL';
% filename = sprintf('%s/model_%s_%s.mat', data_dir, ModelName, Glu);
% load( filename );
% TF = WT.FluxDistr.FCC > AbsTol;
% clear RxnList;
% RxnList(:,2) = WT.FluxDistr.Name( TF );
% RxnList(:,1) = map_full2short( ModelSpecs.EnzNames, RxnList(:,2) );
% clear TF;

nRxn = size(RxnList, 1);
FCC = nan( nMod, nRxn );
EPS = nan( nRxn, nRxn, nMod);

for iMod = 1:nMod
    ModelName = ModuleList{iMod};
    filename = sprintf('%s/model_%s_%s.mat', data_dir, ModelName, Glu);    
    load( filename );
    
    for irxn = 1:length(WT.FluxDistr.Name)
        ix = find(strcmp( RxnList(:,2) , WT.FluxDistr.Name{irxn}));
        if isempty(ix)
            continue;
        end
        FCC(iMod, ix) = WT.FluxDistr.FCC(irxn);
    end
    
    for irxn1 = 1:length(WT.FluxDistr.Name)-1
        ix1 = find(strcmp( RxnList(:,2) , WT.FluxDistr.Name{irxn1}));
        if isempty(ix1)
            continue;
        end
        
        for irxn2 = irxn1+1:length(WT.FluxDistr.Name)
            ix2 = find(strcmp( RxnList(:,2) , WT.FluxDistr.Name{irxn2})); 

            if isempty(ix2)
                continue;
            end

            EPS( ix1, ix2, iMod ) = WT.FluxDistr.FIC( irxn1, irxn2 ) / 2 / FCC( iMod, ix1 ) / FCC( iMod, ix2 );
            EPS( ix2, ix1, iMod ) = EPS( ix1, ix2, iMod );
        end
    end    
end
clear iMod ModelName filename irxn ix irxn1 ix1 irxn2 ix2;


clf;


%%% Specify the following dimensions:
fdim.spwa = 3; % subplotwidth in cm
fdim.spha = 3; % subplotheight in cm

fdim.nx = nRxn; % number of panels along the horizontal dimension
fdim.ny = nRxn; % number of panels along the vertical dimension

fdim.xma = [1 2]; % left right horizontal margin in cm
fdim.yma = [1.1 0.8]; % bottom top vertical margin cm

fdim.dxa = 0.5; % horizontal distance between panels in cm
fdim.dya = 0.5; % vertical distance between panels in cm

fdim.tickfs = 8;
fdim.labelfs = 12;

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

fdim.minx = 1 - 0.2;
fdim.maxx = nMod + 1 - fdim.minx;

%%% First plot FCC for each rxn
for iRxn = 1:nRxn
    subplot('Position', [fdim.spxvec(iRxn) fdim.spyvec(1) fdim.spwr fdim.sphr(1)]),
    hold on, box on;
    set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
    
    ix = find(~isnan(FCC(:,iRxn)));
    X = ix ;
    Y = FCC( ix, iRxn );
    
    plot( X , Y, '-', 'LineWidth', 2, 'Color', 'k');
    plot( X , Y, 'o', 'LineWidth', 1, ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerSize', 10);
    
    set(gca, 'XLim', [fdim.minx fdim.maxx], 'XTick', 1:nMod);
    set(gca, 'XTickLabel', ModuleList );
    set(gca, 'YLim', [-0.1, 1.1], 'YTick', [0 1], 'YGrid', 'on');
    title(RxnList{iRxn}, 'FontName', 'Helvetica', 'FontSize', fdim.labelfs);

    if iRxn == 1
     ylabel('FCC', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs,...
         'Position', [0.53, 0.5]);
    end    
end
clear iRxn ix X Y;


%%% Second plot EPS for each rxn
for iRxn1 = 2:nRxn
    for iRxn2 = 1:iRxn1-1
        
        subplot('Position', [fdim.spxvec(iRxn1) fdim.spyvec(iRxn2+1) fdim.spwr fdim.sphr(1)]),
        hold on, box on;
        set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
        
        ix = find(~isnan( EPS(iRxn1,iRxn2,:)) );
        X = ix ;
        Y = squeeze(EPS(iRxn1,iRxn2,ix));
        
        plot( X , Y, '-', 'LineWidth', 2, 'Color', 'k');
        plot( X , Y, 'o', 'LineWidth', 1, ...
            'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerSize', 10);
        
        set(gca, 'XLim', [fdim.minx fdim.maxx], 'XTick', 1:nMod);
        set(gca, 'XTickLabel', ModuleList );
                
        m = 0.05;
        
        ylim(1) = min(0, floor( min(Y) ));
        ylim(2) = max(0, ceil( max(Y) ) );
        if ylim(1) == ylim(2)
            ylim(1) = ylim(1) - 1;
            ylim(2) = ylim(2) + 1;
        end
        l = ylim(2) - ylim(1);
        
        ytick = ylim;
        
        if ylim(1) < 0 && ylim(2) > 0
            ytick = [ytick(ytick < 0) , 0 , ytick(ytick>0)];
        end
        if ylim(1) < 1 && ylim(2) > 1
            ytick = [ytick(ytick < 1) , 1 , ytick(ytick>1)];
        end
        set(gca, 'YLim', [ylim(1)-l*m ,  ylim(2)+l*m], 'YTick', ytick, 'YGrid', 'on');
        
        if iRxn2 == iRxn1-1
            ylabel('Epistasis', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs,...
                'Position', [0.5, mean(ylim)]);
        end
        
        if iRxn1 == nRxn
            text(4.6, mean(ylim), RxnList{iRxn2}, 'FontName', 'Helvetica', 'FontSize', fdim.labelfs,...
                'FontWeight', 'bold', 'HorizontalAlignment', 'left');%  'Rotation', 270, );
        end
    end
    % title(RxnList{iRxn}, 'FontName', 'Helvetica', 'FontSize', fdim.labelfs);
    %     ylabel('FCC', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs,...
    %         'Position', [0.53, 0.5]);
    %     xlabel('Module', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs,...
    %         'Position', [2, -0.27]);
    
end
clear iRxn ix X Y ylim tmp l m ytick;
filename =  sprintf('%s/eps/FCC EPI [Glu]=%s all.eps', fig_dir, Glu);
% saveas(gcf, filename, 'epsc');







