curr_dir = '/Users/skryazhi/epistasis_theory/data/Chassagnole_etal';
runId = '2019-12-20';

Glu = 'Low'; % 'Low' 'Med' or 'Inf';

IFINIT = false; % initialize models or load previously initialized?
% Should always re-initialized if something changes in the definition of
% any model


MetNames = struct(...
    'g6p', 'Glucose-6-Phosphate', ...
    'f6p', 'Fructose-6-Phosphate', ...
    'gap', 'Glyceraldehyde-3-Phosphate', ...
    'pep','Phosphoenol pyruvate', ...
    'ribu5p', 'Ribulose-5-phosphate', ...
    'glu', 'Extracellular Glucose',...
    'pyr', 'Pyruvate',...
    'g1p', 'Glucose-1-Phosphate',...
    'threepg', '3-Phosphoglycerate', ...
    'pgp', '1,3-diphosphosphoglycerate');

EnzNames = struct(... %%%% Upper glycolysis:
    'PTS', 'Phosphotransferase system', ...
    'PGI', 'Glucose-6-phosphate isomerase',... 
    'PFK', 'Phosphofructokinase', ...
    'ALDO', 'Aldolase', ... %%% Lower glycolysis:
    'TIS', 'Triosephosphate isomerase', ...
    'GAPDH', 'Glyceraldehyde-3-phosphate dehydrogenase', ...
    'PGK', 'Phosphoglycerate kinase', ...
    'PGluMu', 'Phosphoglycerate mutase', ...
    'ENO', 'Enolase', ... %%% PPP:
    'G6PDH', 'Glucose-6-phosphate dehydrogenase', ...
    'PGDH', '6-Phosphogluconate dehydrogenase', ...
    'R5PI', 'Ribose-phosphate isomerase', ...
    'Ru5P', 'Ribulose-phosphate epimerase', ...
    'TKa', 'Transketolase a', ...
    'TA', 'Transaldolase', ...
    'TKb', 'Transketolase b', ...
    'PK', 'Pyruvate kinase',...
    'PEPCxylase', 'PEP carboxylase');



%% Plot the distribution of FCCs in the full model
filename = sprintf('%s/model_FULL_%s.mat', curr_dir, Glu);
if exist(filename, 'file') == 0 || IFINIT
    fprintf('Initializing FULL model\n');
    [WT, ModelSpecs] = initialize_model(curr_dir, Glu, ModelSpecs, MetNames, EnzNames);
else
    load( filename );
end

WT.FluxDistr.ShortName = cell(length(WT.FluxDistr.Name),1);
ShortNameList = fieldnames( EnzNames );
for irxn = 1:length(ShortNameList)
    ShortName = ShortNameList{irxn};
    LongName = EnzNames.( ShortName );
    ix = find( strcmp(WT.FluxDistr.Name, LongName) );
    if ~isempty(ix) 
        WT.FluxDistr.ShortName{ix} = ShortName;
    end
end
clear ShortNameList ShortName LongName ix irxn;

bar(WT.FluxDistr.CtrlCoeff);
set(gca, 'XTick', 1:length(WT.FluxDistr.CtrlCoeff) ,'XTickLabel', WT.FluxDistr.ShortName);
xtickangle(30);
title(sprintf('FULL [Glu] = %s', Glu), 'FontName', 'Helvetica', 'FontSize', 16);
filename =  sprintf('%s/%s/figures [Glu] = %s/eps/FCC Distr.eps', curr_dir, runId, Glu);
saveas(gcf, filename, 'epsc');


%% Calculate epistasis cofficients
 
Pert = 0.9;
AbsTolEps = 1e-6;


% e1 = 'PGI';
% e2 = 'PGDH'; %'TIS'; % 'ALDO'; %'PFK';
% ModuleList = {'UGPPP' ; 'GPPP' ; 'FULL'};


e1 = 'PGDH';
e2 = 'PFK'; %  'ALDO'; %'TIS';
ModuleList = {'UGPPP' ; 'GPPP' ; 'FULL'};


% e1 = 'GAPDH';
% e2 = 'ENO'; % 'PGK', 'PGluMu'
% ModuleList = {'LG' ; 'GPPP' ; 'FULL'};

% e1 = 'PGDH'; %'ALDO'; %'PFK'; %'PGI'; 
% e2 = 'GAPDH';
% ModuleList = {'GPPP' ; 'FULL'};


nModule = length(ModuleList);
DeltaEpsMat = nan( nModule, 3);

% LG = lower glycolysis
% PPPSMALL = G6PDH and PGDH:    
% PPP = pentose phosphate pathway
% UGPPP = upper glycolysis and pentose phosphate pathway
% GPPP = upper and lower glycolysis and pentose phosphate pathway
% FULL = full model (w/o g1p and extreneous rxns)

for iModule = 1:nModule

    clear ModelSpecs WT mutList nMut FullEnzNames RxnNames rxnOutIdVec;
    
    ModelSpecs.Type = ModuleList{iModule};

    if IFINIT
        [WT, ModelSpecs] = initialize_model(curr_dir, Glu, ModelSpecs, MetNames, EnzNames);
    else
        load( sprintf('%s/model_%s_%s.mat', curr_dir, ModelSpecs.Type, Glu) );
    end
    
    RxnNames = get(WT.m.Reactions, 'Name');
    rxnOutIdVec = cellfun(@(EnzShortName) find(strcmp(RxnNames,EnzNames.(EnzShortName))) ,ModelSpecs.rxnOutEnzymes);

    %% Single mutants calculation:
    mutList = [ ...
        {{e1}, Pert};
        {{e2}, Pert};
        ];
    
    
    FullEnzNames = get_full_names(EnzNames, mutList(:,1));
    [ModelL, FluxDistrL] = getMutFlux(WT.m, [FullEnzNames, mutList(:,2)]);
    mutList = [mutList, ModelL, FluxDistrL];
    % Col 1 = enzyme short name list
    % Col 2 = perturbation list
    % Col 3 = perturbed model
    % Col 4 = perturbed steady state fluxes
    % Col 5 = perturbed steady state output flux
    % Col 6 = relative perturbation of output flux

    
    nMut = 2;
    
    %% Single mutants output:
    fprintf('WT:\nfOUT = %.3f \n', WT.OutFlux);
    
    for iMut = 1:nMut
        mutList{iMut, 5} = sum( ModelSpecs.rxnOutEnzymeCoeff .* mutList{iMut,4}.Flux(rxnOutIdVec) );
        mutList{iMut, 6} = (mutList{iMut, 5} - WT.OutFlux)/WT.OutFlux;
        
        DeltaEpsMat(iModule, iMut) = mutList{iMut, 6};
        
        fprintf('\\delta %s = %.3g\n\f=> fOUT = %.3f\n\t=> \\delta f = %.2e\n',...
            mutList{iMut, 1}{1}, mutList{iMut, 2}, ...
            mutList{iMut, 5},...
            mutList{iMut, 6} );
    end
    
    
    %% Double mutants
    
    mut2List = cell(1, 8);
    
    mut2List{1,1} = [mutList{1,1}; mutList{2,1}];
    mut2List{1,2} = [mutList{1,2}; mutList{2,2}];
    mut2List{1,7} = [1; 2];
    
    FullEnzNames = get_full_names(EnzNames, mut2List(:,1));
    [ModelL, FluxDistrL] = getMutFlux(WT.m, [FullEnzNames, mut2List(:,2)]);
    mut2List(:,3:4) = [ModelL, FluxDistrL];
    % Col 1 = enzyme short name list
    % Col 2 = perturbation list
    % Col 3 = perturbed model
    % Col 4 = perturbed steady state fluxes
    % Col 5 = perturbed steady state output flux
    % Col 6 = relative perturbation of output flux
    % Col 7 = indices of perturbed enzymes in mutList
    % Col 8 = epsilon
    
    
    %%% Calculate epistasis coeff:
    mut2List{1, 5} = sum( ModelSpecs.rxnOutEnzymeCoeff .* mut2List{1,4}.Flux(rxnOutIdVec) );
    mut2List{1, 6} = (mut2List{1, 5} - WT.OutFlux)/WT.OutFlux;
    
    fprintf('\\delta %s = %.3g, \\delta %s = %.3g\n\t=> fOUT = %.3f\n\t=> \\delta f = %.4g\n',...
        mut2List{1, 1}{1}, mut2List{1, 2}(1), ...
        mut2List{1, 1}{2}, mut2List{1, 2}(2), ...
        mut2List{1, 5},...
        mut2List{1, 6} );
    
    iMut1 = mut2List{1,7}(1);
    iMut2 = mut2List{1,7}(2);
    
    delta1 = mutList{iMut1,6};
    delta2 = mutList{iMut2,6};
    
    eps = mut2List{1, 6} - delta1 - delta2;
    if abs(delta1) < AbsTolEps || abs(delta2) < AbsTolEps || abs(mut2List{1, 6}) < AbsTolEps % || ...
            % abs(-1-delta1) < AbsTolEps || abs(-1-delta2) < AbsTolEps || abs(-1-mut2List{1, 6}) < AbsTolEps
        eps = nan;
    elseif abs(eps) < AbsTolEps
        eps = 0;
    else
        eps = eps/delta1/delta2;
    end
    mut2List{1,8} = eps;
    
    DeltaEpsMat(iModule, 3) = eps;
    
    fprintf('\t eps f = %.3g\n', eps);
    
    clear iMut1 iMut2 rxnId delta1 delta2 eps;
    
    %% Saving
    filename =  sprintf('%s/%s/data [Glu] = %s/%s %s Pert = %.1f %s.mat', curr_dir, runId, Glu, e1, e2, Pert, ModelSpecs.Type);
    save(filename, 'ModelSpecs', 'WT', 'mutList', 'mut2List');
end

filename =  sprintf('%s/%s/data [Glu] = %s/%s %s Pert = %.1f.mat', curr_dir, runId, Glu, e1, e2, Pert);
save(filename, 'ModelSpecs', 'WT', 'DeltaEpsMat');


%% Visualize epistasis
clf;

%%% Specify the following dimensions:
fdim.spwa = 8; % subplotwidth in cm
fdim.spha = 8; % subplotheight in cm

fdim.nx = 1; % number of panels along the horizontal dimension
fdim.ny = 2; % number of panels along the vertical dimension

fdim.xma = [1.5 0.5]; % left right horizontal margin in cm
fdim.yma = [1.3 1]; % bottom top vertical margin cm

fdim.dxa = 0.3; % horizontal distance between panels in cm
fdim.dya = 0.6; % vertical distance between panels in cm

fdim.tickfs = 10;
fdim.labelfs = 14;

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

cc = [
    0, 114, 178;    % blue
    213, 94, 0;     % vermillion
    86, 180, 233;   % sky blue
    230 159, 0;     % orange
    204, 121, 167;   % raddish purple
    0, 158, 115;    % bluish green
    240, 228, 66   % yellow
    ]./256;

fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );

%%% Plot deltas
subplot('Position', [fdim.spxvec(1) fdim.spyvec(1) fdim.spwr fdim.sphr]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');

bar( (1:nModule)'-0.2, DeltaEpsMat(:,1), 0.4, 'FaceColor', 0.5*[1 1 1], 'EdgeColor', 'none');
bar( (1:nModule)'+0.2, DeltaEpsMat(:,2), 0.4, 'FaceColor', 0.8*[1 1 1], 'EdgeColor', 'none');

set(gca, 'YScale', 'log', 'YLim', [-0.1, -1e-4]);
set(gca, 'XLim', [0.5, nModule+0.5], 'XTick', 1:nModule);
set(gca, 'XTickLabel', {});
ylabel('Effects of mutations on flux', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs);
title(sprintf('%s, %s', e1, e2), 'FontName', 'Helvetica', 'FontSize', 16);


%%% Plot eps
subplot('Position', [fdim.spxvec(1) fdim.spyvec(2) fdim.spwr fdim.sphr]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');

plot( [0.5 nModule+0.5], zeros(1,2), '-', 'LineWidth', 1, 'Color', 0.4*[1 1 1]);

plot( (1:nModule)', DeltaEpsMat(:,3), 'ok-', 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

% set(gca, 'YScale', 'log', 'YLim', [-1, -1e-4]);
set(gca, 'XLim', [0.5, nModule+0.5], 'XTick', 1:nModule, 'XTickLabel', ModuleList);
xlabel('Module', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs);
ylabel('Epistasis for flux', 'FontName', 'Helvetica', 'FontSize', fdim.labelfs);

filename =  sprintf('%s/%s/figures [Glu] = %s/eps/%s %s Pert = %.1f.eps', curr_dir, runId, Glu, e1, e2, Pert);
saveas(gcf, filename, 'epsc');


% fdim.ms = 1;
% 
% nMut = size(epsMat,1);
% hpos = [];
% hneg = [];
% 
% for iMut1 = 1:nMut
%     xpos = iMut1;
%     ypos = nMut-(iMut1-1);
%     
%     delta1 = mutList{iMut1,6};
% 
%     text(xpos, ypos, sprintf('%.1e', delta1),...
%         'HorizontalAlignment','center', 'VerticalAlignment', 'middle',...
%         'FontName', 'Helvetica', 'FontSize', fdim.labelfs); 
% end
% 
% 
% for iMut1 = 1:nMut-1
%     for iMut2 = iMut1+1:nMut
%                 
%         if isnan(epsMat(iMut1,iMut2))
%             c = 0.4 * [1 1 1];
%         elseif epsMat(iMut1,iMut2) < 0
%             c = cc(2,:);
%             % hneg(end+1) = plot( iMut2, nMut-(iMut1-1), 's', 'MarkerSize', 22, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', c);
%         elseif epsMat(iMut1,iMut2) > 0
%             c = cc(6,:);
%             % hpos(end+1) = plot( iMut2, nMut-(iMut1-1), 's', 'MarkerSize', 22, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', c);
%         else 
%             c = cc(7,:);
%         end
%         xpos = iMut2;
%         ypos = nMut-(iMut1-1);
%         
%         rectangle('Position', [xpos-fdim.ms/2, ypos-fdim.ms/2, fdim.ms, fdim.ms],...
%             'EdgeColor', 'none', 'FaceColor', c);
%         text(xpos, ypos, sprintf('%.1f', epsMat(iMut1,iMut2)),...
%             'HorizontalAlignment','center', 'VerticalAlignment', 'middle',...
%             'FontName', 'Helvetica', 'FontSize', fdim.labelfs);
%         
%         delta1 = mutList{iMut1,6};
%         delta2 = mutList{iMut2,6};
% 
%         xpos = iMut1;
%         ypos = nMut-(iMut2-1);
%                 
%         rectangle('Position', [xpos-fdim.ms/2, ypos-fdim.ms/2, fdim.ms, fdim.ms],...
%             'EdgeColor', 'none', 'FaceColor', c);
%         text(xpos, ypos, sprintf('%.1e', epsMat(iMut1,iMut2)*delta1*delta2),...
%             'HorizontalAlignment','center', 'VerticalAlignment', 'middle',...
%             'FontName', 'Helvetica', 'FontSize', fdim.labelfs);
%     end
% end
% 
% for pos = 1.5:1:nMut-0.5
%     plot( [0.5 nMut+0.5], pos*[1 1], '-', 'Color', 0.8*[1 1 1]);
%     plot( pos*[1 1], [0.5 nMut+0.5], '-', 'Color', 0.8*[1 1 1]);
% end
% 
% MutEnzList = cellfun(@(x) x{1}, mutList(:,1), 'UniformOutput', false);
% 
% set(gca, 'XLim', [0.5 nMut+0.5], 'YLim', [0.5 nMut+0.5], ... %'XGrid', 'on', 'YGrid', 'on', ...
%     'XTick', 1:nMut, 'YTick', 1:nMut, 'XTickLabel', MutEnzList, 'YTickLabel', flipud(MutEnzList))
% xtickangle(90);
% % legend([hpos(1), hneg(1)], 'Pos', 'Neg', 'Location', 'SouthWest');
% 
% clear fdim cc hpos hneg iMut1 iMut2 c pos;
% 
% %% Saving
% filename =  sprintf('%s/%s/figures/eps/%s/%s', curr_dir, runId, Glu, ModelSpecs.Type);
% 
% for met = fieldnames(ModelSpecs.ExtMetConc)'
%     filename = sprintf('%s, [%s] = %.2f', filename, met{1}, ModelSpecs.ExtMetConc.(met{1}));
% end
% 
% filename = sprintf('%s, Pert = %.1f', filename, Pert);
% 
% if ModelSpecs.IFREMOVE6PGINH
%     filename = sprintf('%s, no6pginh.eps', filename);
% else
%     filename = sprintf('%s.eps', filename);
% end
% 
% saveas(gcf, filename, 'epsc');


