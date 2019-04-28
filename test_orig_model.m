% This script first importa the original time-dependent model by
% Chassingole and tests that it reproduces the original findings, figures
% 4, 5, 6 in the paper.
%
% Then it modifies the original time-dependent model by Chassingole into a
% steady-state model. In particular, (1) I remove the (time-dependent)
% rules for ATP, ADP, AMP, NADH, NAD, NADPH, and NADP (given in Table VI)
% are removed. Now these metabolite concentrations are parameters with
% their respective time-independent steady-state concentrations given in
% Table V. (2) I also set the external glucose initial condition to 0.0556
% mM (Table V) and constant. The check that everything works properly, I
% calculate the steady-state concentrations of all metabolites and compare
% to Table V.


%% Load the model:
% currdir = '/Users/skryazhi/epistasis_theory/data/Chassagnole_etal';
% filename = sprintf('%s/BIOMD0000000051.xml', currdir);
% 
% OrigModelObj = sbmlimport(filename);
% 
% csObj = getconfigset(OrigModelObj,'active');
% set(csObj,'Stoptime',40);

% [t,x,names] = sbiosimulate(OrigModelObj);


% close all;
% 
% %% PLOT Figure 4:
% 
% ix = nan(7,1);
% ix(1) = find(strcmp('catp', names)); 
% ix(2) = find(strcmp('cadp', names));
% ix(3) = find(strcmp('camp', names));
% ix(4) = find(strcmp('cnad', names)); 
% ix(5) = find(strcmp('cnadh', names)); 
% ix(6) = find(strcmp('cnadp', names));
% ix(7) = find(strcmp('cnadph', names));
% 
% figure
% 
% %%% Specify the following dimensions:
% fdim.spwa = 6; % subplotwidth in cm
% fdim.spha = 5; % subplotheight in cm
% 
% fdim.nx = 1; % number of panels along the horizontal dimension
% fdim.ny = 3; % number of panels along the vertical dimension
% 
% fdim.xma = [1.5 0.5]; % left right horizontal margin in cm
% fdim.yma = [1.3 0.5]; % bottom top vertical margin cm
% 
% fdim.dxa = 0; % horizontal distance between panels in cm
% fdim.dya = 0; % vertical distance between panels in cm
% 
% fdim.tickfs = 8;
% fdim.labelfs = 10;
% 
% %%% These will be computed automatically:
% fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
% fdim.fh = fdim.spha * fdim.ny + fdim.dya * (fdim.ny - 1) + sum(fdim.yma);
% 
% fdim.spwr = fdim.spwa / fdim.fw;
% fdim.sphr = fdim.spha / fdim.fh;
% fdim.xmr = fdim.xma / fdim.fw;
% fdim.ymr = fdim.yma / fdim.fh;
% fdim.dxr = fdim.dxa / fdim.fw;
% fdim.dyr = fdim.dya / fdim.fh;
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', [0 0 fdim.fw fdim.fh]);
% 
% cc = [
%     0, 114, 178;    % blue
%     213, 94, 0;     % vermillion
%     86, 180, 233;   % sky blue
%     230 159, 0;     % orange
%     204, 121, 167;   % raddish purple
%     0, 158, 115;    % bluish green
%     240, 228, 66   % yellow
%     ]./256;
% 
% fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
% fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );
% 
% fdim
% 
% %%% PANEL (a)
% subplot('Position', [fdim.spxvec(1) fdim.spyvec(1) fdim.spwr fdim.sphr]),
% hold on, box on;
% set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
% % ATP
% plot( t, x(:,ix(1)), 'k-', 'LineWidth', 2, 'Marker', 'none');
% % ADP
% plot( t, x(:,ix(2)), 'k:', 'LineWidth', 2, 'Marker', 'none');
% % AMP
% plot( t, x(:,ix(3)), 'k--', 'LineWidth', 2, 'Marker', 'none');
% set(gca, 'YLim', [0 5], 'XLim', [-10 40], 'YGrid', 'on', 'XGrid', 'on', 'XTickLabel', {});
% legend({'ATP', 'ADP', 'AMP'}, 'Location', 'North');
% text(39, 4.9, '(a)', 'FontSize', fdim.labelfs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
% 
% 
% %%% PANEL (b)
% subplot('Position', [fdim.spxvec(1) fdim.spyvec(2) fdim.spwr fdim.sphr]),
% hold on, box on;
% set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
% % NAD
% plot( t, x(:,ix(4)), 'k-', 'LineWidth', 2, 'Marker', 'none');
% % NADH
% plot( t, x(:,ix(5)), 'k:', 'LineWidth', 2, 'Marker', 'none');
% set(gca, 'YLim', [0 2], 'YTick', 0:0.5:1.5, 'XLim', [-10 40], 'YGrid', 'on', 'XGrid', 'on', 'XTickLabel', {});
% legend({'NAD', 'NADH'}, 'Location', 'North');
% text(39, 1.95, '(b)', 'FontSize', fdim.labelfs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
% ylabel('concentration, mM', 'FontSize', fdim.labelfs);
% 
% %%% PANEL (c)
% subplot('Position', [fdim.spxvec(1) fdim.spyvec(3) fdim.spwr fdim.sphr]),
% hold on, box on;
% set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
% % NADP
% plot( t, x(:,ix(6)), 'k-', 'LineWidth', 2, 'Marker', 'none');
% % NADPH
% plot( t, x(:,ix(7)), 'k:', 'LineWidth', 2, 'Marker', 'none');
% set(gca, 'YLim', [0.04 0.24], 'YTick', 0.04:0.04:0.20, 'XLim', [-10 40], 'YGrid', 'on', 'XGrid', 'on');
% legend({'NADP', 'NADPH'}, 'Location', 'North');
% text(39, 0.235, '(c)', 'FontSize', fdim.labelfs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
% xlabel('time, s', 'FontSize', fdim.labelfs);
% 
% clear fdim cc;
% 
% 
% 
% 
% 
% %% PLOT Figure 5:
% 
% ix = nan(9,1);
% ixnames = {'extracellular.Extracellular Glucose', 'cytosol.Fructose-1,6-bisphosphate', 'cytosol.Glucose-1-Phosphate',...
%     'cytosol.Glucose-6-Phosphate', 'cytosol.Phosphoenol pyruvate', 'cytosol.Pyruvate',...
%     'cytosol.Fructose-6-Phosphate', 'cytosol.Glyceraldehyde-3-Phosphate','cytosol.6-Phosphogluconate'};
% 
% for i1 = 1:9
%     ix(i1) = find(strcmp(ixnames{i1}, names));
% end
% 
% figure
% 
% %%% Specify the following dimensions:
% fdim.spwa = 6; % subplotwidth in cm
% fdim.spha = 5; % subplotheight in cm
% 
% fdim.nx = 3; % number of panels along the horizontal dimension
% fdim.ny = 3; % number of panels along the vertical dimension
% 
% fdim.xma = [1.5 0.5]; % left right horizontal margin in cm
% fdim.yma = [1.3 0.5]; % bottom top vertical margin cm
% 
% fdim.dxa = 0; % horizontal distance between panels in cm
% fdim.dya = 0; % vertical distance between panels in cm
% 
% fdim.tickfs = 8;
% fdim.labelfs = 10;
% 
% %%% These will be computed automatically:
% fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
% fdim.fh = fdim.spha * fdim.ny + fdim.dya * (fdim.ny - 1) + sum(fdim.yma);
% 
% fdim.spwr = fdim.spwa / fdim.fw;
% fdim.sphr = fdim.spha / fdim.fh;
% fdim.xmr = fdim.xma / fdim.fw;
% fdim.ymr = fdim.yma / fdim.fh;
% fdim.dxr = fdim.dxa / fdim.fw;
% fdim.dyr = fdim.dya / fdim.fh;
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', [0 0 fdim.fw fdim.fh]);
% 
% cc = [
%     0, 114, 178;    % blue
%     213, 94, 0;     % vermillion
%     86, 180, 233;   % sky blue
%     230 159, 0;     % orange
%     204, 121, 167;   % raddish purple
%     0, 158, 115;    % bluish green
%     240, 228, 66   % yellow
%     ]./256;
% 
% fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
% fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );
% 
% fdim
% 
% for i1 = 1:9
%     ixx = mod(i1-1,3) + 1;
%     ixy = floor((i1-1)/3) + 1;
%     
%     subplot('Position', [fdim.spxvec(ixx) fdim.spyvec(ixy) fdim.spwr fdim.sphr]),
%     hold on, box on;
%     set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
%     plot( t, x(:,ix(i1)), 'k-', 'LineWidth', 2, 'Marker', 'none');
%     set(gca,  'XLim', [-10 40], 'XTick', 0:10:40, 'YGrid', 'on', 'XGrid', 'on');
%     
%     if ixy == 1
%         set(gca, 'YLim', [0, 3.5], 'YTick', 0:1:3);
%         text(-8, 3.4, ixnames{i1}, 'FontSize', fdim.labelfs, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
%     elseif ixy == 2
%         set(gca, 'YLim', [0, 6.5], 'YTick', 0:2:6);
%         text(-8, 6.2, ixnames{i1}, 'FontSize', fdim.labelfs, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
%     elseif ixy == 3
%         set(gca, 'YLim', [0, 1.3], 'YTick', 0:0.4:1.2);
%         text(-8, 1.25, ixnames{i1}, 'FontSize', fdim.labelfs, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
%     end
%     
%     if ixx > 1
%         set(gca, 'YTickLabel', '');
%     end
%     if ixy < 3
%         set(gca, 'XTickLabel', '');
%     end
%     
%     if i1 == 4
%         ylabel('concentration, mM', 'FontSize', fdim.labelfs);
%     elseif i1 == 8
%         xlabel('time, s', 'FontSize', fdim.labelfs);
%     end
%     
% end
% 
% % %%% PANEL (b)
% % subplot('Position', [fdim.spxvec(1) fdim.spyvec(2) fdim.spwr fdim.sphr]),
% % hold on, box on;
% % set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
% % % NAD
% % plot( t, x(:,ix(4)), 'k-', 'LineWidth', 2, 'Marker', 'none');
% % % NADH
% % plot( t, x(:,ix(5)), 'k:', 'LineWidth', 2, 'Marker', 'none');
% % set(gca, 'YLim', [0 2], 'YTick', 0:0.5:1.5, 'XLim', [-10 40], 'YGrid', 'on', 'XGrid', 'on', 'XTickLabel', {});
% % legend({'NAD', 'NADH'}, 'Location', 'North');
% % text(39, 1.95, '(b)', 'FontSize', fdim.labelfs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
% % ylabel('concentration, mM', 'FontSize', fdim.labelfs);
% % 
% % %%% PANEL (c)
% % subplot('Position', [fdim.spxvec(1) fdim.spyvec(3) fdim.spwr fdim.sphr]),
% % hold on, box on;
% % set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
% % % NADP
% % plot( t, x(:,ix(6)), 'k-', 'LineWidth', 2, 'Marker', 'none');
% % % NADPH
% % plot( t, x(:,ix(7)), 'k:', 'LineWidth', 2, 'Marker', 'none');
% % set(gca, 'YLim', [0.04 0.24], 'YTick', 0.04:0.04:0.20, 'XLim', [-10 40], 'YGrid', 'on', 'XGrid', 'on');
% % legend({'NADP', 'NADPH'}, 'Location', 'North');
% % text(39, 0.235, '(c)', 'FontSize', fdim.labelfs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
% % xlabel('time, s', 'FontSize', fdim.labelfs);
% 
% clear fdim cc;





%% MODIFYING THE MODEL FOR STEADY STATE CALCULATION

set(OrigModelObj.Species(1), 'ConstantAmount', true);
set(OrigModelObj.Species(1), 'InitialAmount', 0.0556);

for irule = 1:length(OrigModelObj.Rules)
    delete(OrigModelObj.Rules(1));
end

[success, variant_out, mWT] = sbiosteadystate(OrigModelObj);

mWT.Species
rxnWT = get_fluxes(mWT);




