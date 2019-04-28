%% Load and modify the model:
% currdir = '/Users/skryazhi/epistasis_theory/data/Chassagnole_etal';
% filename = sprintf('%s/BIOMD0000000051.xml', currdir);
% 
% OrigModelObj = sbmlimport(filename);
% 
% % Set a constant amount of extracellular glucose and make it infinite
% delete(OrigModelObj.Reaction(48));
% % delete(OrigModelObj.Reactions(1).Reactants(1));
% set(OrigModelObj.Reactions(1), 'Reaction', 'cytosol.[Phosphoenol pyruvate] <-> cytosol.[Glucose-6-Phosphate] + cytosol.Pyruvate');
% set(OrigModelObj.Reactions(1), 'Stoichiometry', [-1 1 1]);
% set(OrigModelObj.Reactions(1), 'ReactionRate', 'rmaxPTS*(cytosol.[Phosphoenol pyruvate]/cytosol.Pyruvate)/((KPTSa3+(cytosol.[Phosphoenol pyruvate]/cytosol.Pyruvate))*(1+power(cytosol.[Glucose-6-Phosphate],nPTSg6p)/KPTSg6p))');
% delete(OrigModelObj.Species(1));
% 
% % Make the co-metabolite (ATP, ADP, etc) concentrations constant:
% for irule = 1:length(OrigModelObj.Rules)
%     delete(OrigModelObj.Rules(1));
% end
% clear irule;
% 
% % Releasing inhibition of PGI by 6pg
% % set(OrigModelObj.Reactions(2).KineticLaw.Parameters(6), 'Value', 20);
% 
% % We will monitor this flux as the output
% rxnOutId = 24; % PDH
% 
% [success, variant_out, mWT] = sbiosteadystate(OrigModelObj);
% 
% WTrxn = get_fluxes(mWT, 'v');
% 
% WTflux = getMutFlux(OrigModelObj, rxnOutId);

 mutList.Id = [ ... 
    {1}; ... % PTS
      {2};... % PGI
     {5}; ... % PFK
     {10}; ...% ALDO
     {11}; ...% GAPDH
     {15}; ...% PGK
     {17}; ...% PGluMu
     {18}; ...% ENO
     {19}; ...% PK
%     ----------
     {4}; ... % G6PDH
     {26}; ... % PGDH
     {27}; ... % R5P1
     {28}; ... % Ru5p
     {7}; ...  % TKa
     {6}; ... % TA
     {8} ...  % TKb
    ];



% % Works at CextGlu = 1:
% mutList.Pert = [ ... 
%     {0.1}; ...% PTS
%     {0.1}; ...% PGI
%     {0.1};... % PFK
%     {0.1}; ...% ALDO
%     {0.1}; ...% GAPDH
%     {0.1}; ...% PGK
%     {0.1}; ...% PGluMu
%     {0.1}; ...% ENO
%     {10}; ...% PK
%     {10}; ...% G6PDH
%     {0.1}; ...% PGDH
%     {10}; ...% R5PI
%     {0.1}; ...% Ru5p
%     {0.1}; ...% TKa
%     {0.1}; ...% TA
%     {0.1}; ...% TKb
%     ];


% % Works at CextGlu = 2:
% mutList.Pert = [ ... 
%     {0.9}; ...% PTS
%     {0.01}; ...% PGI
%     {0.9};... % PFK
%     {0.025}; ...% ALDO
%     {0.25}; ...% GAPDH
%     {0.01}; ...% PGK
%     {0.01}; ...% PGluMu
%     {0.01}; ...% ENO
%     {0.25}; ...% PK
%     {0.1}; ...% G6PDH
%     {0.075}; ...% PGDH
%     {10}; ...% R5PI
%     {0.1}; ...% Ru5p
%     {0.01}; ...% TKa
%     {0.05}; ...% TA
%      {0.1}; ...% TKb
%     ];

% Works at CextGlu = 10:
mutList.Pert = [ ... 
   {0.9}; ...% PTS
    {0.2}; ...% PGI
    {0.9};... % PFK
    {10}; ...% ALDO
    {0.1}; ...% GAPDH
    {10}; ...% PGK
    {10}; ...% PGluMu
    {10}; ...% ENO
    {0.25}; ...% PK
    {0.1}; ...% G6PDH
    {0.075}; ...% PGDH
    {0.01}; ...% R5PI
    {10}; ...% Ru5p
    {0.001}; ...% TKa
    {0.001}; ...% TA
     {10}; ...% TKb
    ];

% Works at CextGlu = 1e6:
% mutList.Pert = [ ... 
%     {0.9}; ...% PTS
%     {0.1}; ...% PGI
%     {0.9};... % PFK
%     {0.05}; ...% ALDO
%     {0.25}; ...% GAPDH
%     {0.05}; ...% PGK
%     {0.01}; ...% PGluMu
%     {0.01}; ...% ENO
%     {0.25}; ...% PK
%     {0.1}; ...% G6PDH
%     {0.075}; ...% PGDH
%     {0.01}; ...% R5PI
%     {10}; ...% Ru5p
%     {0.001}; ...% TKa
%     {0.001}; ...% TA
%      {100}; ...% TKb
%     ];


close all;

plotFluxDistr(WTrxn);
fprintf('WT:\nfIN = %.3f \n', WTrxn.Flux(1));


% Single mut A:
rxnId1 = 15;
Pert1 = 0.001;

fprintf('Mut A:\n');
MutModelObj = copyobj( OrigModelObj );
x = get( MutModelObj.Reactions(rxnId1).KineticLaw.Parameters(1), 'Value') * Pert1;
set(MutModelObj.Reactions(rxnId1).KineticLaw.Parameters(1), 'Value', x  );
[success, variant_out, m2] = sbiosteadystate(MutModelObj);
MUTrxn = get_fluxes(m2);
MUTflux = MUTrxn.Flux(rxnOutId);
plotFluxDistr(MUTrxn);
fprintf('%s\t%.2g\n', get( OrigModelObj.Reactions(rxnId1).KineticLaw.Parameters(1), 'Name'), Pert1 );
deltaFlux1 = (MUTflux - WTflux) / WTflux;
fprintf('fIN = %.3f, \\delta f = %.4g\n', MUTrxn.Flux(1), deltaFlux1 );


% Single mut B:
rxnId2 = 11;
Pert2 = 0.1;

fprintf('Mut B:\n');
MutModelObj = copyobj( OrigModelObj );
x = get( MutModelObj.Reactions(rxnId2).KineticLaw.Parameters(1), 'Value') * Pert2;
set(MutModelObj.Reactions(rxnId2).KineticLaw.Parameters(1), 'Value', x  );
[success, variant_out, m2] = sbiosteadystate(MutModelObj);
MUTrxn = get_fluxes(m2);
MUTflux = MUTrxn.Flux(rxnOutId);
plotFluxDistr(MUTrxn);
fprintf('%s\t%.2g\n', get( OrigModelObj.Reactions(rxnId2).KineticLaw.Parameters(1), 'Name'), Pert2 );
deltaFlux2 = (MUTflux - WTflux) / WTflux;
fprintf('fIN = %.3f, \\delta f = %.4g\n', MUTrxn.Flux(1), deltaFlux2 );

% Double mut AB:
fprintf('Mut AB:\n');
MutModelObj = copyobj( OrigModelObj );
x = get( MutModelObj.Reactions(rxnId1).KineticLaw.Parameters(1), 'Value') * Pert1;
set(MutModelObj.Reactions(rxnId1).KineticLaw.Parameters(1), 'Value', x  );
x = get( MutModelObj.Reactions(rxnId2).KineticLaw.Parameters(1), 'Value') * Pert2;
set(MutModelObj.Reactions(rxnId2).KineticLaw.Parameters(1), 'Value', x  );
[success, variant_out, m2] = sbiosteadystate(MutModelObj);
MUTrxn = get_fluxes(m2);
MUTflux = MUTrxn.Flux(rxnOutId);
plotFluxDistr(MUTrxn);
fprintf('%s\t%.2g\n', get( OrigModelObj.Reactions(rxnId1).KineticLaw.Parameters(1), 'Name'), Pert1 );
fprintf('%s\t%.2g\n', get( OrigModelObj.Reactions(rxnId2).KineticLaw.Parameters(1), 'Name'), Pert2 );
deltaFlux = (MUTflux - WTflux) / WTflux;
fprintf('fIN = %.3f, \\delta f = %.4g\n', MUTrxn.Flux(1), deltaFlux );
epsFlux = (deltaFlux-deltaFlux1-deltaFlux2)/deltaFlux1/deltaFlux2;
fprintf('\\eps = %.4g\n', epsFlux );


%% Single mutants
% mutList.Flux = getMutFlux(OrigModelObj, rxnOutId, mutList);
% 
% nMut = length( mutList.Id );
% 
% mutList.deltaFlux = cell(nMut, 1);
% 
% for iMut = 1:nMut
%     rxnIdVec = mutList.Id{iMut};
%     nRxn = size( mutList.Pert{iMut}, 1 );
%     nPert = size( mutList.Pert{iMut}, 2 );
% 
%     mutList.deltaFlux{iMut} = nan(1, nPert);
% 
%     for iPert = 1:nPert
%         fprintf('===\n');
%         for iParam = 1:nRxn
%             rxnId = rxnIdVec(iParam);
%             fprintf('%s\t%.2g\t', get( OrigModelObj.Reactions(rxnId).KineticLaw.Parameters(1), 'Name'), mutList.Pert{iMut}(iParam, iPert) );
%         end
%         mutList.deltaFlux{iMut}(iPert) = (mutList.Flux{iMut}(iPert) - WTflux) / WTflux;
%         fprintf('\\delta f = %.4g\n', mutList.deltaFlux{iMut}(iPert));
%     end
% end
% 
% 
% %% Double mutants
% 
% n2Mut = nMut * (nMut - 1)/2;
% mut2List.Id = cell(n2Mut, 1);
% mut2List.MutIx = cell(n2Mut, 1);
% mut2List.Pert = cell(n2Mut, 1);
% 
% i2Mut = 1;
% for iMut1 = 1:nMut
%     for iMut2 = (iMut1+1):nMut
%         mut2List.MutIx{i2Mut} = [iMut1; iMut2];
%         mut2List.Id{i2Mut} = [mutList.Id{iMut1} ; mutList.Id{iMut2}];
%         mut2List.Pert{i2Mut} = [mutList.Pert{iMut1} ; mutList.Pert{iMut2}];
%         i2Mut = i2Mut + 1;
%     end
% end
% 
% mut2List.Flux = getMutFlux(OrigModelObj, rxnOutId, mut2List);
% 
% mut2List.deltaFlux = cell(nMut, 1);
% mut2List.epsFlux = cell(nMut, 1);
% 
% epsMat = nan(nMut, nMut);
% 
% for iMut = 1:n2Mut
%     rxnIdVec = mut2List.Id{iMut};
%     nParam = size( mut2List.Pert{iMut}, 1 );
%     nPert = size( mut2List.Pert{iMut}, 2 );
% 
%     mut2List.deltaFlux{iMut} = nan(1, nPert);
% 
%     for iPert = 1:nPert
%         fprintf('===\n');
%         for iParam = 1:nParam
%             rxnId = rxnIdVec(iParam);
%             fprintf('%s\t%.2g\n', get( OrigModelObj.Reactions(rxnId).KineticLaw.Parameters(1), 'Name'), mut2List.Pert{iMut}(iParam, iPert) );
%         end
%         mut2List.deltaFlux{iMut}(iPert) = (mut2List.Flux{iMut}(iPert) - WTflux) / WTflux;
%         fprintf('\\delta f = %.4g\n', mut2List.deltaFlux{iMut}(iPert));
%         
%         iMut1 = mut2List.MutIx{iMut}(1);
%         iMut2 = mut2List.MutIx{iMut}(2);
%         
%         delta1 = mutList.deltaFlux{iMut1}(iPert);
%         delta2 = mutList.deltaFlux{iMut2}(iPert);
%         
%         eps = (mut2List.deltaFlux{iMut}(iPert) - delta1 - delta2)/delta1/delta2;
%         epsMat(iMut1, iMut2) = eps;
%         mut2List.epsFlux{iMut}(iPert) = eps;
%         
%         fprintf('\\eps f = %.2g\n', eps);
%     end
% end
% clear i2Mut iMut1 iMut2 iPert rxnId delta1 delta2 eps rxnIdVec
% 
% filename = sprintf('%s/data Cextglu = Inf.mat', currdir);
% save(filename, 'OrigModelObj', 'mWT', 'WTflux', 'mutList', 'mut2List', 'epsMat');
% 
% 
% 
% %% 
% 
% RxnTopo = nan(16,16);
% % +1 = serial, -1 = parallel
% RxnTopo(1,2:end) = +1; % PTS vs everything
% RxnTopo(2,3:9) = +1;   % PGI vs glycolysis downstream
% RxnTopo(2,10:16) = -1; % PGI vs pentose phosphate
% RxnTopo(3,4:9) = +1;   % PFK vs glycolysis downstream
% RxnTopo(3,10:13) = +1; % PFK vs pentose phosphate up to xyl5p and rib5p
% RxnTopo(3,14:16) = -1; % PFK vs pentose phosphate (TKa, TKb, TA)
% RxnTopo(4,5:9) = +1;   % ALDO vs glycolysis downstream
% RxnTopo(4,10:13) = +1; % ALDO vs pentose phosphate up to xyl5p and rib5p
% RxnTopo(4,14:16) = -1; % ALDO vs pentose phosphate (TKa, TKb, TA)
% RxnTopo(5,6:16) = +1;   % GAPDH vs everything
% RxnTopo(6,7:16) = +1;   % PGK vs everything
% RxnTopo(7,8:16) = +1;   % PGluMu vs everything
% RxnTopo(8,9:16) = +1;   % ENO vs everything
% RxnTopo(9,10:16) = +1;  % PK vs everything
% RxnTopo(10,11:16) = +1;  % G6PDH vs pentose phosphate downstream
% RxnTopo(11,12:16) = +1;  % PGDH vs pentose phosphate downstream
% RxnTopo(12,13) = -1;     % R5PI vs Ru5P
% RxnTopo(12,14:16) = +1;   % R5PI vs pentose phosphate downstream
% RxnTopo(13,14:16) = +1;   % Ru5P vs pentose phosphate downstream
% RxnTopo(14,15) = +1;   % TKa vs TA
% RxnTopo(15,16) = +1;   % TA vs TKb
