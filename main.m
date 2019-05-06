%%% Investigating epistasis in the lower glycolysis pathway

clear;

currdir = '/Users/skryazhi/epistasis_theory/data/Chassagnole_etal';
runId = '2019-05-02';

Pert = 0.9;
Glu = 'Low'; % 'Low' 'Med' or 'Inf';
ModelSpecs.Type = 'GPPP';
% LG = lower glycolysis
% PPPSMALL = G6PDH and PGDH:    
% PPP = pentose phosphate pathway
% UGPPP = upper glycolysis and pentose phosphate pathway
% GPPP = upper and lower glycolysis and pentose phosphate pathway
% FULL = full model (w/o g1p and extreneous rxns)

AbsTolMut = 1e-4;
AbsTolEps = 1e-6;

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
    'PTS', 'Phosphotransferase system', ...
    'PK', 'Pyruvate kinase');

% Cext = Inf:
% 1         cytosol         Phosphoenol pyruvate          2.15868
% 2         cytosol         Glucose-6-Phosphate           3.7979
% 3         cytosol         Pyruvate                      2.67
% 4         cytosol         Fructose-6-Phosphate          0.653706
% 5         cytosol         6-Phosphogluconate            0.867364
% 6         cytosol         Fructose-1,6-bisphosphate     1.09056
% 7         cytosol         sedoheptulose-7-phosphate     0.244278
% 8         cytosol         Glyceraldehyde-3-Phosphate    0.423085
% 9         cytosol         Erythrose-4-phosphate         0.157633
% 10        cytosol         Xylulose-5-phosphate          0.178106
% 11        cytosol         Ribose-5-phosphate            0.517141
% 12        cytosol         Dihydroxyacetonephosphate     0.341033
% 13        cytosol         1,3-diphosphosphoglycerate    0.00667297
% 14        cytosol         3-Phosphoglycerate            1.74827
% 15        cytosol         2-Phosphoglycerate            0.325327
% 16        cytosol         Ribulose-5-phosphate          0.140996

% Cext = 10:
% 1         extracellular    Extracellular Glucose         10
% 2         cytosol          Phosphoenol pyruvate          3.51196
% 3         cytosol          Glucose-6-Phosphate           4.04183
% 4         cytosol          Pyruvate                      2.67
% 5         cytosol          Fructose-6-Phosphate          0.696757
% 6         cytosol          6-Phosphogluconate            0.911921
% 7         cytosol          Fructose-1,6-bisphosphate     0.722251
% 8         cytosol          sedoheptulose-7-phosphate     0.282199
% 9         cytosol          Glyceraldehyde-3-Phosphate    0.359287
% 10        cytosol          Erythrose-4-phosphate         0.144751
% 11        cytosol          Xylulose-5-phosphate          0.176212
% 12        cytosol          Ribose-5-phosphate            0.514247
% 13        cytosol          Dihydroxyacetonephosphate     0.275267
% 14        cytosol          1,3-diphosphosphoglycerate    0.010578
% 15        cytosol          3-Phosphoglycerate            2.81128
% 16        cytosol          2-Phosphoglycerate            0.525739
% 17        cytosol          Ribulose-5-phosphate          0.140553

% Cext = 5:
%    1         extracellular    Extracellular Glucose         5.02              
%    2         cytosol          Phosphoenol pyruvate          3.88449           
%    3         cytosol          Glucose-6-Phosphate           3.96602           
%    4         cytosol          Pyruvate                      2.67              
%    5         cytosol          Fructose-6-Phosphate          0.683989          
%    6         cytosol          6-Phosphogluconate            0.898188          
%    7         cytosol          Fructose-1,6-bisphosphate     0.464476          
%    8         cytosol          sedoheptulose-7-phosphate     0.297704          
%    9         cytosol          Glyceraldehyde-3-Phosphate    0.290969          
%    10        cytosol          Erythrose-4-phosphate         0.125116          
%    11        cytosol          Xylulose-5-phosphate          0.163073          
%    12        cytosol          Ribose-5-phosphate            0.478687          
%    13        cytosol          Dihydroxyacetonephosphate     0.219815          
%    14        cytosol          1,3-diphosphosphoglycerate    0.0116159         
%    15        cytosol          3-Phosphoglycerate            3.09897           
%    16        cytosol          2-Phosphoglycerate            0.58036           
%    17        cytosol          Ribulose-5-phosphate          0.131176          



%%% Lower glycolysis:
if strcmp(ModelSpecs.Type, 'LG') 
    mutList = [ ...
        {{'GAPDH'}, Pert};
        {{'PGK'}, Pert};
        {{'PGluMu'},  Pert};
        {{'ENO'}, Pert};
        ];
    ModelSpecs.rxnOutEnzymes = {'ENO'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1];
    
    if strcmp(Glu, 'Inf')
        % For Cext = Inf
        ModelSpecs.ExtMetConc.gap = 0.423085; % 0.218;
        ModelSpecs.ExtMetConc.pep = 2.15868;  % 2.67;
    elseif strcmp(Glu, 'Med')
        % For Cext = 10
        ModelSpecs.ExtMetConc.gap = 0.359287; % 0.432208; % 0.218;
        ModelSpecs.ExtMetConc.pep = 3.51196; % 2.13742;  % 2.67;
    elseif strcmp(Glu, 'Low')
        % For Cext = 5.02
        ModelSpecs.ExtMetConc.gap = 0.290969; % 0.218;
        ModelSpecs.ExtMetConc.pep = 3.88449;  % 2.67;
    end
    
    ModelSpecs.IFREMOVE6PGINH = false;

%%% G6PDH and PGDH:    
elseif strcmp(ModelSpecs.Type, 'PPPSMALL') 
    mutList = [ ...        
        {{'G6PDH'}, Pert};
        {{'PGDH'},  Pert};
        ];
    ModelSpecs.rxnOutEnzymes = {'PGDH'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1]';

    if strcmp(Glu, 'Inf')
        % For Cext = Inf
        ModelSpecs.ExtMetConc.g6p = 3.7979 ; % 3.48        
        ModelSpecs.ExtMetConc.ribu5p = 0.140996; % 0.111
        % ModelSpecs.ExtMetConc.f6p = 0.653706; % 0.6;
        % ModelSpecs.ExtMetConc.gap = 0.423085; % 0.218;
    elseif strcmp(Glu, 'Med')
        % For Cext = 10
        ModelSpecs.ExtMetConc.g6p = 4.04183 ; % 3.79382 ; % 3.48        
        ModelSpecs.ExtMetConc.ribu5p = 0.140553; % 0.111
        % ModelSpecs.ExtMetConc.f6p = 0.696757; % 0.6;
        % ModelSpecs.ExtMetConc.gap = 0.359287; % 0.218;
    elseif strcmp(Glu, 'Low')
        % For Cext = 5.02
        ModelSpecs.ExtMetConc.g6p = 3.96602;  % 3.48        
        ModelSpecs.ExtMetConc.ribu5p = 0.131176; % 0.111
        % ModelSpecs.ExtMetConc.f6p = 0.683989; % 0.6;
        % ModelSpecs.ExtMetConc.gap = 0.290969; % 0.218;
    end
    
    ModelSpecs.IFREMOVE6PGINH = false;    

%%% PPP:    
elseif strcmp(ModelSpecs.Type, 'PPP') 
    mutList = [ ...        
        {{'G6PDH'}, Pert};
        {{'PGDH'},  Pert};
%         {{'R5PI'}, Pert};
%         {{'Ru5P'}, Pert};
%         {{'TKa'},  Pert};
%         {{'TA'},  Pert};
%         {{'TKb'}, Pert};
        ];
    ModelSpecs.rxnOutEnzymes = {'TKb', 'TKa', 'TA'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1, +1, -1]';

    if strcmp(Glu, 'Inf')
        % For Cext = Inf
        ModelSpecs.ExtMetConc.g6p = 3.7979 ; % 3.48        
        % ModelSpecs.ExtMetConc.ribu5p = 0.140996; % 0.111
        ModelSpecs.ExtMetConc.f6p = 0.653706; % 0.6;
        ModelSpecs.ExtMetConc.gap = 0.423085; % 0.218;
    elseif strcmp(Glu, 'Med')
        % For Cext = 10
        ModelSpecs.ExtMetConc.g6p = 4.04183 ; % 3.79382 ; % 3.48        
        % ModelSpecs.ExtMetConc.ribu5p = 0.140553; % 0.111
        ModelSpecs.ExtMetConc.f6p = 0.696757; % 0.6;
        ModelSpecs.ExtMetConc.gap = 0.359287; % 0.218;
    elseif strcmp(Glu, 'Low')
        % For Cext = 5.02
        ModelSpecs.ExtMetConc.g6p = 3.96602;  % 3.48        
        % ModelSpecs.ExtMetConc.ribu5p = 0.131176; % 0.111
        ModelSpecs.ExtMetConc.f6p = 0.683989; % 0.6;
        ModelSpecs.ExtMetConc.gap = 0.290969; % 0.218;
    end
    
    ModelSpecs.IFREMOVE6PGINH = false;

%%% Upper glycolysis and PPP:    
elseif strcmp(ModelSpecs.Type, 'UGPPP')
    
    mutList = [ ...
        %%% Upper glycolysis
        {{'PGI'}, Pert};
        {{'PFK'}, Pert};
%         {{'ALDO'}, Pert};
%         {{'TIS'}, Pert};
%         %%% PPP:
        {{'G6PDH'}, Pert};
        {{'PGDH'},  Pert};
%         {{'R5PI'}, Pert};
%         {{'Ru5P'}, Pert};
%         {{'TKa'},  Pert};
%         {{'TA'},  Pert};
%         {{'TKb'}, Pert};
        ];
    ModelSpecs.rxnOutEnzymes = {'ALDO', 'TIS', 'TKb', 'TKa', 'TA'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1, +1, +1, +1, -1]';

    if strcmp(Glu, 'Inf')
        % For Cext = Inf
        ModelSpecs.ExtMetConc.g6p = 3.7979 ; % 3.48
        ModelSpecs.ExtMetConc.gap = 0.423085; % 0.218;
        ModelSpecs.ExtMetConc.pep = 2.15868;  % 2.67;
    elseif strcmp(Glu, 'Med')
        % For Cext = 10
        ModelSpecs.ExtMetConc.g6p = 4.04183 ; % 3.79382 ; % 3.48
        ModelSpecs.ExtMetConc.gap = 0.359287; % 0.432208; % 0.218;
        ModelSpecs.ExtMetConc.pep = 3.51196; % 2.13742;  % 2.67;
    elseif strcmp(Glu, 'Low')
        % For Cext = 5.02
        ModelSpecs.ExtMetConc.g6p = 3.96602;  % 3.48
        ModelSpecs.ExtMetConc.gap = 0.290969; % 0.218;
        ModelSpecs.ExtMetConc.pep = 3.88449;  % 2.67;
    end
        
    ModelSpecs.IFREMOVE6PGINH = false;

%%% Upper and lower glycolysis and PPP:    
elseif strcmp(ModelSpecs.Type, 'GPPP')
    
    mutList = [ ...
%        %%% Upper glycolysis
        {{'PGI'}, Pert};
        {{'PFK'}, Pert};
%         {{'ALDO'}, Pert};
%         {{'TIS'}, Pert};
%       %%% Lower glycolysis
%       {{'GAPDH'}, Pert};
%        {{'PGK'}, Pert};
%       {{'PGluMu'},  Pert};
%       {{'ENO'}, Pert};
%         %%% PPP:
        {{'G6PDH'}, Pert};
         {{'PGDH'},  Pert};
%         {{'R5PI'}, Pert};
%         {{'Ru5P'}, Pert};
%         {{'TKa'},  Pert};
%         {{'TA'},  Pert};
%         {{'TKb'}, Pert};
        ];
    ModelSpecs.rxnOutEnzymes = {'ENO'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1];

    if strcmp(Glu, 'Inf')
        % For Cext = Inf
        ModelSpecs.ExtMetConc.g6p = 3.7979;  % 3.48
        ModelSpecs.ExtMetConc.pep = 2.15868;  % 2.67;
    elseif strcmp(Glu, 'Med')
        % For Cext = 10
        ModelSpecs.ExtMetConc.g6p = 4.04183; % 3.48
        ModelSpecs.ExtMetConc.pep = 3.51196; % 2.67;
    elseif strcmp(Glu, 'Low')
        % For Cext = 5.02
        ModelSpecs.ExtMetConc.g6p = 3.96602; % 3.48
        ModelSpecs.ExtMetConc.pep = 3.88449; % 2.67;
    end          
    
    ModelSpecs.IFREMOVE6PGINH = false;
    
elseif strcmp(ModelSpecs.Type, 'FULL')
    
    mutList = [ ...
        %        %%% Upper glycolysis
        {{'PGI'}, Pert};
        {{'PFK'}, Pert};
%         {{'ALDO'}, Pert};
%         {{'TIS'}, Pert};
%        %%% Lower glycolysis
%        {{'GAPDH'}, Pert};
%         {{'PGK'}, Pert};
%         {{'PGluMu'},  Pert};
%         {{'ENO'}, Pert};
%        %%% PPP:
        {{'G6PDH'}, Pert};
        {{'PGDH'},  Pert};
%        {{'R5PI'}, Pert};
%        {{'Ru5P'}, Pert};
%        {{'TKa'},  Pert};
%        {{'TA'},  Pert};
%        {{'TKb'}, Pert};
        ];
    ModelSpecs.rxnOutEnzymes = {'PK', 'PTS'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1 +1]';
    
    if strcmp(Glu, 'Inf')
        % For Cext = Inf
        ModelSpecs.ExtMetConc.glu = Inf; %0.0556 ;
        ModelSpecs.ExtMetConc.pyr = 2.67; % 2.67;
    elseif strcmp(Glu, 'Med')
        % For Cext = 10
        ModelSpecs.ExtMetConc.glu = 10 ; %0.0556 ;
        ModelSpecs.ExtMetConc.pyr = 2.67; % 2.67;
    elseif strcmp(Glu, 'Low')
        % For Cext = 5.02
        ModelSpecs.ExtMetConc.glu = 5.02; % 5.02; %0.0556 ;
        ModelSpecs.ExtMetConc.pyr = 2.67; % 2.67;
    end
    
    % ModelSpecs.ExtMetConc.g1p = 0.653;
    
    ModelSpecs.IFREMOVE6PGINH = false;
else
    error('Unknown model');
end


%% Load model:
filename = sprintf('%s/BIOMD0000000051.xml', currdir);

OrigModelObj = sbmlimport(filename);

%% Modify model

%%% 1. Make the co-metabolite (ATP, ADP, etc) concentrations constant:
for irule = 1:length(OrigModelObj.Rules)
    delete(OrigModelObj.Rules(1));
end
clear irule;

%%% 2. Remove dilution by growth:
for irxn = 48:-1:31
    delete(OrigModelObj.Reaction(irxn));
end
clear irxn;

if strcmp(ModelSpecs.Type, 'PGMENO')
    
    %%% 3. Remove extreneous rxns:
    for irxn = sort([1 3 30 29 9 14 13 16 25 24 19:23 2 5 10 12 4 26:28 6:8 11 15], 'descend')
        delete(OrigModelObj.Reaction(irxn));
    end
    clear irxn;
        
    %%% 4. Remove extreneous metabolites
    for isp = sort([1 4 6, 3 5 7 8 9 11 12 13 14 18 10 15], 'descend')
        delete(OrigModelObj.Species(isp));
    end
    
    %%% 5. Delete extracellular compartment
    set(OrigModelObj.Compartments(2), 'Owner', []);
    delete(OrigModelObj.Compartments(1));
    
elseif strcmp(ModelSpecs.Type, 'LG')
    
    %%% 3. Remove extreneous rxns:
    for irxn = sort([1 3 30 29 9 14 13 16 25 24 19:23 2 5 10 12 4 26:28 6:8], 'descend')
        delete(OrigModelObj.Reaction(irxn));
    end
    clear irxn;
        
    %%% 4. Remove extreneous metabolites
    for isp = sort([1 4 6, 3 5 7 8 9 11 12 13 14 18], 'descend')
        delete(OrigModelObj.Species(isp));
    end
    
    %%% 5. Delete extracellular compartment
    set(OrigModelObj.Compartments(2), 'Owner', []);
    delete(OrigModelObj.Compartments(1));
        
    
elseif strcmp(ModelSpecs.Type, 'PPPSMALL')    
    %%% 3. Remove extreneous rxns:
    for irxn = sort( setdiff(1:30, [4 26:28]) , 'descend')
        delete(OrigModelObj.Reaction(irxn));
    end
    clear irxn;
        
    %%% 4. Remove extreneous metabolites
    for isp = sort(setdiff(1:18, [3, 7, 18 12 13]), 'descend')
        delete(OrigModelObj.Species(isp));
    end
    
    %%% 5. Delete extracellular compartment
    set(OrigModelObj.Compartments(2), 'Owner', []);
    delete(OrigModelObj.Compartments(1));
    
elseif strcmp(ModelSpecs.Type, 'PPP')    
    %%% 3. Remove extreneous rxns:
    for irxn = sort([1 3 30 29 9 14 13 16 25 24 19:23, 2 5 10 12 11 15 17 18], 'descend')
        delete(OrigModelObj.Reaction(irxn));
    end
    clear irxn;
        
    %%% 4. Remove extreneous metabolites
    for isp = sort([1 4 6, 2 8 14:17], 'descend')
        delete(OrigModelObj.Species(isp));
    end
    
    %%% 5. Delete extracellular compartment
    set(OrigModelObj.Compartments(2), 'Owner', []);
    delete(OrigModelObj.Compartments(1));

elseif strcmp(ModelSpecs.Type, 'UGPPP')
    
    %%% 3. Releasing inhibition of PGI by 6pg
    if ModelSpecs.IFREMOVE6PGINH
        set(OrigModelObj.Reactions(2), 'ReactionRate',...
            'cytosol*rmaxPGI*(cytosol.[Glucose-6-Phosphate]-cytosol.[Fructose-6-Phosphate]/KPGIeq)/(KPGIg6p*(1+cytosol.[Fructose-6-Phosphate]/(KPGIf6p*(1+cytosol.[6-Phosphogluconate]/KPGIf6ppginh)) )+cytosol.[Glucose-6-Phosphate])');
        delete(OrigModelObj.Reactions(2).KineticLaw.Parameters(6));
    end

    %%% 4. Remove extreneous rxns:
    for irxn = sort([1 3 30 29 9 14 13 16 25 24 19:23 11 15 17 18], 'descend')
        delete(OrigModelObj.Reaction(irxn));
    end
    clear irxn;
        
    %%% 5. Remove extreneous metabolites
    for isp = sort([1 4 6 15:17], 'descend')
        delete(OrigModelObj.Species(isp));
    end

    %%% 5. Delete extracellular compartment
    set(OrigModelObj.Compartments(2), 'Owner', []);
    delete(OrigModelObj.Compartments(1));
    
elseif strcmp(ModelSpecs.Type, 'GPPP')
    
    %%% 3. Releasing inhibition of PGI by 6pg
    if ModelSpecs.IFREMOVE6PGINH
        set(OrigModelObj.Reactions(2), 'ReactionRate',...
            'cytosol*rmaxPGI*(cytosol.[Glucose-6-Phosphate]-cytosol.[Fructose-6-Phosphate]/KPGIeq)/(KPGIg6p*(1+cytosol.[Fructose-6-Phosphate]/(KPGIf6p*(1+cytosol.[6-Phosphogluconate]/KPGIf6ppginh)) )+cytosol.[Glucose-6-Phosphate])');
        delete(OrigModelObj.Reactions(2).KineticLaw.Parameters(6));
    end
    
    %%% 4. Remove extreneous rxns:
    for irxn = sort([1 3 30 29 9 14 13 16 25 24 19:23], 'descend')
        delete(OrigModelObj.Reaction(irxn));
    end
    clear irxn;
    
    %%% 5. Remove extreneous metabolites
    for isp = sort([1 4 6], 'descend')
        delete(OrigModelObj.Species(isp));
    end

    %%% 6. Delete extracellular compartment
    set(OrigModelObj.Compartments(2), 'Owner', []);
    delete(OrigModelObj.Compartments(1));

elseif strcmp(ModelSpecs.Type, 'FULL')
    
    %%% 3. Releasing inhibition of PGI by 6pg
    if ModelSpecs.IFREMOVE6PGINH
        set(OrigModelObj.Reactions(2), 'ReactionRate',...
            'cytosol*rmaxPGI*(cytosol.[Glucose-6-Phosphate]-cytosol.[Fructose-6-Phosphate]/KPGIeq)/(KPGIg6p*(1+cytosol.[Fructose-6-Phosphate]/(KPGIf6p*(1+cytosol.[6-Phosphogluconate]/KPGIf6ppginh)) )+cytosol.[Glucose-6-Phosphate])');
        delete(OrigModelObj.Reactions(2).KineticLaw.Parameters(6));
    end
    
    %%% 4. Remove extreneous rxns:
    for irxn = sort([3 30 9 13 14 16 21 22 24 25 29], 'descend')
        delete(OrigModelObj.Reaction(irxn));
    end
    clear irxn;
    
       
    %% 5. Remove extreneous metabolites
    for isp = sort([6], 'descend')
        delete(OrigModelObj.Species(isp));
    end

    % 6. Changing the weird stoichiometry of the PTS rxn:
    set(OrigModelObj.Reactions(1), 'Stoichiometry', [-1 -1 1 1]);
    
    % ModelSpecs.ExtMetConc.glu == Inf
end
clear irxn isp;

%%% 7. Set initial amounts
s = [fieldnames(MetNames), struct2cell( MetNames )];
for imet = length(OrigModelObj.Species):-1:1
    ix = find( strcmp( s(:,2), get(OrigModelObj.Species(imet), 'Name') ) );
    if isempty(ix)
        continue;
    end
    
    if isfield(ModelSpecs.ExtMetConc, s{ix,1})
        if strcmp(s{ix,1}, 'glu') && (ModelSpecs.ExtMetConc.(s{ix,1}) == Inf)
            set(OrigModelObj.Reactions(1), 'Reaction', 'cytosol.[Phosphoenol pyruvate] <-> cytosol.[Glucose-6-Phosphate] + cytosol.Pyruvate');
            set(OrigModelObj.Reactions(1), 'Stoichiometry', [-1 1 1]);
            set(OrigModelObj.Reactions(1), 'ReactionRate',...
                'rmaxPTS*(cytosol.[Phosphoenol pyruvate]/cytosol.Pyruvate)/((KPTSa3+(cytosol.[Phosphoenol pyruvate]/cytosol.Pyruvate))*(1+power(cytosol.[Glucose-6-Phosphate],nPTSg6p)/KPTSg6p))');
            delete(OrigModelObj.Species(1));
            
            %%% 6. Delete extracellular compartment
            set(OrigModelObj.Compartments(2), 'Owner', []);
            delete(OrigModelObj.Compartments(1));
        else
            set(OrigModelObj.Species(imet), 'InitialAmount', ModelSpecs.ExtMetConc.(s{ix,1}) );
            set(OrigModelObj.Species(imet), 'ConstantAmount', true);
        end
    end
end
clear s imet ix;

%     %%% 6. Set all coefficients to 1
%     for irxn = 1:length(OrigModelObj.Reactions)
%         for iparam = 1:length(OrigModelObj.Reactions(irxn).KineticLaw.Parameters)
%             set(OrigModelObj.Reactions(irxn).KineticLaw.Parameters(iparam), 'Value', 1);
%         end
%     end


% We will monitor this flux as the output
RxnNames = arrayfun(@(x) get(x, 'Name'), OrigModelObj.Reactions, 'UniformOutput', false);
rxnOutIdVec = cellfun(@(EnzShortName) find(strcmp(RxnNames,EnzNames.(EnzShortName))) ,ModelSpecs.rxnOutEnzymes);


%% Find wildtype steady state

% [success, variant_out, ModelL] = sbiosteadystate(OrigModelObj, 'MaxStopTime', 1e7);
% 
% if ~success
%     fprintf('Failed to converge to steady state!\n');
% end

[WT.m, WT.FluxDistr] = getMutFlux(OrigModelObj, []);
WT.OutFlux = sum( ModelSpecs.rxnOutEnzymeCoeff .* WT.FluxDistr.Flux(rxnOutIdVec) );

WT.m.Species
WT.FluxDistr.Flux

%% Single mutants calculation:

FullEnzNames = get_full_names(EnzNames, mutList(:,1));
[ModelL, FluxDistrL] = getMutFlux(WT.m, [FullEnzNames, mutList(:,2)]);
mutList = [mutList, ModelL, FluxDistrL];
% Col 1 = enzyme short name list
% Col 2 = perturbation list
% Col 3 = perturbed model
% Col 4 = perturbed steady state fluxes
% Col 5 = perturbed steady state output flux
% Col 6 = relative perturbation of output flux

nMut = size( mutList, 1);


%% Single mutants output:
fprintf('WT:\nfOUT = %.3f \n', WT.OutFlux);

for iMut = 1:nMut
    mutList{iMut, 5} = sum( ModelSpecs.rxnOutEnzymeCoeff .* mutList{iMut,4}.Flux(rxnOutIdVec) );
    mutList{iMut, 6} = (mutList{iMut, 5} - WT.OutFlux)/WT.OutFlux;
        
    fprintf('\\delta %s = %.3g\n\f=> fOUT = %.3f\n\t=> \\delta f = %.2e\n',...
        mutList{iMut, 1}{1}, mutList{iMut, 2}, ...
        mutList{iMut, 5},...
        mutList{iMut, 6} );
end


%% Double mutants

n2Mut = nMut * (nMut - 1)/2;
mut2List = cell(n2Mut, 8);

i2Mut = 1;
for iMut1 = 1:nMut
    for iMut2 = (iMut1+1):nMut
        mut2List{i2Mut,1} = [mutList{iMut1,1}; mutList{iMut2,1}];
        mut2List{i2Mut,2} = [mutList{iMut1,2}; mutList{iMut2,2}];
        mut2List{i2Mut,7} = [iMut1; iMut2];
        i2Mut = i2Mut + 1;
    end
end

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


%% Write epistasis matrix
epsMat = nan(nMut, nMut);

for iMut = 1:n2Mut
    mut2List{iMut, 5} = sum( ModelSpecs.rxnOutEnzymeCoeff .* mut2List{iMut,4}.Flux(rxnOutIdVec) );
    mut2List{iMut, 6} = (mut2List{iMut, 5} - WT.OutFlux)/WT.OutFlux;
        
    fprintf('\\delta %s = %.3g, \\delta %s = %.3g\n\t=> fOUT = %.3f\n\t=> \\delta f = %.4g\n',...
        mut2List{iMut, 1}{1}, mut2List{iMut, 2}(1), ...
        mut2List{iMut, 1}{2}, mut2List{iMut, 2}(2), ...
        mut2List{iMut, 5},...
        mut2List{iMut, 6} );
    
    iMut1 = mut2List{iMut,7}(1);
    iMut2 = mut2List{iMut,7}(2);
    
    delta1 = mutList{iMut1,6};
    delta2 = mutList{iMut2,6};
    
    eps = mut2List{iMut, 6} - delta1 - delta2;
	if abs(delta1) < AbsTolEps || abs(delta2) < AbsTolEps || abs(mut2List{iMut, 6}) < AbsTolEps || ...
        abs(-1-delta1) < AbsTolEps || abs(-1-delta2) < AbsTolEps || abs(-1-mut2List{iMut, 6}) < AbsTolEps
        eps = nan;
    elseif abs(eps) < AbsTolEps
        eps = 0;
    else
        eps = eps/delta1/delta2;
    end
    epsMat(iMut1, iMut2) = eps;
    mut2List{iMut,8} = eps;
    fprintf('\t eps f = %.3g\n', eps);
end
clear i2Mut iMut1 iMut2 iPert rxnId delta1 delta2 eps





%% Saving
filename =  sprintf('%s/data %s/%s/data %s', currdir, runId, Glu, ModelSpecs.Type);

for met = fieldnames(ModelSpecs.ExtMetConc)'
    filename = sprintf('%s, [%s] = %.2f', filename, met{1}, ModelSpecs.ExtMetConc.(met{1}));
end

filename = sprintf('%s, Pert = %.1f', filename, Pert);

if ModelSpecs.IFREMOVE6PGINH
    filename = sprintf('%s, no6pginh.mat', filename);
else
    filename = sprintf('%s.mat', filename);
end

save(filename, 'ModelSpecs', 'OrigModelObj', 'WT', 'mutList', 'mut2List', 'epsMat');





%% Visualize epistasis matrix
clf;

%%% Specify the following dimensions:
fdim.spwa = 8; % subplotwidth in cm
fdim.spha = 8; % subplotheight in cm

fdim.nx = 1; % number of panels along the horizontal dimension
fdim.ny = 1; % number of panels along the vertical dimension

fdim.xma = [1.5 0.5]; % left right horizontal margin in cm
fdim.yma = [1.3 0.5]; % bottom top vertical margin cm

fdim.dxa = 0.3; % horizontal distance between panels in cm
fdim.dya = 0.4; % vertical distance between panels in cm

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

subplot('Position', [fdim.spxvec(1) fdim.spyvec(1) fdim.spwr fdim.sphr]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');

fdim.ms = 1;

nMut = size(epsMat,1);
hpos = [];
hneg = [];

for iMut1 = 1:nMut
    xpos = iMut1;
    ypos = nMut-(iMut1-1);
    
    delta1 = mutList{iMut1,6};

    text(xpos, ypos, sprintf('%.1e', delta1),...
        'HorizontalAlignment','center', 'VerticalAlignment', 'middle',...
        'FontName', 'Helvetica', 'FontSize', fdim.labelfs); 
end


for iMut1 = 1:nMut-1
    for iMut2 = iMut1+1:nMut
                
        if isnan(epsMat(iMut1,iMut2))
            c = 0.4 * [1 1 1];
        elseif epsMat(iMut1,iMut2) < 0
            c = cc(2,:);
            % hneg(end+1) = plot( iMut2, nMut-(iMut1-1), 's', 'MarkerSize', 22, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', c);
        elseif epsMat(iMut1,iMut2) > 0
            c = cc(6,:);
            % hpos(end+1) = plot( iMut2, nMut-(iMut1-1), 's', 'MarkerSize', 22, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', c);
        else 
            c = cc(7,:);
        end
        xpos = iMut2;
        ypos = nMut-(iMut1-1);
        
        rectangle('Position', [xpos-fdim.ms/2, ypos-fdim.ms/2, fdim.ms, fdim.ms],...
            'EdgeColor', 'none', 'FaceColor', c);
        text(xpos, ypos, sprintf('%.1f', epsMat(iMut1,iMut2)),...
            'HorizontalAlignment','center', 'VerticalAlignment', 'middle',...
            'FontName', 'Helvetica', 'FontSize', fdim.labelfs);
        
        delta1 = mutList{iMut1,6};
        delta2 = mutList{iMut2,6};

        xpos = iMut1;
        ypos = nMut-(iMut2-1);
                
        rectangle('Position', [xpos-fdim.ms/2, ypos-fdim.ms/2, fdim.ms, fdim.ms],...
            'EdgeColor', 'none', 'FaceColor', c);
        text(xpos, ypos, sprintf('%.1e', epsMat(iMut1,iMut2)*delta1*delta2),...
            'HorizontalAlignment','center', 'VerticalAlignment', 'middle',...
            'FontName', 'Helvetica', 'FontSize', fdim.labelfs);
    end
end

for pos = 1.5:1:nMut-0.5
    plot( [0.5 nMut+0.5], pos*[1 1], '-', 'Color', 0.8*[1 1 1]);
    plot( pos*[1 1], [0.5 nMut+0.5], '-', 'Color', 0.8*[1 1 1]);
end

MutEnzList = cellfun(@(x) x{1}, mutList(:,1), 'UniformOutput', false);

set(gca, 'XLim', [0.5 nMut+0.5], 'YLim', [0.5 nMut+0.5], ... %'XGrid', 'on', 'YGrid', 'on', ...
    'XTick', 1:nMut, 'YTick', 1:nMut, 'XTickLabel', MutEnzList, 'YTickLabel', flipud(MutEnzList))
xtickangle(90);
% legend([hpos(1), hneg(1)], 'Pos', 'Neg', 'Location', 'SouthWest');

clear fdim cc hpos hneg iMut1 iMut2 c pos;

%% Saving
filename =  sprintf('%s/figures_%s/eps/%s/%s', currdir, runId, Glu, ModelSpecs.Type);

for met = fieldnames(ModelSpecs.ExtMetConc)'
    filename = sprintf('%s, [%s] = %.2f', filename, met{1}, ModelSpecs.ExtMetConc.(met{1}));
end

filename = sprintf('%s, Pert = %.1f', filename, Pert);

if ModelSpecs.IFREMOVE6PGINH
    filename = sprintf('%s, no6pginh.eps', filename);
else
    filename = sprintf('%s.eps', filename);
end

saveas(gcf, filename, 'epsc');


