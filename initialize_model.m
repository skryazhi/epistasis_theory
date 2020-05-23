function [WT, ModelSpecs] = initialize_model(modelfile, data_dir, ModelName, Glu)
%
% initialize_model initializes the model
%
% [WT, ModelSpecs] = initializde_model(data_dir, ModelName, Glu) loads the
% original model by Chassagnole et al, trims and modifies it according to
% the input variables to this function, calculates steady state and control
% coefficient for the output flux.
% 
% INPUT VARIABLES:
%
% modelfile is the pathway to the file containing the original model
%
% data_dir is a string that contains the directory where the model data
% are/will be stored
%
% ModelName is a strong that contains what model ot use. Must be one of the
% following.
%
% 'LG' is lower glycolysis: external metabolites are gap and pep; output
% flux is the flux towards pep,  i.e., flux through ENO;
% 'PPSMALL' is the linear pathway containing G6PDH and PGDH: external
% metabolites are g6p and ribu5p; output flux is the flux towards ribu5p,
% i.e., flux through PGDH;
% 'PP' is the pentose phosphage pathway: external metabolites are g6p, f6p
% and gap; output flux is flux towards gap, i.e., + TKb flux + TKa flux -
% TA flux;
% 'UGPP' is upper glycolysis plus pentose phosphate pathway: external
% metabolites are g6p, gap and pep; output flux is the flux towards gap,
% i.e, + ALDO flux + TIS flux + TKb flux + TKa flux - TA flux;
% 'GPP' is upper and lower glycolysis and pentose phosphate pathway:
% external metabolites are g6p and pep; output flux is flux towards pep,
% i.e., flux through ENO;
% 'FULL' is the full model: external metabolites are glu and pyr; output
% flux is the flux towards pyr, i.e. PK flux + PTS flux;
%    
%
% Glu is a string that contains the concentration of the external glucose.
% Either 'Low', 'Med' or 'Inf'. For the Low model, [glu] = 2 然, [pyr] = 10
% 然; for the Med model, [glu] = 20 然, [pyr] = 100 然; for the Inf model,
% [glu] = infinity, [pyr] = 1 mM
% 
% OUTPUT VARIABELS:
%
% WT is a struct with the following fields. 'm' is the initialized
% simbiology model where the initial concentrations of metabolites are at
% steady state. 'OutFlux' is the steady-state output flux. 'FluxDistr' is a
% struct with fields 'Flux' (values of steady-state flux through each rxn);
% 'Name' (list of full rxn names); 'FCC' (vector flux control coeficients);
% and 'FH' (normalized Hessian matrix)
%
% ModelSpecs is a struct with the following fields. 'Name' is the type of
% the model, identical to ModelName. 'MetNames' is a cell array with two
% columns: first column has short metabolite names; second column has full
% metabolite names. 'EnzNames' is a cell array with two columns: first
% column has short enzyme names; second column has full enzyme names.
% 'ExtMetConc' is a struct that contains the concenctrations of external
% metabolites. 'rxnOutEnzymes' is the cell array of enzymes definining the
% output flux. 'rxnOutEnzymeCoeff' is the list of coefficients with which
% the output fluxes must summed to get output flux. 'rxnOutIX' is the
% vector of indices of output reactions in the WT.m.Reactions array.


ModelSpecs.Name = ModelName;

ModelSpecs.MetNames = {...
    'g6p', 'Glucose-6-Phosphate'; ...
    'f6p', 'Fructose-6-Phosphate'; ...
    'gap', 'Glyceraldehyde-3-Phosphate'; ...
    'pep','Phosphoenol pyruvate'; ...
    'ribu5p', 'Ribulose-5-phosphate'; ...
    'glu', 'Extracellular Glucose'; ...
    'pyr', 'Pyruvate'; ...
    'g1p', 'Glucose-1-Phosphate'; ...
    'threepg', '3-Phosphoglycerate'; ...
    'pgp', '1,3-diphosphosphoglycerate'};

ModelSpecs.EnzNames = {... %%%% Upper glycolysis:
    'PGI', 'Glucose-6-phosphate isomerase'; ... 
    'PFK', 'Phosphofructokinase'; ...
    'ALDO', 'Aldolase'; ... %%% Lower glycolysis:
    'TIS', 'Triosephosphate isomerase'; ...
    'GAPDH', 'Glyceraldehyde-3-phosphate dehydrogenase'; ...
    'PGK', 'Phosphoglycerate kinase'; ...
    'PGluMu', 'Phosphoglycerate mutase'; ...
    'ENO', 'Enolase'; ... %%% PP:
    'G6PDH', 'Glucose-6-phosphate dehydrogenase'; ...
    'PGDH', '6-Phosphogluconate dehydrogenase'; ...
    'R5PI', 'Ribose-phosphate isomerase'; ...
    'Ru5P', 'Ribulose-phosphate epimerase'; ...
    'TKa', 'Transketolase a'; ...
    'TA', 'Transaldolase'; ...
    'TKb', 'Transketolase b'; ... %%% Other rxns:
    'PTS', 'Phosphotransferase system'; ...
    'PK', 'Pyruvate kinase'; ...
    'PEPCxylase', 'PEP carboxylase'};


%% Designating external metabolites:
%%% Lower glycolysis:
if strcmp(ModelSpecs.Name, 'LG') 
    ExtMetShortNameList = {'gap', 'pep'};
        
%%% G6PDH and PGDH:    
elseif strcmp(ModelSpecs.Name, 'PPSMALL') 
    ExtMetShortNameList = {'g6p', 'ribu5p'};
    
%%% PP:    
elseif strcmp(ModelSpecs.Name, 'PP') 
    ExtMetShortNameList = {'g6p', 'f6p', 'gap'};

%%% Upper glycolysis and PP:    
elseif strcmp(ModelSpecs.Name, 'UGPP')
    ExtMetShortNameList = {'g6p', 'gap', 'pep'};
    
%%% Upper and lower glycolysis and PP:    
elseif strcmp(ModelSpecs.Name, 'GPP')
        ExtMetShortNameList = {'g6p', 'pep'};
    
%%% Full model    
elseif strcmp(ModelSpecs.Name, 'FULL')
    ExtMetShortNameList = {'glu', 'pyr'};

else
    error('Unknown model');
end


%% In the FULL model, set external metabolite concentrations (in mM). In other models, check that the FULL file exists
if strcmp(ModelSpecs.Name, 'FULL')
    
    if strcmp(Glu, 'Inf')
        ModelSpecs.ExtMetConc.glu = Inf;
        ModelSpecs.ExtMetConc.pyr = 1;
        %%%  Taken from Rohwer et al JBC 2000 (Table I), measured at 22 mM of
        %%%  external glucose (effectively infinite, given that the Km of PTS ~ 0.02 mM, see below);
        %%%  also 0.390 mM from https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=7&id=101192,

    elseif strcmp(Glu, 'Med')
        ModelSpecs.ExtMetConc.glu = 0.02; %0.0556 ;
        ModelSpecs.ExtMetConc.pyr = 0.1;
        
    elseif strcmp(Glu, 'Low')

        ModelSpecs.ExtMetConc.glu = 0.002;
        ModelSpecs.ExtMetConc.pyr = 0.01;
    end
    
else
    filename = sprintf('%s/model_FULL_%s.mat', data_dir, Glu);
    
    if exist(filename, 'file')
        load(filename, 'WT');
        WTFULL = WT;
        clear WT;
    else
        error('Before generating model %s, need to obtain external metabolite concentrations from FULL model. Please, initialize the FULL model\n',...
            ModelSpecs.Name);
    end
end



%% Designate output fluxes
%%% Lower glycolysis:
if strcmp(ModelSpecs.Name, 'LG') 
    ModelSpecs.rxnOutEnzymes = {'ENO'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1];

%%% G6PDH and PGDH:    
elseif strcmp(ModelSpecs.Name, 'PPSMALL') 
    ModelSpecs.rxnOutEnzymes = {'PGDH'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1]';        

%%% PP:    
elseif strcmp(ModelSpecs.Name, 'PP') 
    ModelSpecs.rxnOutEnzymes = {'TKb', 'TKa', 'TA'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1, +1, -1]';

%%% Upper glycolysis and PP:    
elseif strcmp(ModelSpecs.Name, 'UGPP')
    ModelSpecs.rxnOutEnzymes = {'ALDO', 'TIS', 'TKb', 'TKa', 'TA'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1, +1, +1, +1, -1]';

%%% Upper and lower glycolysis and PP:    
elseif strcmp(ModelSpecs.Name, 'GPP')
    ModelSpecs.rxnOutEnzymes = {'ENO'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1];
    
%%% FULL:
elseif strcmp(ModelSpecs.Name, 'FULL')
    ModelSpecs.rxnOutEnzymes = {'PK', 'PTS'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1 +1]';
else
    error('Model %s unknown', ModelSpecs.Name);
end




%% Load model:
OrigModelObj = sbmlimport(modelfile);

%%% ===== Modify model ======

%% 1. Make the co-metabolite (ATP, ADP, etc) concentrations constant:
for irule = 1:length(OrigModelObj.Rules)
    delete(OrigModelObj.Rules(1));
end
clear irule;

%% 2. Remove dilution by growth:
for irxn = 48:-1:31
    delete(OrigModelObj.Reaction(irxn));
end
clear irxn;

%% 3. Remove extreneous rxns:
if strcmp(ModelSpecs.Name, 'PGMENO')
    RemoveRxnsIX = sort([1 3 30 29 9 14 13 16 25 24 19:23 2 5 10 12 4 26:28 6:8 11 15], 'descend');
elseif strcmp(ModelSpecs.Name, 'LG')
    RemoveRxnsIX = sort([1 3 30 29 9 14 13 16 25 24 19:23 2 5 10 12 4 26:28 6:8], 'descend');    
elseif strcmp(ModelSpecs.Name, 'PPSMALL')
    RemoveRxnsIX = sort( setdiff(1:30, [4 26:28]) , 'descend');
elseif strcmp(ModelSpecs.Name, 'PP')    
    RemoveRxnsIX = sort([1 3 30 29 9 14 13 16 25 24 19:23, 2 5 10 12 11 15 17 18], 'descend');
elseif strcmp(ModelSpecs.Name, 'UGPP')
    RemoveRxnsIX = sort([1 3 30 29 9 14 13 16 25 24 19:23 11 15 17 18], 'descend');
elseif strcmp(ModelSpecs.Name, 'GPP')
    RemoveRxnsIX = sort([1 3 30 29 9 14 13 16 25 24 19:23], 'descend');
elseif strcmp(ModelSpecs.Name, 'FULL')
    RemoveRxnsIX = sort([3 30 9 13 14 16 21:25 29], 'descend');
end

for irxn = RemoveRxnsIX
    delete(OrigModelObj.Reaction(irxn));
end
clear irxn RemoveRxnsIX;


%% 4. Remove extreneous metabolites:
if strcmp(ModelSpecs.Name, 'PGMENO')
    RemoveMetIX = sort([1 4 6, 3 5 7 8 9 11 12 13 14 18 10 15], 'descend');    
elseif strcmp(ModelSpecs.Name, 'LG')
    RemoveMetIX = sort([1 4 6, 3 5 7 8 9 11 12 13 14 18], 'descend');
elseif strcmp(ModelSpecs.Name, 'PPSMALL')    
    RemoveMetIX = sort(setdiff(1:18, [3, 7, 18 12 13]), 'descend');
elseif strcmp(ModelSpecs.Name, 'PP')    
    RemoveMetIX = sort([1 4 6, 2 8 14:17], 'descend');
elseif strcmp(ModelSpecs.Name, 'UGPP')
    RemoveMetIX = sort([1 4 6 15:17], 'descend');    
elseif strcmp(ModelSpecs.Name, 'GPP')
    RemoveMetIX = sort([1 4 6], 'descend');
elseif strcmp(ModelSpecs.Name, 'FULL')
    RemoveMetIX = sort([6], 'descend');
end

for isp = RemoveMetIX
    delete(OrigModelObj.Species(isp));
end
clear isp RemoveMetIX;


%% 5. In all models other than FULL, remove extracellular compartment and set concentrations of external metabolites
if ~strcmp(ModelSpecs.Name, 'FULL')
    
    %%% Delete extracellular compartment
    set(OrigModelObj.Compartments(2), 'Owner', []);
    delete(OrigModelObj.Compartments(1));        
    
    %%% Set initial amounts of external metabolites:    
    for imet = 1:length( ExtMetShortNameList )
        CurrMetShortName = ExtMetShortNameList{ imet };
        CurrMetFullName = get_full_name( ModelSpecs.MetNames, CurrMetShortName) ;
        ixFULL = find( strcmp( get(WTFULL.m.Species, 'Name'), CurrMetFullName ) );

        if isempty(ixFULL)
            error('External metabolite %s (%s) not found in the FULL model', CurrMetShortName, CurrMetFullName);
        end

        CurrMetConc = get(WTFULL.m.Species(ixFULL), 'InitialAmount');
        ModelSpecs.ExtMetConc.( CurrMetShortName ) = CurrMetConc;
        
        ix = find( strcmp( get(OrigModelObj.Species, 'Name'), CurrMetFullName ) );
        if isempty(ix)
            error('Unknown external metabolite');
        end

        set(OrigModelObj.Species(ix), 'InitialAmount', CurrMetConc );
        set(OrigModelObj.Species(ix), 'ConstantAmount', true);
    end
    clear ExtMetShortNameList imet ix;
end


%% 6. In the FULL model, change stoichiometry of the PTS rxn, value of the Km in the PTS system and PTS kinetics if [Glu] = Inf
if strcmp(ModelSpecs.Name, 'FULL')
    
    %%% Change the stoichiometry of the PTS rxn:
    set(OrigModelObj.Reactions(1), 'Stoichiometry', [-1 -1 1 1]);
    
    set(OrigModelObj.Reactions(1).KineticLaw.Parameters(2), 'Value', 0.02);
    %%% The value of Km for the PTS system seems to be unrealistically
    %%% large in the original model (~3 M). But kinetic studies show that
    %%% it is closer to 20 然 (see Stock et al, JBC 1982) or even lower (5
    %%% 然, see Natarajan, Srienc, Metabolic Eng (1999))
        
    %%% Set initial amounts of external metabolites:
    ExtMetShortNameList = fieldnames(ModelSpecs.ExtMetConc);
    for imet = 1:length( ExtMetShortNameList )
        CurrMetShortName = ExtMetShortNameList{ imet };
        CurrMetFullName = get_full_name( ModelSpecs.MetNames,  CurrMetShortName) ;
        CurrMetConc = ModelSpecs.ExtMetConc.( CurrMetShortName );
        
        ix = find( strcmp( get(OrigModelObj.Species, 'Name'), CurrMetFullName ) );
        if isempty(ix)
            error('Unknown external metabolite');
        end
        
        if strcmp(Glu, 'Inf') && strcmp(CurrMetShortName, 'glu')
            continue;
        end
        
        set(OrigModelObj.Species(ix), 'InitialAmount', CurrMetConc );
        set(OrigModelObj.Species(ix), 'ConstantAmount', true);
    end
    clear ExtMetShortNameList imet ix;
    
    %%% If [Glu] = Inf, re-write PTS kinetics and delete extracellular compartment
    if strcmp(Glu, 'Inf')
        set(OrigModelObj.Reactions(1), 'Reaction', 'cytosol.[Phosphoenol pyruvate] <-> cytosol.[Glucose-6-Phosphate] + cytosol.Pyruvate');
        set(OrigModelObj.Reactions(1), 'ReactionRate',...
            'rmaxPTS*(cytosol.[Phosphoenol pyruvate]/cytosol.Pyruvate)/((KPTSa3+(cytosol.[Phosphoenol pyruvate]/cytosol.Pyruvate))*(1+power(cytosol.[Glucose-6-Phosphate],nPTSg6p)/KPTSg6p))');
        
        %%% Change stoichoimetry to reflect deletion of glucose
        set(OrigModelObj.Reactions(1), 'Stoichiometry', [-1 1 1]);
        
        %%% Delete glucose
        delete(OrigModelObj.Species(1));
        
        %%% Delete extracellular compartment
        set(OrigModelObj.Compartments(2), 'Owner', []);
        delete(OrigModelObj.Compartments(1));
    end
    
end


% We will monitor this flux as the output
RxnNames = get(OrigModelObj.Reactions, 'Name');
ModelSpecs.rxnOutIX = nan * ModelSpecs.rxnOutEnzymeCoeff;
for irxn = 1:length( ModelSpecs.rxnOutEnzymes )
    CurrShortRxnName = ModelSpecs.rxnOutEnzymes{ irxn };
    CurrLongRxnName = get_full_name( ModelSpecs.EnzNames, CurrShortRxnName ) ;
    ModelSpecs.rxnOutIX(irxn) = find( strcmp(RxnNames , CurrLongRxnName) );
end



%% Find wildtype steady state
[WT.m, WT.FluxDistr] = getMutFlux(OrigModelObj);
WT.OutFlux = sum( ModelSpecs.rxnOutEnzymeCoeff .* WT.FluxDistr.Flux(ModelSpecs.rxnOutIX) );


%% Calculate control coefficients and the diagonal interaction coefficients for all rxns:
fprintf('Calculating control coefficients ...\n');
[WT.FluxDistr.FCC, WT.FluxDistr.FIC] = getFCC( WT.m, ModelSpecs, WT.OutFlux);
save(sprintf('%s/model_%s_%s.mat', data_dir, ModelSpecs.Name, Glu), 'WT', 'ModelSpecs');


%% Calculate the non-diagonal interaction coefficients for all rxns:
load(sprintf('%s/model_%s_%s.mat', data_dir, ModelSpecs.Name, Glu), 'WT', 'ModelSpecs');
fprintf('Calculating interaction coefficients ...\n');
WT.FluxDistr.FIC = getFIC( WT.m, ModelSpecs, WT.OutFlux, WT.FluxDistr.FCC, WT.FluxDistr.FIC);

save(sprintf('%s/model_%s_%s.mat', data_dir, ModelSpecs.Name, Glu), 'WT', 'ModelSpecs');


