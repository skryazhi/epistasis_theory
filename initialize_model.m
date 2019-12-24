function [WT, ModelSpecs] = initialize_model(curr_dir, Glu, ModelSpecs, MetNames, EnzNames)


%% Designating external metabolites:
%%% Lower glycolysis:
if strcmp(ModelSpecs.Type, 'LG') 
    ExtMetShortNameList = {'gap', 'pep'};
        
%%% G6PDH and PGDH:    
elseif strcmp(ModelSpecs.Type, 'PPPSMALL') 
    ExtMetShortNameList = {'g6p', 'ribu5p'};
    
%%% PPP:    
elseif strcmp(ModelSpecs.Type, 'PPP') 
    ExtMetShortNameList = {'g6p', 'f6p', 'gap'};

%%% Upper glycolysis and PPP:    
elseif strcmp(ModelSpecs.Type, 'UGPPP')
    ExtMetShortNameList = {'g6p', 'gap', 'pep'};
    
%%% Upper and lower glycolysis and PPP:    
elseif strcmp(ModelSpecs.Type, 'GPPP')
        ExtMetShortNameList = {'g6p', 'pep'};
    
%%% Full model    
elseif strcmp(ModelSpecs.Type, 'FULL')
    ExtMetShortNameList = {'glu', 'pyr'};

else
    error('Unknown model');
end


%% In the FULL model, set external metabolite concentrations (in mM). In other models, checke that the FULL file exists
if strcmp(ModelSpecs.Type, 'FULL')
    
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
    filename = sprintf('%s/model_FULL_%s.mat', curr_dir, Glu);
    
    if exist(filename, 'file')
        load(filename, 'WT');
        WTFULL = WT;
        clear WT;
    else
        error('Before generating model %s, need to obtain external metabolite concentrations from FULL model. Please, initialize the FULL model\n',...
            ModelSpecs.Type);
    end
end



%% Designate output fluxes
%%% Lower glycolysis:
if strcmp(ModelSpecs.Type, 'LG') 
    ModelSpecs.rxnOutEnzymes = {'ENO'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1];

%%% G6PDH and PGDH:    
elseif strcmp(ModelSpecs.Type, 'PPPSMALL') 
    ModelSpecs.rxnOutEnzymes = {'PGDH'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1]';        

%%% PPP:    
elseif strcmp(ModelSpecs.Type, 'PPP') 
    ModelSpecs.rxnOutEnzymes = {'TKb', 'TKa', 'TA'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1, +1, -1]';

%%% Upper glycolysis and PPP:    
elseif strcmp(ModelSpecs.Type, 'UGPPP')
    ModelSpecs.rxnOutEnzymes = {'ALDO', 'TIS', 'TKb', 'TKa', 'TA'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1, +1, +1, +1, -1]';

%%% Upper and lower glycolysis and PPP:    
elseif strcmp(ModelSpecs.Type, 'GPPP')
    ModelSpecs.rxnOutEnzymes = {'ENO'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1];
    
%%% FULL:
elseif strcmp(ModelSpecs.Type, 'FULL')
    ModelSpecs.rxnOutEnzymes = {'PK', 'PTS'};
    ModelSpecs.rxnOutEnzymeCoeff = [+1 +1]';
else
    error('Model %s unknown', ModelSpecs.Type);
end




%% Load model:
filename = sprintf('%s/BIOMD0000000051.xml', curr_dir);

OrigModelObj = sbmlimport(filename);


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
if strcmp(ModelSpecs.Type, 'PGMENO')
    RemoveRxnsIX = sort([1 3 30 29 9 14 13 16 25 24 19:23 2 5 10 12 4 26:28 6:8 11 15], 'descend');
elseif strcmp(ModelSpecs.Type, 'LG')
    RemoveRxnsIX = sort([1 3 30 29 9 14 13 16 25 24 19:23 2 5 10 12 4 26:28 6:8], 'descend');    
elseif strcmp(ModelSpecs.Type, 'PPPSMALL')
    RemoveRxnsIX = sort( setdiff(1:30, [4 26:28]) , 'descend');
elseif strcmp(ModelSpecs.Type, 'PPP')    
    RemoveRxnsIX = sort([1 3 30 29 9 14 13 16 25 24 19:23, 2 5 10 12 11 15 17 18], 'descend');
elseif strcmp(ModelSpecs.Type, 'UGPPP')
    RemoveRxnsIX = sort([1 3 30 29 9 14 13 16 25 24 19:23 11 15 17 18], 'descend');
elseif strcmp(ModelSpecs.Type, 'GPPP')
    RemoveRxnsIX = sort([1 3 30 29 9 14 13 16 25 24 19:23], 'descend');
elseif strcmp(ModelSpecs.Type, 'FULL')
    RemoveRxnsIX = sort([3 30 9 13 14 16 21:25 29], 'descend');
end

for irxn = RemoveRxnsIX
    delete(OrigModelObj.Reaction(irxn));
end
clear irxn RemoveRxnsIX;


%% 4. Remove extreneous metabolites:
if strcmp(ModelSpecs.Type, 'PGMENO')
    RemoveMetIX = sort([1 4 6, 3 5 7 8 9 11 12 13 14 18 10 15], 'descend');    
elseif strcmp(ModelSpecs.Type, 'LG')
    RemoveMetIX = sort([1 4 6, 3 5 7 8 9 11 12 13 14 18], 'descend');
elseif strcmp(ModelSpecs.Type, 'PPPSMALL')    
    RemoveMetIX = sort(setdiff(1:18, [3, 7, 18 12 13]), 'descend');
elseif strcmp(ModelSpecs.Type, 'PPP')    
    RemoveMetIX = sort([1 4 6, 2 8 14:17], 'descend');
elseif strcmp(ModelSpecs.Type, 'UGPPP')
    RemoveMetIX = sort([1 4 6 15:17], 'descend');    
elseif strcmp(ModelSpecs.Type, 'GPPP')
    RemoveMetIX = sort([1 4 6], 'descend');
elseif strcmp(ModelSpecs.Type, 'FULL')
    RemoveMetIX = sort([6], 'descend');
end

for isp = RemoveMetIX
    delete(OrigModelObj.Species(isp));
end
clear isp RemoveMetIX;


%% 5. In all models other than FULL, remove extracellular compartment and set concentrations of external metabolites
if ~strcmp(ModelSpecs.Type, 'FULL')
    
    %%% Delete extracellular compartment
    set(OrigModelObj.Compartments(2), 'Owner', []);
    delete(OrigModelObj.Compartments(1));        
    
    %%% Set initial amounts of external metabolites:    
    for imet = 1:length( ExtMetShortNameList )
        CurrMetShortName = ExtMetShortNameList{ imet };
        CurrMetFullName = MetNames.( CurrMetShortName );
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
if strcmp(ModelSpecs.Type, 'FULL')
    
    %%% Change the stoichiometry of the PTS rxn:
    set(OrigModelObj.Reactions(1), 'Stoichiometry', [-1 -1 1 1]);
    
    set(OrigModelObj.Reactions(1).KineticLaw.Parameters(2), 'Value', 0.02);
    %%% The value of Km for the PTS system seems to be unrealistically
    %%% large in the original model (~3 M). But kinetic studies show that
    %%% it is closer to 20 µM (see Stock et al, JBC 1982) or even lower (5
    %%% µM, see Natarajan, Srienc, Metabolic Eng (1999))
        
    %%% Set initial amounts of external metabolites:
    ExtMetShortNameList = fieldnames(ModelSpecs.ExtMetConc);
    for imet = 1:length( ExtMetShortNameList )
        CurrMetShortName = ExtMetShortNameList{ imet };
        CurrMetFullName = MetNames.( CurrMetShortName );
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



%     %%% 6. Set all coefficients to 1
%     for irxn = 1:length(OrigModelObj.Reactions)
%         for iparam = 1:length(OrigModelObj.Reactions(irxn).KineticLaw.Parameters)
%             set(OrigModelObj.Reactions(irxn).KineticLaw.Parameters(iparam), 'Value', 1);
%         end
%     end


% We will monitor this flux as the output
RxnNames = get(OrigModelObj.Reactions, 'Name');
rxnOutIdVec = cellfun(@(EnzShortName) find(strcmp(RxnNames,EnzNames.(EnzShortName))) ,ModelSpecs.rxnOutEnzymes);


%% Find wildtype steady state

[WT.m, WT.FluxDistr] = getMutFlux(OrigModelObj, []);
WT.OutFlux = sum( ModelSpecs.rxnOutEnzymeCoeff .* WT.FluxDistr.Flux(rxnOutIdVec) );


%% Calculate control coefficients for all rxns:

WT.FluxDistr.CtrlCoeff = getCtrlCoeff( WT.m, ModelSpecs, rxnOutIdVec, WT.OutFlux);

save(sprintf('%s/model_%s_%s.mat', curr_dir, ModelSpecs.Type, Glu), 'WT', 'ModelSpecs');




