function rxn = get_fluxes(mObj, verbose)
% GET_FLUXES calculates reaction fluxes at the initial state of the model
%
% RXN = GET_FLUXES( MOBJ, '-v' ) takes the SimBio object MOBJ, retreaves
% the Reaction attribute of each reaction as a string and replaces all
% values in this string with the model's parameter values and initial
% species concentrations, then evaluates this expression. If the optional
% tag '-v' is given, the verbose version is called. RXN is a structure with
% two fields. RXN.NAME is an n by 1 cell array containing the names of
% reactions. RXN.FLUX is an n by 1 vector containing the values of the
% reaction rates.

IFVERBOSE = false;

if nargin > 1
    if strcmp(verbose, 'v')
        IFVERBOSE = true;
    end
end

SpeciesDict = getDict(mObj.Species);
ParamDict = getDict(mObj.Parameters);

nRxn = length(mObj.Reactions);

rxn.Flux = nan(nRxn,1);
rxn.Name = cell(nRxn,1);

for irxn = 1:nRxn
    
    CurrRxn = mObj.Reactions(irxn);    
    
    rxn.Name{irxn} = CurrRxn.Name;
    
    RxnParamDict = getDict(CurrRxn.KineticLaw.Parameters);
    
    RxnRateString = CurrRxn.ReactionRate;

    % Replacing species concentrations with their current values:
    for isp = 1:length(SpeciesDict.Key)
        RxnRateString = strrep(RxnRateString, SpeciesDict.Key{isp}, SpeciesDict.Value{isp});
    end

    % Replacing parameters of the current reaction with their current values:
    for ip = 1:length(RxnParamDict.Key)
        RxnRateString = strrep(RxnRateString, RxnParamDict.Key{ip}, RxnParamDict.Value{ip});
    end

    % Replacing global parameters with their current values:
    for ip = 1:length(ParamDict.Key)
        RxnRateString = strrep(RxnRateString, ParamDict.Key{ip}, ParamDict.Value{ip});
    end

    if length(mObj.Compartments) == 2
        RxnRateString = strrep(RxnRateString, 'extracellular', sprintf('%.3e', mObj.Compartments(1).Capacity));
        RxnRateString = strrep(RxnRateString, 'cytosol', sprintf('%.3e', mObj.Compartments(2).Capacity));
    elseif length(mObj.Compartments) == 1
        RxnRateString = strrep(RxnRateString, 'cytosol', sprintf('%.3e', mObj.Compartments(1).Capacity));
    end
        
            
    rxn.Flux(irxn) = eval( RxnRateString );
    
    if IFVERBOSE
        fprintf('rnx #%d, %s:\n', irxn, rxn.Name{irxn});
        fprintf('\t%s\n', RxnRateString);
        fprintf('\t%s', RxnRateString);
        fprintf(' = %.3g\n', rxn.Flux(irxn));
    end
end

%% Checking that the fluxes actually balance for each metabolite:
for isp = 1:length(mObj.Species)
    spName = get(mObj.Species(isp), 'Name');
    
    % Find rxns in which this species participates as a substrate:
    substrList = [];
    for irxn = 1:length(mObj.Reactions)
        for isubstr = 1:length(mObj.Reactions(irxn).Reactants)
            if strcmp( get(mObj.Reactions(irxn).Reactants(isubstr), 'Name'), spName )
               substrList = [substrList irxn];
               break;
            end
        end
    end
    
    % Find rxns in which this species participates as a product:
    prodList = [];
    for irxn = 1:length(mObj.Reactions)
        for iprod = 1:length(mObj.Reactions(irxn).Products)
            if strcmp( get(mObj.Reactions(irxn).Products(iprod), 'Name'), spName )
                prodList = [prodList irxn];
                break;
            end
        end
    end
    
    x = sum( rxn.Flux( prodList  ) ) - sum( rxn.Flux( substrList ) );

    if IFVERBOSE
        fprintf('Species #%d (%s) derivative = %.3g\n', isp, get(mObj.Species(isp), 'Name'), x );
    end
end







function dict = getDict(ObjL)
% Create a dictionary for the list of objects ObjL (parameters or species).
% Dictionary dict is a structure with two fields. dict.Key contains the
% name of the object, dict.Value contains the value of the parameter or
% species.

nObj = length( ObjL );

dict.Key = cell(nObj,1);
dict.Value = cell(nObj,1);

for isp = 1:nObj
    CurrObj = ObjL(isp);

    if strcmp(CurrObj.Type, 'species')
        if isempty(regexp(CurrObj.Name, '[\d -]', 'once'))
            dict.Key{isp} = sprintf('%s.%s', CurrObj.Parent.Name, CurrObj.Name);
        else
            dict.Key{isp} = sprintf('%s.[%s]', CurrObj.Parent.Name, CurrObj.Name);
        end
        dict.Value{isp} = sprintf('%.16e', CurrObj.InitialAmount);
    elseif strcmp(CurrObj.Type, 'parameter')
        dict.Key{isp} = sprintf('%s', CurrObj.Name);
        dict.Value{isp} = sprintf('%.16e', CurrObj.Value);
    end
end
[tmp, ix] = sort( cellfun(@length, dict.Key), 'descend' );
dict.Key = dict.Key(ix);
dict.Value = dict.Value(ix);

