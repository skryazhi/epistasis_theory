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

    RxnRateString = strrep(RxnRateString, 'extracellular', sprintf('%.3e', mObj.Compartments(1).Capacity));
    RxnRateString = strrep(RxnRateString, 'cytosol', sprintf('%.3e', mObj.Compartments(2).Capacity));
            
    rxn.Flux(irxn) = eval( RxnRateString );
    
    if IFVERBOSE
        fprintf('rnx #%d, %s:\n', irxn, rxn.Name{irxn});
        fprintf('\t%s\n', RxnRateString);
        fprintf('\t%s', RxnRateString);
        fprintf(' = %.3g\n', rxn.Flux(irxn));
    end
end

%% Checking that the fluxes actually balance for each metabolite:
% pep
isp = 2;
species_rxn.plus{isp} = [18]; species_rxn.minus{isp} = [1 19 20 21 23 39];
% g6p
isp = 3;
species_rxn.plus{isp} = [1]; species_rxn.minus{isp} = [2 3 4 31];
% pyr
isp = 4;
species_rxn.plus{isp} = [1 13 19 25]; species_rxn.minus{isp} = [22 24 44];
% f6p
isp = 5;
species_rxn.plus{isp} = [2 6 8]; species_rxn.minus{isp} = [5 9 32]; % MurSynt (Rxn #9) has stoichiometry coefficient 2!
% g1p
isp = 6;
species_rxn.plus{isp} = [3]; species_rxn.minus{isp} = [30 47];
% 6pg
isp = 7;
species_rxn.plus{isp} = [4]; species_rxn.minus{isp} = [26 45];
% fdp
isp = 8;
species_rxn.plus{isp} = [5]; species_rxn.minus{isp} = [10 33];
% s7p
isp = 9;
species_rxn.plus{isp} = [7]; species_rxn.minus{isp} = [6 43];
% gap
isp = 10;
species_rxn.plus{isp} = [7 8 10 12 13]; species_rxn.minus{isp} = [6 11 34];
% e4p
isp = 11;
species_rxn.plus{isp} = [6]; species_rxn.minus{isp} = [8 23 46];
% xyl5p
isp = 12;
species_rxn.plus{isp} = [28]; species_rxn.minus{isp} = [7 8 42];
% rib5p
isp = 13;
species_rxn.plus{isp} = [27]; species_rxn.minus{isp} = [7 29 41];
% dhap
isp = 14;
species_rxn.plus{isp} = [10]; species_rxn.minus{isp} = [12 14 35];
% pgp
isp = 15;
species_rxn.plus{isp} = [11]; species_rxn.minus{isp} = [15 36];
% 3pg
isp = 16;
species_rxn.plus{isp} = [15]; species_rxn.minus{isp} = [16 17 37];
% 2pg
isp = 17;
species_rxn.plus{isp} = [17]; species_rxn.minus{isp} = [18 38];
% ribu5p
isp = 18;
species_rxn.plus{isp} = [26]; species_rxn.minus{isp} = [27 28 40];

for isp = 2:length(mObj.Species)
    x = sum( rxn.Flux( species_rxn.plus{isp} ) ) - sum( rxn.Flux( species_rxn.minus{isp} ) );
    if isp == 5
        x = x - rxn.Flux(9);
    end
    if IFVERBOSE
        fprintf('Species #%d (%s) derivative = %.3g\n', isp, get(mObj.Species(isp), 'Name'), x );
    end
end
clear isp x;




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
        dict.Value{isp} = sprintf('%.6e', CurrObj.InitialAmount);
    elseif strcmp(CurrObj.Type, 'parameter')
        dict.Key{isp} = sprintf('%s', CurrObj.Name);
        dict.Value{isp} = sprintf('%.6e', CurrObj.Value);
    end
end
[tmp, ix] = sort( cellfun(@length, dict.Key), 'descend' );
dict.Key = dict.Key(ix);
dict.Value = dict.Value(ix);

