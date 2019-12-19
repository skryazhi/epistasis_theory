function [ModelL, FluxDistrL] = getMutFlux(OrigModelObj, mutList)
%
% getMutFlux calculates the steady-state fluxes for a list mutants
%
% [ModelL, FluxDistrL] = getMutFlux(OrigModelObj) takes the simbiology
% model OrigModelObj, finds the steady state, calculates the steady-state
% flux distribution. ModelL is the perturbed model at steady-state,
% FluxDistrL is the resulting flux distribution.
% 
% ... = getMutFlux(OrigModelObj, []) is equivalent to
% getMutFlux(OrigModelObj)
%
% ... = getMutFlux(OrigModelObj, mutList) perturbes the
% original model it accoring to perturbations given in mutList, finds the
% steady states for these perturbed models, calculates the steady-state
% fluxes. ModelL is the cell array of perturbed models at steady-state,
% FluxDistrL is the cell array of resulting flux distributions.
%
% mutList is an M by 2 cell array where each row corresponds to a mutant;
% mutList{irow, 1} is an N by 1 cell array of shorthand names the enzymes
% being pertubed by mutations in this mutant; mutList{irow, 2} is an N by 1
% cell array of long name of these enzymes; mutList{irow,2} is a N by 1
% array of perturbations. Each perturbation is the relative effect of the
% mutation on that rMAX of the corresponding reaction. EXAMPLE: one row of
% mutList might look like this: {{'PGI'}, [0.8]}. This means that rMAX of
% PGI in the mutant will be 80% of the value of the wildtype.

warning('off', 'all');

if nargin < 2 || isempty(mutList)
    [success, variant_out, ModelL] = sbiosteadystate(OrigModelObj, 'MaxStopTime', 1e7, 'RelTol', 1e-8, 'AbsTol', 1e-10);
    if ~success
        fprintf('Failed to converge to steady state!\n');
    end
    FluxDistrL = get_fluxes(ModelL, 'v');
    return;
end

RxnNames = arrayfun(@(x) get(x, 'Name'), OrigModelObj.Reactions,...
    'UniformOutput', false);

nMut = size( mutList, 1 );
ModelL = cell( nMut, 1);
FluxDistrL = cell( nMut, 1);

for imut = 1:nMut
    
    % fprintf('Mut #%d out of %d (%.1f%%)\n', imut, nMut, imut/nMut*100 );
    
    % which reactions/parameters are affected in this mutant?
    EnzList = mutList{imut,1};
    PertList = mutList{imut,2};
    EnzIdList = cellfun(@(EnzName) find( strcmp(RxnNames, EnzName) ), EnzList);
    nEnz = length( mutList{imut,1} );
    
    MutModelObj = copyobj( OrigModelObj );
    
    for iEnz = 1:nEnz
        EnzId = EnzIdList(iEnz);
        x = get( MutModelObj.Reactions(EnzId).KineticLaw.Parameters(1), 'Value');
        x = x * PertList(iEnz);
        set(MutModelObj.Reactions(EnzId).KineticLaw.Parameters(1), 'Value', x  );

%         fprintf('\t%s\t%.2g\n',...
%             get( OrigModelObj.Reactions(EnzId).KineticLaw.Parameters(1), 'Name'),...
%             PertList(iEnz) );
    end
    
    [success, variant_out,  ModelL{imut}] = sbiosteadystate(MutModelObj, 'MaxStopTime', 1e7, 'RelTol', 1e-8, 'AbsTol', 1e-10);
    if ~success
        fprintf('Failed to converge to steady state!\n');
    end

    FluxDistrL{imut} = get_fluxes(ModelL{imut});    
end

warning('on');

