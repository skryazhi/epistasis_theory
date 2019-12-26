function [MutModel, FluxDistr] = getMutFlux(OrigModelObj, mutRxnNames, mutPerts)
%
% getMutFlux calculates the steady-state fluxes for a mutant
%
% [MutModel, FluxDistr] = getMutFlux(OrigModelObj) takes the simbiology
% model OrigModelObj, finds the steady state, calculates the steady-state
% flux distribution. MutModel is the model at steady-state (i.e.,
% wildtype), FluxDistr is the resulting flux distribution.
%
% ... = getMutFlux(OrigModelObj, mutRxnNames, mutPerts) perturbes the
% original model it accoring to the perturbation given in mutRxnNames and
% mutPerts, finds the steady state for the perturbed model, calculates the
% steady-state fluxes. MutModel is the perturbed model at steady-state;
% FluxDistr is the resulting flux distribution
%
% mutRxnNames is an N by 1 cell array of full names of reacations being
% pertubed by the mutation. mutPerts is a N by 1 array of perturbations.
% Each perturbation is the relative effect of the mutation on that rMAX of
% the corresponding reaction. EXAMPLE: mutRxnNames = {'PGI'}, mutPerts =
% 0.8. This means that rMAX of PGI in the mutant will be 80% of the value
% of the wildtype.

warning('off', 'all');

RelTol = 1e-8;
AbsTol = 1e-10;
MaxStop = 1e7; 

if nargin < 2
    
    [success, variant_out, MutModel] = sbiosteadystate(OrigModelObj,...
        'MaxStopTime', MaxStop, 'RelTol', RelTol, 'AbsTol', AbsTol);
    if ~success
        fprintf('Failed to converge to steady state!\n');
    end
    FluxDistr = get_fluxes(MutModel);
    return;
end

RxnNames = get( OrigModelObj.Reactions, 'Name');

PertRxnIX = cellfun(@(EnzName) find( strcmp(RxnNames, EnzName) ), mutRxnNames);
nPerRxn = length( mutRxnNames );

MutModelObj = copyobj( OrigModelObj );

for iRxn = 1:nPerRxn
    x = get( MutModelObj.Reactions( PertRxnIX(iRxn)  ).KineticLaw.Parameters(1), 'Value');
    y = x * mutPerts(iRxn);
    
    set(MutModelObj.Reactions( PertRxnIX(iRxn) ).KineticLaw.Parameters(1), 'Value', y  );    
end

[success, variant_out,  MutModel] = sbiosteadystate(MutModelObj,...
    'MaxStopTime', MaxStop, 'RelTol', RelTol, 'AbsTol', AbsTol);

if ~success
    fprintf('Failed to converge to steady state!\n');
end

FluxDistr = get_fluxes(MutModel);

warning('on');

