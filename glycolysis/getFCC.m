function [FCC, FICD] = getFCC(ModelObj, ModelSpecs, WTOutFlux)
%
% getFCC calculates the control coefficients in the model
%
% [FCC, FICD] = getFCC(ModelObj, ModelSpecs, WTOutFlux) takes the
% simbiology model ModelObj and calculates the control coefficients and the
% diagonal interaction coefficients of all reactions
% 
% INPUT VARIABLES:
%
% ModelObj is a simbiology kinetic model of metabolism with initial
% concentrations of metabolites at steady state
% 
% ModelSpecs is a structure that must have three fields:
% ModelSpecs.rxnOutIX is the vector of indices of the reactions that are
% considered as output fluxes;  ModelSpecs.rxnOutEnzymeCoeff are the
% coefficients with which these fluxes should be summed up to get the total
% output flux. ModelSpecs.EnzNames is a cell array with two columns: first
% column has short enzyme names; second column has full enzyme names.
%
% WTOutFlux specifies the value of the output flux for the wildtype
% (unperturbed) model, if it was previously calculated (which saves time on
% this calculation)
%
% OUTPUT VARIABELS:
%
% FCC is a vector of control coefficients corresponding to
% ModelObj.Reactions array
%
% FICD is the diagnoal matrix of flux interaction coefficients, i.e., d^2
% F/ dv^2 (\delta v)^2/F, where F is flux and v is the reaction rate

RxnNames = arrayfun(@(x) get(x, 'Name'), ModelObj.Reactions, 'UniformOutput', false);

AbsTol = 1e-4;

nRxn = length(RxnNames);
FCC = zeros(nRxn,1);
FICD = zeros( nRxn, nRxn );

nPoints = 10;
minPert = 0.75;
maxPert = 1.25;

PertVec = linspace(minPert, maxPert, nPoints)';
FluxVec = nan(nPoints, 1);

for irxn = 1:nRxn
    CurrFullRxnName = RxnNames{irxn};
    
    fprintf('%s: ', CurrFullRxnName );
    
    
    for ipoint = 1:nPoints
        
        CurrPert = PertVec(ipoint);
        
        [MutModel, FluxDistr] = getMutFlux(ModelObj, {CurrFullRxnName}, CurrPert);
        FluxVec(ipoint) = sum( ModelSpecs.rxnOutEnzymeCoeff .* FluxDistr.Flux( ModelSpecs.rxnOutIX) );
    end
    
    Y = FluxVec/WTOutFlux-1;
    X1 = PertVec-1;
    X2 = 1/2*X1.^2;
   
    
    b = regress(Y, [X1 , X2] );

    %%% Diagonose problems by plotting flux vs perturbation
    % clf; hold on; box on;
    %     plot( X1, Y, 'ok');
    %     plot( X1, b(1)*X1 + b(2)*X2, 'k-', 'LineWidth', 2);
    
    if abs(b(1)) > AbsTol
        FCC(irxn) = b(1);
    end
    
    if abs(b(2)) > AbsTol
        FICD( irxn, irxn ) = b(2);
    end
    
    fprintf('FCC = %.2f%%, FIC = %.4f \n', FCC(irxn)*100, FICD(irxn,irxn) );
end