function FIC = getFIC(ModelObj, ModelSpecs, WTOutFlux, FCC, FICD)
%
% getFIC calculates the non-diagonal flux interaction coefficients
%
% FIC = getFIC(ModelObj, ModelSpecs, WTOutFlux, FCC, FICD) takes the
% simbiology model ModelObj and calculates the non-diagonal control
% coefficients for all reaction pairs
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
% FCC is a vector of control coefficients corresponding to
% ModelObj.Reactions array
%
% FICD is the diagnoal matrix of flux interaction coefficients
%
% OUTPUT VARIABELS:
%
% FIC is the full matrix of the flux interaction coefficients. The
% diagonal of this matrix is equal to that of FICD

RxnNames = arrayfun(@(x) get(x, 'Name'), ModelObj.Reactions, 'UniformOutput', false);

AbsTol = 1e-4;

nRxn = length(RxnNames);
FIC = FICD;

nPert = 4;
nPoints = nPert^2 ;
minPert = 0.75;
maxPert = 1.25;

PertVec = linspace(minPert, maxPert, nPert)';
PointVec(:,1) = reshape( repmat(PertVec', nPert, 1), nPoints, 1);
PointVec(:,2) = repmat( PertVec, nPert, 1);

for irxn1 = 1:nRxn-1
    
    % Only try to calculate FIC if both FCCs are non-zero
    if abs(FCC(irxn1)) < AbsTol
        continue;
    end
    
    CurrFullRxnName1 = RxnNames{irxn1};
    
    for irxn2 = irxn1+1:nRxn
        
        % Only try to calculate FIC if both FCCs are non-zero
        if abs(FCC(irxn2)) < AbsTol
            continue;
        end
        
        CurrFullRxnName2 = RxnNames{irxn2};
        
        fprintf('%s and %s: \n', CurrFullRxnName1, CurrFullRxnName2);
        
        FluxVec = nan( nPoints, 1);
        
        for ipoint = 1:nPoints
            mutPert = PointVec(ipoint,:)' ;
            [MutModel, FluxDistr] = getMutFlux(ModelObj, {CurrFullRxnName1; CurrFullRxnName2}, mutPert);
            FluxVec(ipoint) = sum( ModelSpecs.rxnOutEnzymeCoeff .* FluxDistr.Flux( ModelSpecs.rxnOutIX) );
        end
        
        X1 = PointVec(:,1)-1;
        X2 = PointVec(:,2)-1;
        
        Y = (FluxVec/WTOutFlux - 1)...
            - FCC(irxn1) * X1 ...
            - FCC(irxn2) * X2 ...
            - 1/2 * FICD(irxn1,irxn1) * X1.^2 ...
            - 1/2 * FICD(irxn2,irxn2) * X2.^2 ;
                
        b = regress(Y , X1 .* X2 );
        
        %%% Diagonose problems by plotting perdicted vs observed
        % clf; hold on; box on;
        %         plot( b(1) * X1 .* X2, Y, 'ok');
        %         plot( [min(Y), max(Y)], [min(Y), max(Y)], 'k--');
        %         xlabel('Predicted');
        %         ylable('Observed');
                
        if abs(b(1)) > AbsTol
            FIC(irxn1, irxn2) = b(1);
            FIC(irxn2, irxn1) = b(1);            
        end
        
        fprintf('\tFIC = %.4f\n', FIC(irxn1,irxn2) );
    end
end