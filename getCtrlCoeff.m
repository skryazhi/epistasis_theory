function [CtrlCoeffs, RxnNames] = getCtrlCoeff(ModelObj, ModelSpecs, rxnOutIdVec, OrigOutFlux)
%
% getCtrlCoeff calculates the control coefficients in the model
%
% CtrlCoeff = getCtrlCoeff(ModelObj) takes the simbiology model ModelObj
% and calculates the control coefficients of all reactions

RxnNames = arrayfun(@(x) get(x, 'Name'), ModelObj.Reactions, 'UniformOutput', false);

AbsTolCC = 1e-4;

nRxn = length(RxnNames);
CtrlCoeffs = zeros(nRxn,1);

npoints = 10;
minPert = 0.9;
maxPert = 1.1;


PertVec = linspace(minPert, maxPert, npoints)';
FluxVec = nan(npoints, 1);


for irxn = 1:nRxn
    CurrRxnName = RxnNames{irxn};
    
    for iiter = 1:npoints
        mutList = cell(1,6);
        
        CurrPert = PertVec(iiter);
        mutList(1:2) = [{{CurrRxnName}}, {CurrPert}];
        [ModelL, FluxDistrL] = getMutFlux(ModelObj, mutList(1:2));
        mutList(3:4) = [ModelL, FluxDistrL];
        % Col 1 = enzyme short name list
        % Col 2 = perturbation list
        % Col 3 = perturbed model
        % Col 4 = perturbed steady state fluxes
        % Col 5 = perturbed steady state output flux
        % Col 6 = relative perturbation of output flux
        
        mutList{5} = sum( ModelSpecs.rxnOutEnzymeCoeff .* mutList{1,4}.Flux(rxnOutIdVec) );
        mutList{6} = (mutList{5} - OrigOutFlux)/OrigOutFlux;
        FluxVec(iiter) = mutList{5};
        
        % CCvec(iiter) = mutList{6} / (CurrPert-1);
    end
    
    clf; hold on; box on;
    b = regress(FluxVec/OrigOutFlux-1, [(PertVec-1) , (PertVec-1).^2] );

    plot( PertVec-1, FluxVec/OrigOutFlux-1, 'ok');
    plot( PertVec-1, b(1)*(PertVec-1) + b(2)*(PertVec-1).^2, 'k-', 'LineWidth', 2);
    
    if abs(b(1)) > AbsTolCC
        CtrlCoeffs(irxn) = b(1);
    end
    fprintf('FCC of %s = %.2f%%\n', CurrRxnName, CtrlCoeffs(irxn)*100);
end