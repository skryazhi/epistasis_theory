function mutFlux = getMutFlux(OrigModelObj, rxnOutId, mutList)

warning('off', 'all');

if nargin < 3 || isempty(mutList)
    [success, variant_out, m2] = sbiosteadystate(OrigModelObj);
    rxn = get_fluxes(m2);
    mutFlux = rxn.Flux(rxnOutId);
    return;
end



nMut = length( mutList.Id );
mutFlux = cell( nMut, 1);

for imut = 1:nMut
    
    fprintf('Mut #%d out of %d (%.1f%%)\n', imut, nMut, imut/nMut*100 );
    
    % which reactions/parameters are affected in this mutant?
    rxnIdVec = mutList.Id{imut};
    nRxn = size( mutList.Pert{imut}, 1 );
    nPert = size( mutList.Pert{imut}, 2 );
    
    mutFlux{imut} = nan(1, nPert);

    for iPert = 1:nPert

        MutModelObj = copyobj( OrigModelObj );

        for iRxn = 1:nRxn

            rxnId = rxnIdVec(iRxn);
                        
            x = get( MutModelObj.Reactions(rxnId).KineticLaw.Parameters(1), 'Value') * mutList.Pert{imut}(iRxn, iPert);
            
            set(MutModelObj.Reactions(rxnId).KineticLaw.Parameters(1), 'Value', x  );
        end
        
        [success, variant_out, m2] = sbiosteadystate(MutModelObj);
        rxn = get_fluxes(m2);
        mutFlux{imut}(iPert) = rxn.Flux(rxnOutId);
    end
end

warning('on');

