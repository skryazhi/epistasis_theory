function [flux, delta] = get_epistasis(OrigModelObj, rxnOutId, mutA, mutB, ExtGluVec, Tol)

warning('off', 'all');

if nargin < 6
    Tol = 1e-5;
end

ixA = mutA.rxnId;
ixB = mutB.rxnId;

nA = length(mutA.relRmax);
nB = length(mutB.relRmax);
nGlu = length(ExtGluVec);

flux.WT = nan( 1, 1, nGlu);
flux.A = nan( nA, 1, nGlu);
flux.B = nan( 1, nB, nGlu);
flux.AB = nan( nA, nB, nGlu);

delta.A = nan( nA, 1, nGlu);
delta.B = nan( 1, nB, nGlu);
delta.AB = nan( nA, nB, nGlu);
delta.eps = nan( nA, nB, nGlu);

for iglu = 1:nGlu
    fprintf('=======\n');
    fprintf('C_glu = %.2f:\n', ExtGluVec(iglu));
    set(OrigModelObj.Species(1), 'InitialAmount', ExtGluVec(iglu) );
    [success, variant_out, mWT] = sbiosteadystate(OrigModelObj);
    
    rxnWT = get_fluxes(mWT);
    flux.WT(1,1,iglu) = rxnWT.Flux(rxnOutId);
    
    
    for iA = 1:nA
        deltaXa = mutA.relRmax(iA);
        
        % Mut A
        MutModelObj = copyobj( OrigModelObj );
        x = get( MutModelObj.Reactions(ixA).KineticLaw.Parameters(1), 'Value');
        set(MutModelObj.Reactions(ixA).KineticLaw.Parameters(1), 'Value', x*deltaXa );
        [success, variant_out, m2] = sbiosteadystate(MutModelObj);
        rxn = get_fluxes(m2);
        flux.A(iA,1,iglu) = rxn.Flux(rxnOutId);
        delta.A(iA,1,iglu) = (flux.A(iA,1,iglu) - flux.WT(1,1,iglu))/flux.WT(1,1,iglu);
        
        if abs(flux.A(iA,1,iglu)-flux.WT(1,1,iglu)) < Tol
            delta.A(iA,1,iglu) = 0;
        end
        
        fprintf('Mutation A affects %s, delta^A f = %.3g\n',...
            get( MutModelObj.Reactions(ixA).KineticLaw.Parameters(1), 'Name'), delta.A(iA,1,iglu) );
    end
    
    for iB = 1:nB
        deltaXb = mutB.relRmax(iB);
        
        % Mut B
        MutModelObj = copyobj( OrigModelObj );
        x = get( MutModelObj.Reactions(ixB).KineticLaw.Parameters(1), 'Value');
        set(MutModelObj.Reactions(ixB).KineticLaw.Parameters(1), 'Value', x*deltaXb );
        [success, variant_out, m2] = sbiosteadystate(MutModelObj);
        rxn = get_fluxes(m2);
        flux.B(1,iB,iglu) = rxn.Flux(rxnOutId);
        delta.B(1,iB,iglu) = (flux.B(1,iB,iglu) - flux.WT(1,1,iglu))/flux.WT(1,1,iglu);
        if abs(flux.B(1,iB,iglu)-flux.WT(1,1,iglu)) < Tol
            delta.B(1,iB,iglu) = 0;
        end
        fprintf('Mutation B affects %s, delta^B f = %.3g\n',...
            get( MutModelObj.Reactions(ixB).KineticLaw.Parameters(1), 'Name'), delta.B(1,iB,iglu) );
    end
            
    for iA = 1:nA
        for iB = 1:nB
            deltaXa = mutA.relRmax(iA);
            deltaXb = mutA.relRmax(iB);
            
            % Mut AB
            MutModelObj = copyobj( OrigModelObj );
            
            x = get( MutModelObj.Reactions(ixA).KineticLaw.Parameters(1), 'Value');
            set(MutModelObj.Reactions(ixA).KineticLaw.Parameters(1), 'Value', x*deltaXa );
            
            x = get( MutModelObj.Reactions(ixB).KineticLaw.Parameters(1), 'Value');
            set(MutModelObj.Reactions(ixB).KineticLaw.Parameters(1), 'Value', x*deltaXb );
            
            [success, variant_out, m2] = sbiosteadystate(MutModelObj);
            rxn = get_fluxes(m2);
            flux.AB(iA,iB,iglu) = rxn.Flux(rxnOutId);
            delta.AB(iA,iB,iglu) = (flux.AB(iA,iB,iglu) - flux.WT(1,1,iglu))/flux.WT(1,1,iglu);
            if abs(flux.AB(iA,iB,iglu)-flux.WT(1,1,iglu)) < Tol
                delta.AB(iA,iB,iglu) = 0;
            end
            
            delta.eps(iA,iB,iglu) =  (delta.AB(iA,iB,iglu) - delta.A(iA,1,iglu) - delta.B(1,iB,iglu) )/ delta.A(iA,1,iglu) / delta.B(1,iB,iglu);
            
            if abs(delta.AB(iA,iB,iglu) - delta.A(iA,1,iglu) - delta.B(1,iB,iglu)) < Tol
                delta.eps(iA,iB,iglu) = 0;
            end
            
            fprintf('delta^{AB} f = %.3g, eps f = %.4f\n', delta.AB(iA,iB,iglu), delta.eps(iA,iB,iglu) );
        end
    end
end

