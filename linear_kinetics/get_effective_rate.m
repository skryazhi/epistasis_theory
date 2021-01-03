function y = get_effective_rate( Xmat, Kmat )
% GET_EFFECTIVE_RATE theoretical calculation of the effective parameter
% rate (after eliminnation of internal metabolites)

% Y = GET_EFFECTIVE_RATE(XMAT, KMAT) takes an arbitrary matrix of rate
% constants XMAT and the corresponding matrix of equilibrium constants KMAT
% (XMAT must obey Haldane relationships with respect to KMAT). Metabolites
% 1 and 2 are considered I/O, all other metabolites are internal. Carry out
% the recursive metabolite elimination procedure for all internal
% metabolites/

y = eliminate_node( Xmat );

return;

    function u = eliminate_node( x )
        n = size(x,1);
        if n == 2
            u = x(1,2);
            return;
        else
            x1 = zeros(n-1,n-1);
            D = sum( x(n,:) );
            for i1 = 1:(n-2)
                for i2 = (i1+1):(n-1)
                    x1(i1,i2) = x(i1,i2) + x(i1,n) * x(n,i2) / D;
                    x1(i2,i1) = x1(i1,i2) / Kmat(i1,i2);
                end
            end
            u = eliminate_node( x1 );
        end
    end

end
