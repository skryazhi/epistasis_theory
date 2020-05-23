function [y, ey] = get_u( Xmat, Kmat )
% GET_U theoretical calculation of the effective parameter u (after
% eliminnation of nodes 3 and 4) and epistasis for u in a simple module
% depicted in Figure 5

% [Y, EY] = GET_U(XMAT, KMAT) takes the 4 by 4 matrix of rate constants
% XMAT and a 4 by 4 matrix of equilibrium constants KMAT (XMAT must obey
% Haldane relationships with respect to KMAT) that corresponds to the
% module depicted in Figure 2A that consists of metabolites 1, 5, 3, 4 (in
% that order in KMAT and XMAT). Metabolites 1 and 5 are I/O, metabolites 3
% and 4 are internal. Y is the rate constant Y for the effective reaction
% between metabolites 1 and 5 after elimination of metabolites 3 and 4,
% calculated analytically through the application of the recursive node
% elimination. EY is the epistasis for Y, calculated using the formula
% (XXXXXXXX) in the Supplementary Information.

%% Calculate y through recursive node elimination:
y = eliminate_node( Xmat );


%% Calculate epistasis using the analytically obtained formula:
a = 1/Kmat(2,4) * sum( Xmat(1,3:4) ) * sum( Xmat(2,3:4) );
b = sum( Xmat(3,1:2) ) * Xmat(1,4) * Xmat(4,2) + sum( Xmat(4,1:2) ) * Xmat(1,3) * Xmat(3,2) ;
c1 = 1/Kmat(2,4) * sum( Xmat(2,3:4) );
d1 = Xmat(3,2) * sum( Xmat(4,1:2) ) ;
c2 = sum( Xmat(1,3:4) );
d2 = Xmat(1,4) * sum( Xmat(3,1:2) );

z = Xmat(3,4);
ey = z * (a * z + b)/( c1 * z + d1) / (c2 * z + d2 );

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
