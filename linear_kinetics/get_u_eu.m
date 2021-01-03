function [u, eu] = get_u_eu( Xmat, Kmat )
% GET_U_EU theoretical calculation of the effective parameter u (after
% elimination of nodes 3 and 4) and epistasis eu for u in sub-module mu of
% module nu (Figure 3A). Note that metabolites 2 and 5 have switched
% labels here relative to Figure 3A.

% [U, EU] = GET_U(XMAT, KMAT) takes the 4 by 4 matrix of rate constants
% XMAT and a 4 by 4 matrix of equilibrium constants KMAT (XMAT must obey
% Haldane relationships with respect to KMAT) that corresponds to module mu
% depicted in Figure 3A that consists of metabolites 1, 5, 3, 4 (in that
% order in KMAT and XMAT). Metabolites 1 and 5 are I/O, metabolites 3 and 4
% are internal. U is the rate constant for the effective reaction between
% metabolites 1 and 5 after elimination of metabolites 3 and 4, calculated
% analytically through the application of the recursive node elimination.
% EU is the epistasis for U, calculated using the Equation (61).


%% Calculate y through recursive node elimination:
u = get_effective_rate( Xmat, Kmat );


%% Calculate epistasis using the analytically obtained formula:
c1 = 1/Kmat(2,4) * sum( Xmat(2,3:4) );
c2 = sum( Xmat(1,3:4) );
a = c1 * c2;

b = sum( Xmat(3,1:2) ) * Xmat(1,4) * Xmat(4,2) + sum( Xmat(4,1:2) ) * Xmat(1,3) * Xmat(3,2) ;
d1 = Xmat(3,2) * sum( Xmat(4,1:2) ) ;
d2 = Xmat(1,4) * sum( Xmat(3,1:2) );

z = Xmat(3,4);
eu = z * (a * z + b)/( c1 * z + d1) / (c2 * z + d2 );

end
