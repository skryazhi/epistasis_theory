function [u, eu] = get_u( Xmat, Kmat )
% GET_U theoretical calculation of the effective parameter u (after
% eliminnation of nodes 3 and 4) and epistasis for u in a simple module
% depicted in Figure 5

% [U, EU] = GET_U(XMAT, KMAT) takes the 4 by 4 matrix of rate constants
% XMAT and a 4 by 4 matrix of equilibrium constants KMAT (XMAT must obey
% Haldane relationships with respect to KMAT) that corresponds to the
% module depicted in Figure 2A that consists of metabolites 1, 5, 3, 4 (in
% that order in KMAT and XMAT). Metabolites 1 and 5 are I/O, metabolites 3
% and 4 are internal. U is the rate constant for the effective reaction
% between metabolites 1 and 5 after elimination of metabolites 3 and 4,
% calculated analytically through the application of the recursive node
% elimination. EU is the epistasis for U, calculated using the formula
% (A4.6) in Appendix 4.


%% Calculate y through recursive node elimination:
u = get_effective_rate( Xmat, Kmat );


%% Calculate epistasis using the analytically obtained formula:
a = 1/Kmat(2,4) * sum( Xmat(1,3:4) ) * sum( Xmat(2,3:4) );
b = sum( Xmat(3,1:2) ) * Xmat(1,4) * Xmat(4,2) + sum( Xmat(4,1:2) ) * Xmat(1,3) * Xmat(3,2) ;
c1 = 1/Kmat(2,4) * sum( Xmat(2,3:4) );
d1 = Xmat(3,2) * sum( Xmat(4,1:2) ) ;
c2 = sum( Xmat(1,3:4) );
d2 = Xmat(1,4) * sum( Xmat(3,1:2) );

z = Xmat(3,4);
eu = z * (a * z + b)/( c1 * z + d1) / (c2 * z + d2 );

end
