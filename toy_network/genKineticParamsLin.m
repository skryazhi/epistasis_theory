function [VmMat, KeqMat] = genKineticParamsLin( AdjMat, DistrType, DistrParam )

%genKineticParamsLin generates kinetic parameter matrix for network of
%linear reactions
%
% VmMat = genKineticParamsLin( AdjMat, DistrType, DistrParam ) takes the adjacency
% matrix AdjMat of metabolites and generates the matrix of corresponding
% equilibrium constants KeqMat and kinetic parameters VmMat. AdjMat(i,j) =
% 1 indicates that there is a reaction between metabolites i and j. AdjMat
% is symmetric.
%
% DistrType is the type of the distribution. Current options are 'Delta'
% (point-measure), 'Exponential' (exponential distribution) or 'Gamma'
% (gamma distribution). Default is 'Exponential'.
%
% DistrParam is the parameter or parameters for the distribution from
% which the kinetic parameters and equilibrium constants are drawn. Default
% DistrParam = 1. For the gamma distribution, DistrParam must be a 1 by 2
% vector with DistrParam(1) specifying the shape parameter k and
% DistrParam(2) specifying the scale parameter theta. y=f(x|k,theta)= 1/theta^k/Gamma(k)
% x^{k-1} exp(-x/theta)
%
% KeqMat is the matrix of equilibrium constants, such that KeqMat(i,j) is
% the equilibrium constant with i as the source metabolite and j as the
% sink metabolite, which means that, at equilibrium Si = Sj/Kij. KeqMat =
% 1./KeqMat';
%
% VmMat is the matrix of Vmax values for the reactions. VmMat has non-zero
% values only in those entries where AdjMat is non-zero. VmMat(j,i) =
% VmMat(i,j)/KeqMat(i,j) (Haldane's relationship).

if nargin < 2
    DistrType = 'Exponential';
elseif ~strcmp(DistrType, 'Delta') && ~strcmp(DistrType, 'Exponential') && ~strcmp(DistrType, 'Gamma')
        error('%s distribution is not supported', DistrType);
end

if nargin < 3
    DistrParam = 1;
end


nm = size(AdjMat, 1); % number of metabolites
nrxn = nnz(AdjMat) / 2; % number of rxns

if strcmp(DistrType, 'Delta')
    if DistrParam(1) <= 0
        error('Parameter of the Point-measure must be non-negative');
    end
    Keqvec = DistrParam(1) * ones(nm, 1);
    Vmvec = DistrParam(1) * ones(nrxn, 1);
elseif strcmp(DistrType, 'Exponential')
    if DistrParam(1) <= 0
        error('Parameter of the Exponential distribution must be non-negative');
    end
    Keqvec = random('Exponential', DistrParam(1), [nm, 1]);
    Vmvec = random('Exponential', DistrParam(1), [nrxn, 1]);
elseif strcmp(DistrType, 'Gamma')
    if length(DistrParam) < 2
        error('DistrParam must have 2 parameters for the Gamma distribution');
    elseif any(DistrParam <= 0)
        error('Both parameters of the Gamma distribution must be non-negative');
    end
    Keqvec = random('Gamma', DistrParam(1), DistrParam(2), [nm, 1]);
    Vmvec = random('Gamma', DistrParam(1), DistrParam(2), [nrxn, 1]);
end

KeqMat = repmat( Keqvec, 1, nm) ./ repmat(Keqvec', nm, 1); % Matrix of equilibrium constants
% KeqMat = KeqMat - diag(diag(KeqMat)) + eye( nm ); 
clear Keqvec;

Tvec = squareform(AdjMat);
Tvec( Tvec > 0 ) = Vmvec;

VmMat = triu(squareform(Tvec, 'tomatrix'), 1);
VmMat = VmMat + tril(VmMat'./KeqMat' ,-1);
clear Vmvec Tvec;