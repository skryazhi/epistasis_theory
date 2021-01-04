function [ISMODULE, plIO ] = isModule( RS )
% Evaluates whether the rxn set specified by RS (which is a 1 by
% nrxn boolean vector) induces a module. ISMODULE is TRUE if it
% does. plIO is the cell array of all simple paths connecting the
% IO metabolites (if it does in fact induce a module)

A = squareform( RS );
nTot = size(A,1);
ISMODULE = false;
plIO = {};

% First check: each IO metabolite must have at least one other node to which
% it is connected. Otherwise this reaction set does not induce a module
if sum(A(1,:)) < 1 || sum(A(2,:)) < 1
    return;
end

% Second check: each I metabolite must have at least two other
% nodes to which it is connected. Otherwise this reaction set does
% not induce a modules
if any( sum(A(3:end,:),2) < 2 )
    return;
end

% Third check. Each metabolite must be part of at least one simple path
% connecting the IO metabolites
plIO = getAllSimplePaths( A , 1, 2);

for iMet = 3:nTot
   TF = cellfun(@(p) ismember(iMet,p), plIO);
   if ~any(TF)
       return;
   end
end


%%% If we are here, the current set of reactions does induce a module
ISMODULE = true;
% plIO = pl{1,2};
end
