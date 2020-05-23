function [S, P] = getSerParGenSets( nint, a, b )
% getSerParGenSets produces the generating reaction sets for all serial and
% parallel modules within a particular topological class
%
% [S, P] = getSerParGenSets( nint , a, b ) takes as input the number of
% internal vertices (metabolites) nint in the module and the two marked
% reactions a and b. Each a and b are 1 by 2 arrays of vertices
% (metabolites).
%
% S and P are structures that correspond to serial and parallel generating
% sets, respectively. S.RS is a binary M by NRXN array representing the
% generating set for serial modules. Here, NRS is the number of rxn sets in
% the generating set and NRXN is the total number of possible reactions in
% the module. S.RS(iRS, jRXN) is TRUE if in the reaction set iRS reaction
% jRXN is present. S.rmrxns is the cell array of rxns that have been
% removed from the full set. Analogously for P.RS and P.rxns

nTot = nint + 2;
nRxn = nTot * (nTot - 1)/2; % total number of rxns

MarkedRxns = [a ; b];

B = zeros( nTot );
B(MarkedRxns(1,1), MarkedRxns(1,2)) = 1;
B(MarkedRxns(1,2), MarkedRxns(1,1)) = 1;
B(MarkedRxns(2,1), MarkedRxns(2,2)) = 1;
B(MarkedRxns(2,2), MarkedRxns(2,1)) = 1;
MarkedRxnsIx = find( squareform( B ) ); % ids of the marked rxns (these cannot be removed)
clear B;

%%% Adjacency matrices (in vector form); use squareform to convert to
%%% matrix form
S.RS = nan( 0, nRxn );
P.RS = nan( 0, nRxn );

RSFULL = true(1, nRxn); % The full Reaction Set for the current topological class

plIO = getAllSimplePaths( squareform( RSFULL ) , 1, 2 );
[TFS, TFP] = isSP( plIO , MarkedRxns );

if TFS && ~TFP % SERIAL
    S.RS(end+1,:) = RSFULL;
elseif ~TFS && TFP % PARALLEL
    P.RS(end+1,:) = RSFULL;
else % Neither serial nor parallel (i.e., both S and P conditions are true)
    RSList = rmRxn( RSFULL, MarkedRxnsIx );
    ReduceRxnSet( RSList );
end

S.rmrxns = cell( size( S.RS, 1), 1 );
P.rmrxns = cell( size( P.RS, 1), 1 );

for iRS = 1:size( S.RS, 1)
    X = squareform( S.RS(iRS,:) );
    ix = find( X == 0 );
    [ii, jj] = ind2sub([nTot, nTot], ix);
    TF = ii < jj;
    S.rmrxns{iRS} = [ii(TF), jj(TF)];
end
for iRS = 1:size( P.RS, 1)
    X = squareform( P.RS(iRS,:) );
    ix = find( X == 0 );
    [ii, jj] = ind2sub([nTot, nTot], ix);
    TF = ii < jj;
    P.rmrxns{iRS} = [ii(TF), jj(TF)];
end
clear i1 X ix;

    
    function [] = ReduceRxnSet( RSList )
        % This function takes a list of reaction sets RS (which are binary
        % adjacency matrices in the vector form), evaluates whether
        % these reaction sets are serial or parallel (w.r.t. the marked
        % rxns a and b) and populates S.A and P.A
        
        newRSList = [];
        
        for iRS = 1:size(RSList,1)
            
            % Current RS
            CurrRS = RSList( iRS, :);
            
            %%% Does current RS induce a module?
            [ISMODULE, plIO ] = isModule( CurrRS );
            
            % If it does not, then move on to the next one
            if ~ISMODULE
                continue;
            end
            
            IFPP = false;
            IFPS = false;
            
            IFPS1 = false;
            IFPP1 = false;
            
            % Is current RS a generating RS? To answer that, test whether
            % all "parent" RS are neither serial nor parallel. Parent RSs
            % of CurrRS is CurrRS + any rxn
            for ixRmRxn = find( ~CurrRS )
                ParentRS = CurrRS;
                ParentRS(ixRmRxn) = true;
                                               
                % These are all rxns present in the parent for checking if
                % it is serial
                TFS = repmat(ParentRS, size(S.RS,1), 1);

                % These are all rxns present in the parent for checking if
                % it is parallel
                TFP = repmat(ParentRS, size(P.RS,1), 1);
                               
                % If parent RS can be created from any S.RS by removing
                % rxns, it is serial. If parent RS can be created from any
                % P.RS by removing rxns, it is parallel
                if any( all( S.RS | ~xor(S.RS, TFS) , 2 ), 1)
                    IFPS = true;
                end
                if any( all( P.RS | ~xor(P.RS, TFP) , 2 ), 1)
                    IFPP = true;
                end
                
                % If one parent is either parallel or serial, then stop
                % looking at other parents
                if IFPS || IFPP
                    break;
                end
            end
            
            % If any parent is either serial or parallel, then the current
            % RS inherits this property. Move on to the next RS
            if IFPS || IFPP
                continue;
            end
            
            % Otherwise, check if the current RS is parallel or serial
            [TFS, TFP] = isSP( plIO , MarkedRxns );
            
            if TFS && ~TFP % SERIAL
                S.RS(end+1,:) = CurrRS;
            elseif ~TFS && TFP % PARALLEL
                P.RS(end+1,:) = CurrRS;
            else % Neither serial nor parallel (i.e., both S and P conditions are true)
                % Then, generate all of its children and add them to the
                % new RS list
                newRSList = [newRSList; rmRxn( CurrRS , MarkedRxnsIx )];
            end
        end
        
        if isempty( newRSList )
            return;
        else
            newRSList = flipud( unique(newRSList, 'rows') );
            ReduceRxnSet( newRSList );
        end
    end
end



function RSList = rmRxn( RS , MarkedRxnsIx )
% remove each rxn from the rxn set RS and generate the list of reduced
% reaction sets RSList

ixRxns = find( RS ); % reactions that are in the current RS

% Can't have less than 2 rxns in the module:
if length(ixRxns) == 2
    RSList = [];
    return;
end

RSList = repmat(RS, length(ixRxns)-2, 1);

k = 1;
for ix = 1:length(ixRxns)
    ixRxn = ixRxns(ix);
    % Cannot remove marked reactions
    if ixRxn == MarkedRxnsIx(1) || ixRxn == MarkedRxnsIx(2)
        continue;
    end
    
    if k > size(RSList, 1)
        error('Oops! This should not happen!');
    end
    
    RSList(k, ixRxn ) = false;
    k = k+1;
end

RSList = flipud( RSList );

end




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



function [TFS, TFP] = isSP( plIO , MarkedRxns  )
% Evaluates whether marked rxn are serial and/or parallel. plIO is the cell
% array of simple paths connecting the IO metabolites. MarkedRxns (a 2 by 2
% matrix) is the list of two marked rxns. TFS is TRUE if there is at
% least one path where the marked reactions are serial (condition S is
% satisfied), TFP is TRUE if there are paths where the marked rxns are
% parallel (condition P is satisfied)

nPathsIO = length( plIO ); % number of paths between the IO metabolites
ISINPATH = false( nPathsIO, 2); % whether the two marked rxns are in each path
TFS = false;
TFP = false;

% Determine whether each of the marked rxns is in each path
for ip = 1:nPathsIO
    ISINPATH(ip, :) = isinpath( plIO{ip} , MarkedRxns)';
end

if any( all( ISINPATH, 2 ) )
    TFS = true;
end

if any( ISINPATH(:,1) & ~ISINPATH(:,2) ) && any( ~ISINPATH(:,1) & ISINPATH(:,2) )
    TFP = true;
end

end



function TF = isinpath( p , l )
% ISINPATH checks if the path passes through marked reactions
%
% TF = ISINPATH( P , L ) checks if path P passes through marked reactions
% specified in L. P is an 1 by m array which is an ordered list of
% metabolites through which the path passes. TF is an 2 by 1 boolean array
% containing TRUE if the corresponding rxn is in the pathway and FALSE
% otherwise.

p_rxns1 = [p(1:end-1) ; p(2:end)]';
p_rxns2 = [p(2:end) ; p(1:end-1)]';

TF = ismember( l , p_rxns1, 'rows' ) | ismember( l , p_rxns2, 'rows' );
end



function PathList = getAllSimplePaths(A, iB, iE)
% GETALLSIMPLEPATHS computes all simple paths between two nodes in a graph
%
% PathList = GETALLSIMPLEPATHS( A, IB, IE ) takes as input the adjacency
% matrix A and computes all simple paths between nodes IB and IE. PATHLIST
% is a cell array where each cell is a list of nodes (metabolites) starting
% with IB and ending in IE

PathList = extend_path( iB );

% ISDEADEND = cellfun(@(p) p(end) == -1, pl, 'UniformOutput', true);
% PathList = pl( ~ISDEADEND );

p_len = cellfun(@length, PathList);
[tmp, ix] = sort(p_len, 'ascend');

PathList = PathList( ix );



    function pl = extend_path(p)
        % Extends path p in all possible directions. Returns the cell array
        % of paths that share p and then diverge
        
        curr_node = p(end);
        next_potential = find( A( curr_node, :) ); % potential nodes which path p can visit from the current node

        % Check which nodes have not been visited yet by path p
        ISNEW = ~ismember( next_potential, p );    
        
        % If there are no new nodes to visit, then stop
        if all(~ISNEW)
            pl = {};
            return;
        end
        
        next_possible = next_potential( ISNEW ); % List of nodes that have not been visited yet by this path
            
        % If the terminal node is not 1 or 2, then the path cannot pass
        % through 1 or 2
        if iE ~= 1 && iE ~= 2
            TF = next_possible ~= 1 & next_possible ~= 2;
            next_possible = next_possible( TF );
        end
        
        n = length(next_possible); % number of possible path extensions
        
        % Add the next node
        pl1 = [repmat(p, n, 1), next_possible'];
        % pl1 is the array of paths extended from p by one node
                
        pll = cell(n, 1);
        % pll is a cell array, such that each cell is a cell array of paths
        % that share the first part up to (and including) the currently
        % added node and subsequently diverge               
        
        % Extending each of these paths in pl1 further
        for i1 = 1:n
            if pl1(i1,end) == iE %|| pl1(i1,end) == 1 || pl1(i1,end) == 2
                % if the path reached node iE
                pll{i1} = { pl1(i1,: ) };
            else % otherwise, continue adding nodes
                pll{i1} = extend_path( pl1(i1,:) );
                
                % At this point, only paths that terminate at iE should be retained.
                % Paths that terminate at other nodes, are removed
            end
        end        
        
        np = cellfun(@length, pll); % vector or the number of paths in each extension
        
        pl = cell( sum(np) , 1); % cell array of all paths that share p
        
        % "Flattening" pll into pl
        k = 1;
        for i1 = 1:n
            pl(k:k+np(i1)-1) = pll{i1};
            k = k + np(i1);
        end        
    end
end