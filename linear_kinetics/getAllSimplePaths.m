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