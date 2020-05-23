function new_rss = remSym( rss , SymList )
%
% REMSYM takes a set of reaction sets (RS) and removes those RSs
% that are symmetric with respect to permutations of node labels
%
% NEWSS = REMSYM( RSS, SYMLIST ) takes as in put RSS, the set of
% reaction sets, which is a cell array of reaction matrices such that each
% cell contains an N by 2 array of rxns (rows = rxns, columns = metabolite
% nodes), and the cell array of symmetries SYMLIST, where each cell is a 2 by k
% list of metabolite indices that can be swapped. Specifically, the first
% row of SYMLIST{i} specifies a list of metabolite labels that can be swapped
% (as a group) with the corresponding list of metabolite labels in the
% second row of SYMLIST{i}.
%
% EXAMPLE 1. SYMLIST{1} = [1;2] means that metabolite
% labels 1 and 2 are exchangable;
%
% EXAMPLE 2. SYMLIST{1} = [3 4; 5 6] means that the pair of metabolite labels
% 3, 4 can be swapped with labels 5 and 6, respectively.
%
% RSS is the split into equivalence classes, such that all RSs from the same
% equivalence class can be transformed into each other by relabelling the
% metabolite nodes according to SYMLIST. One member of each equivalence class
% is retained to generate NEW_RSS (also a cell array). The rule is to retain
% that RS which contains RS with the lowest lexicographic rank


nrs = length(rss); % number of RSs

if nrs < 2
    new_rss = rss;
    return;
end

tmp_ss = rss;

% Sort rows of each RS in a lexicographic order
for is = 1:nrs
    tmp_ss{is} = sortrows( rss{is} );
end
rss = tmp_ss;
clear tmp_ss;

% Start by assuming that each RS belongs to its own equivalence class
EquivClass = (1:nrs)';
for isym = 1:length(SymList)
    EquivClass = ApplySym( rss, rss, SymList(isym:end), EquivClass);
end


% Populate new_rss with a single representative from each equivalence class
unqEquivClasses = unique( EquivClass );
new_rss = cell( length(unqEquivClasses), 1 );

for iec = 1:length(unqEquivClasses)
    ec = unqEquivClasses(iec);
    new_rss{iec} = getEquivClassRep( rss(EquivClass == ec) );
end
end


function newEquivClass = ApplySym( rssOrig, rssCurr, SymList, EquivClass )
% Apply all symmetries recursively and update equivalence class memberships
% accordingly

if isempty( SymList )
    newEquivClass = EquivClass;
    return;
end

rssNew = premute_labels( rssCurr , SymList{1} );
newEquivClass = UpdateEquivClass( rssOrig, rssNew, EquivClass );
newEquivClass = ApplySym( rssOrig, rssNew, SymList(2:end), newEquivClass );

end




function newEquivClass = UpdateEquivClass( rssOrig, rssNew, EquivClass )
% Compare original set of RS rssOrig with the new set of RS rssNew and
% update the equivalence classes

% Obtain the index of each new RS in rssOrig
ixNewInOrig = zeros( size( rssNew ));
for ixNew = 1:size(rssNew,1)
    for ixOrig = 1:size(rssOrig,1)
        if isequal( rssNew{ixNew} , rssOrig{ixOrig} )
            ixNewInOrig( ixNew ) = ixOrig;
        end
    end
end

% Now compare the orignal equivalence class of each RS with the equivalence
% class of its new (permuted) version. Retain the lower of the two
newEquivClass = EquivClass( ixNewInOrig );
newEquivClass = min( newEquivClass, EquivClass );

end






function rssNew = premute_labels( rssCurr , symCurr )
% Permute labels of all metabolites in the RS set according to the current
% symmetry

rssNew = cell( size(rssCurr) );

for irs = 1:length(rssNew)
    rssNew{irs} = rssCurr{irs};
    
    for iMetIx = 1:size(symCurr,2)
        TF = rssCurr{irs} == symCurr(1,iMetIx);
        rssNew{irs}(TF) = symCurr(2,iMetIx);
        
        TF = rssCurr{irs} == symCurr(2,iMetIx);
        rssNew{irs}(TF) = symCurr(1,iMetIx);
    end
    rssNew{irs} = sort( rssNew{irs}, 2);
    rssNew{irs} = sortrows( rssNew{irs} );
end

end






function RSsRep = getEquivClassRep( EquivRSs )
% Select a single representative from the equivalence class

% Convert to 3D matrix
M = nan( size(EquivRSs{1},1), size(EquivRSs{1},2), length(EquivRSs) );
for irs = 1:length(EquivRSs)    
    M(:,:,irs) = EquivRSs{irs};
end

% Sorting RSs in descending lexicographic order
for irxn = size(M,1):-1:1
    tmp = permute( M(irxn,:,:), [3 2 1] );
    [tmp, ix] = sortrows(tmp);
    M = M(:,:,flipud(ix));
    
    CurrRxnFirst = permute( M(irxn,:,1), [3 2 1] );
    CurrRxnAll = permute(M(irxn,:,:), [3 2 1]);
    
    % Find all RSs that share the same first irs rxns as the top RS
    TF = all( CurrRxnAll == repmat( CurrRxnFirst, size(CurrRxnAll,1), 1), 2);
    
    % If the top RS is unique along the first irs rxns, then stop
    if nnz(TF) == 1
        break;
    end
    
    % Otherwise, sort those RSs that share the share the same first irs
    % rxns as the top RS along rxn irs+1
    M = M(:,:,TF);
end

RSsRep = M(:,:,1);

end
