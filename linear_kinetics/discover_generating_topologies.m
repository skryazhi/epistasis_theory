% This code implements the algorithm for discovering all strictly serial
% and strictly parallel generating reaction sets in each of the nine
% topological classes.

ClassList = {'bioIO', 'bi0', 'ioioI', 'ioioIO', 'ioio0', 'ioiI', 'ioi0', 'iiI', 'ii0'};

%% Definition of topological classes

% Class M^{b,io,IO}
ClassDef.bioIO.n = 1;
ClassDef.bioIO.a = [1 2];
ClassDef.bioIO.b = [1 3];
ClassDef.bioIO.sym = {};

% Class M^{b,i,0}
ClassDef.bi0.n = 2;
ClassDef.bi0.a = [1 2];
ClassDef.bi0.b = [3 4];
ClassDef.bi0.sym = {[1;2], [3;4]};

% Class M^{io,io,I}
ClassDef.ioioI.n = 1;
ClassDef.ioioI.a = [1 3];
ClassDef.ioioI.b = [2 3];
ClassDef.ioioI.sym = {[1; 2]};

% Class M^{io,io,IO}
ClassDef.ioioIO.n = 2;
ClassDef.ioioIO.a = [1 3];
ClassDef.ioioIO.b = [1 4];
ClassDef.ioioIO.sym = {[3;4]};

% Class M^{io,io,0}
ClassDef.ioio0.n = 2;
ClassDef.ioio0.a = [1 3];
ClassDef.ioio0.b = [2 4];
ClassDef.ioio0.sym = {[1 3; 2 4]};

% Class M^{io,i,I}
ClassDef.ioiI.n = 2;
ClassDef.ioiI.a = [1 3];
ClassDef.ioiI.b = [3 4];
ClassDef.ioiI.sym = {};

% Class M^{io,i,0}
ClassDef.ioi0.n = 3;
ClassDef.ioi0.a = [1 3];
ClassDef.ioi0.b = [4 5];
ClassDef.ioi0.sym = {[4;5]};

% Class M^{i,i,I}
ClassDef.iiI.n = 3;
ClassDef.iiI.a = [3 4];
ClassDef.iiI.b = [3 5];
ClassDef.iiI.sym = {[1;2], [4;5]};

% Class M^{i,i,0}
ClassDef.ii0.n = 4;
ClassDef.ii0.a = [3 4];
ClassDef.ii0.b = [5 6];
ClassDef.ii0.sym = {[1;2], [3;4], [5;6], [3 4; 5 6]};

for iClass = 1:length(ClassList)
    
    fprintf('Class %s ...\n', ClassList{iClass} );
    
    CurrClass = ClassDef.( ClassList{iClass} );
    
    [S.( ClassList{iClass} ) , P.( ClassList{iClass} )]...
        = getSerParGenSets( CurrClass.n , CurrClass.a, CurrClass.b );
    
    if isempty( CurrClass.sym )
        S.( ClassList{iClass} ).rmrxns_nosym = S.( ClassList{iClass} ).rmrxns;
        P.( ClassList{iClass} ).rmrxns_nosym = P.( ClassList{iClass} ).rmrxns;
        continue;
    end
    
    fprintf('removing symmetries ... \n' );
    S.( ClassList{iClass} ).rmrxns_nosym = remSym( S.( ClassList{iClass} ).rmrxns , CurrClass.sym );
    P.( ClassList{iClass} ).rmrxns_nosym = remSym( P.( ClassList{iClass} ).rmrxns , CurrClass.sym );
end

save('AllSerParGenRxnSets.mat', 'ClassDef', 'ClassList', 'P', 'S');











