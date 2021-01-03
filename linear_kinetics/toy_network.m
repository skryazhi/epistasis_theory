%%% This code makes the plot of epistasis in a toy metabolic network in
%%% Figure 3B. Note that in this code, metabolites the labels of
%%% metabolites 2 and 5 are switched relative to Figure 3B.

%%% Update 21 Dec 2020 (SK). Minor modifications to the code and
%%% annotations: using get_u_eu function, which itself now relies on
%%% get_effective_rate function

%%% Update 22 Jan 2020 (SK). This code calculates the epistasis coefficient
%%% for the simple 5-node module depicted in Figure 2A.

%% Parameter definitions
S1 = 100;   % Concentration of metabolite 1
S2 = 0;     % Concentration of metabolite 2
% AbsTol = 1e-8;  % Tolerance for numerical calculations

IFCALCULATE = true; % carry out a new calculation or plot only?
IFGENPARAMS = false; % generate new parameters or use the ones shown in the paper?

%%% Adjacency matrix:
AdjMat = [ ...
    0 1 1 1 0;  % 1 is connected to 2, 3, 4
    1 0 0 0 1;  % 2 is connected to 1, 5
    1 0 0 1 1;  % 3 is connected to 1, 4, 5
    1 0 1 0 1;  % 4 is connected to 1, 3, 5
    0 1 1 1 0;  % 5 is connected to 2, 3, 4
    ];

if IFGENPARAMS
    DistrType = 'Exponential';
    DistrParam  = 1;
    
    %%% Draw Vmax and Keq from an exponential distribution
    [VmMat0, KeqMat] = genKineticParamsLin( AdjMat, DistrType, DistrParam );
    clear DistrType DistrParam;
else
    %%% Parmeters used in the paper (where z has a lot of control over flux):
    VmMat0 = [0,0.378294875654394,0.514262129487158,0.236527792446535,0;1.81015374316592,0,0,0,1.00061594228389;42.2317207833885,0,0,1.57974387902022,2.44568731284106;7.95705388943905,0,0.647146442888450,0,0.258856332233772;0,6.98219008225122,0.994388174497328,0.256920007806055,0];
    KeqMat = [1,0.208984942346811,0.0121771530960074,0.0297255486431310,0.0299495806539366;4.78503373865327,1,0.0582680884051414,0.142237753157361,0.143309753887604;82.1210008707108,17.1620526324279,1,2.44109180600492,2.45948953895935;33.6410948038492,7.03048225806603,0.409652761743770,1,1.00753668211461;33.3894491397014,6.97789210345227,0.406588433965495,0.992519694569541,1];
end


%% 01. Calculate
if IFCALCULATE
    
    zvec = (0:0.5:10)'; % vector of z values
    nz = length(zvec);  % total number of z values used
    yvec = nan(nz,1);   % effective rate constant y for the whole module
    epsyvec = nan(nz,1); % epistasis for y
    
    %%% Auxiliary variables:
    uvec = nan(nz,1);  % effective rate constant after coarse-graining nodes 3 and 4
    epsuvec = nan(nz,1);    % epistasis for u
        
    for iz = 1:nz
        z = zvec(iz);
        VmMat = VmMat0;
        VmMat(3,4) = z;
        VmMat(4,3) = z/KeqMat(3,4);
                
        % get u and epsu  through analytical metabolite elimination:
        [uvec(iz), epsuvec(iz)] = get_u_eu( VmMat([1 5 3 4], [1 5 3 4]) , KeqMat([1 5 3 4], [1 5 3 4]) );
        
        % get y by eliminating metabolite 5
        yvec(iz) = VmMat(1,2) + uvec(iz) * VmMat(5,2) / (uvec(iz)/KeqMat(1,5) + VmMat(5,2) );
        
        C = uvec(iz)/yvec(iz) * ( VmMat(5,2) / (uvec(iz)/KeqMat(1,5) + VmMat(5,2) ))^2;
        H = -2 * uvec(iz)^2/yvec(iz)/KeqMat(1,5) * VmMat(5,2)^2 / ((uvec(iz)/KeqMat(1,5) + VmMat(5,2) ))^3;
        
        epsyvec(iz) = 1/C * (epsuvec(iz) + H/2/C );
    end
    
    clear iz z VmMat C H;
end

%% 02. Plot

%%% First, plot epistasis versus z

clf;

%%% Specify the following dimensions:
fdim.spwa = 5.5; % subplotwidth in cm
fdim.spha = 3.5; % subplotheight in cm

fdim.nx = 1; % number of panels along the horizontal dimension
fdim.ny = 1; % number of panels along the vertical dimension

fdim.xma = [0.75 0.2]; % left right horizontal margin in cm
fdim.yma = [0.75 0.2]; % bottom top vertical margin cm

fdim.dxa = 0.3; % horizontal distance between panels in cm
fdim.dya = 0.4; % vertical distance between panels in cm

fdim.tickfs = 8; % tick font size
fdim.labelfs = 10; % axis label font size

%%% These will be computed automatically:
fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
fdim.fh = fdim.spha * fdim.ny + fdim.dya * (fdim.ny - 1) + sum(fdim.yma);

fdim.spwr = fdim.spwa / fdim.fw;
fdim.sphr = fdim.spha / fdim.fh;
fdim.xmr = fdim.xma / fdim.fw;
fdim.ymr = fdim.yma / fdim.fh;
fdim.dxr = fdim.dxa / fdim.fw;
fdim.dyr = fdim.dya / fdim.fh;

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 fdim.fw fdim.fh]);

cc = [
    0, 114, 178;    % blue
    213, 94, 0;     % vermillion
    86, 180, 233;   % sky blue
    230 159, 0;     % orange
    204, 121, 167;   % raddish purple
    0, 158, 115;    % bluish green
    240, 228, 66   % yellow
    ]./256;

fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );

fdim;

subplot('Position', [fdim.spxvec(1) fdim.spyvec(1) fdim.spwr fdim.sphr]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');

plot(zvec, epsyvec, 'k-', 'LineWidth', 3);

set(gca, 'XLim', [0 8], 'Xtick', 0:2:10, 'XGrid', 'off');
set(gca, 'YLim', [-3 3], 'Ytick', [-3 0 1 3], 'YGrid', 'on');

xlabel('z', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');%, 'Position', [55, -0.4]);
ylabel('\epsilon y', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');%, 'Position', [-8 3.1]);

clear fdim cc;



