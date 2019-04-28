function h = plotFluxDistr(rxn)

% close all;

h = figure;

% Margins
h = 0.95;   % height of the plot
w = 0.95;   % width of the plot
mx = (1-w)/2; % x margin;
my = (1-h)/2; % y margin;

% Text Box
TB.w = 0.075; % width
TB.h = 0.10;  % height
TB.nx = 5;    % number of boxes along x axis
TB.ny = 4;    % number of boxes along y axis
TB.bx  = 0.01; % x buffer
TB.by  = TB.bx / w * h ; % y buffer

% Arrow
A.lx = (w - (TB.w*TB.nx + 2*TB.bx*(TB.nx-1)) ) / (TB.nx - 1);
% A.ly = A.lx / w * h;
A.ly = (h - (TB.h*TB.ny + 2*TB.by*(TB.ny-1)) ) / (TB.ny - 1);
A.w100 = 5;

% lm = 0.05; % Left
c = 0.5;   % center

if TB.nx == 1
    xpos = mx;
else
    xpos = mx + (TB.w + TB.bx + A.lx + TB.bx) * (0:TB.nx-1);
end

if TB.ny == 1
    ypos = my;
else
    ypos = my + (TB.h + TB.by + A.ly + TB.by) * (0:TB.ny-1);
end

tb(1,3,'Glu');
tb(2,3,'g6p');
tb(3,3,'f6p');
tb(4,3,'fdp');
tb(5,3,'gap');

tb(3,4,'ppp');
tb(3,2,'pep');
tb(3,1,'pyr');

tb(4,1,'acc');
tb(4,2,'oaa');

% glu -> g6p:
% F = 1/2;
F = 1;
arr([1 3], [2 3], F); 

% glu -> pyr:
% F = 1/2;
F = 1;
arr([1 3], [3 1], F); 

% pep -> g6p:
% F = 1/2;
F = 1;
arr([3 2], [2 3], F);

% pep -> pyr:
% F = 1/2 + rxn.Flux( 19 )/rxn.Flux(1) ;
F = 1 + rxn.Flux( 19 )/rxn.Flux(1) ;
arr([3 2], [3 1], F);

% g6p -> f6p:
F = rxn.Flux(2)/rxn.Flux(1) ;
arr([2 3], [3 3], F); 

% g6p -> ppp:
F = rxn.Flux(4)/rxn.Flux(1) ;
arr([2 3], [3 4], F); 

% ppp -> f6p:
% F = (rxn.Flux(8)/2 + rxn.Flux(6)/2)/rxn.Flux(1);
F = (rxn.Flux(8) + rxn.Flux(6))/rxn.Flux(1);
arr([3 4], [3 3], F ); 

% ppp -> gap:
% F = rxn.Flux(8)/2/rxn.Flux(1);
F = rxn.Flux(8)/rxn.Flux(1);
arr([3 4], [5 3], F ); 

% f6p -> fdp:
F = rxn.Flux(5)/rxn.Flux(1);
arr([3 3], [4 3], F ); 

% fdp -> gap:
% F = (rxn.Flux(10)/2+rxn.Flux(12))/rxn.Flux(1);
F = (rxn.Flux(10)+rxn.Flux(12))/rxn.Flux(1);
arr([4 3], [5 3], F); 

% gap -> pep:
F = rxn.Flux(18)/rxn.Flux(1) ;
arr([5 3], [3 2], F);

% pep -> oaa:
F = rxn.Flux(20)/rxn.Flux(1);
arr([3 2], [4 2], F);

% pyr -> accoa:
F = rxn.Flux(24)/rxn.Flux(1);
arr([3 1], [4 1], F);


    function [] = tb(ix, iy, s)
        annotation('textbox', [xpos(ix) ypos(iy) TB.w TB.h], 'String', s,...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
    end


    function [] = arr( pos1, pos2, frac)
        if pos2(1) > pos1(1) && pos2(2) > pos1(2)
            % arrow points up and right
            x1 = xpos(pos1(1)) + TB.w + TB.bx;
            y1 = ypos(pos1(2)) + TB.h + TB.by;
            x2 = xpos(pos2(1)) - TB.bx;
            y2 = ypos(pos2(2)) - TB.by;
            tx = (x1+x2)/2 - TB.bx;
            ty = (y1+y2)/2 + TB.by/2/2;
            ha = 'right';
            va = 'bottom';
        elseif pos2(1) > pos1(1) && pos2(2) < pos1(2)
            % arrow points down and right
            x1 = xpos(pos1(1)) + TB.w + TB.bx;
            y1 = ypos(pos1(2)) - TB.by;
            x2 = xpos(pos2(1)) - TB.bx;
            y2 = ypos(pos2(2)) + TB.h + TB.by;
            tx = (x1+x2)/2 + TB.bx;
            ty = (y1+y2)/2 ;
            ha = 'left';
            va = 'middle';
        elseif pos2(1) < pos1(1) && pos2(2) < pos1(2)
            % arrow points down and left
            x1 = xpos(pos1(1)) - TB.bx;
            y1 = ypos(pos1(2)) - TB.by;
            x2 = xpos(pos2(1)) + TB.w + TB.bx;
            y2 = ypos(pos2(2)) + TB.h + TB.by;
            tx = (x1+x2)/2 - TB.bx;
            ty = (y1+y2)/2 ;
            ha = 'right';
            va = 'middle';            
        elseif pos2(1) < pos1(1) && pos2(2) > pos1(2)
            % arrow points up and left
            x1 = xpos(pos1(1)) - TB.bx;
            y1 = ypos(pos1(2)) + TB.h + TB.by;
            x2 = xpos(pos2(1)) + TB.w + TB.bx;
            y2 = ypos(pos2(2)) - TB.by;
            tx = (x1+x2)/2 - TB.bx;
            ty = (y1+y2)/2 ;
            ha = 'left';
            va = 'middle';                        
        elseif pos2(1) > pos1(1) && pos2(2) == pos1(2)
            % arrow points right
            x1 = xpos(pos1(1)) + TB.w + TB.bx;
            y1 = ypos(pos1(2)) + TB.h/2;
            x2 = xpos(pos2(1)) - TB.bx;
            y2 = ypos(pos1(2)) + TB.h/2;
            tx = (x1+x2)/2 - TB.w/4;
            ty = (y1+y2)/2 + TB.by;
            ha = 'center';
            va = 'bottom';                                    
        elseif pos2(1) < pos1(1) && pos2(2) == pos1(2)
            % arrow points left
            x1 = xpos(pos1(1)) - TB.bx;
            y1 = ypos(pos1(2)) + TB.h/2;
            x2 = xpos(pos2(1)) + TB.w + TB.bx;
            y2 = ypos(pos1(2)) + TB.h/2;
            tx = (x1+x2)/2 - TB.w/4;
            ty = (y1+y2)/2 + TB.by;
            ha = 'center';
            va = 'bottom';                                    
        elseif pos2(1) == pos1(1) && pos2(2) > pos1(2)
            % arrow points up 
            x1 = xpos(pos1(1)) + TB.w/2;
            y1 = ypos(pos1(2)) + TB.h + TB.by;
            x2 = xpos(pos1(1)) + TB.w/2;
            y2 = ypos(pos2(2)) - TB.by;
            tx = (x1+x2)/2 + TB.bx ;
            ty = (y1+y2)/2 - TB.h/2/2;
            ha = 'left';
            va = 'bottom';                                                
        elseif pos2(1) == pos1(1) && pos2(2) < pos1(2)
            % arrow points down
            x1 = xpos(pos1(1)) + TB.w/2;
            y1 = ypos(pos1(2)) - TB.by;
            x2 = xpos(pos1(1)) + TB.w/2;
            y2 = ypos(pos2(2)) + TB.h + TB.by;
            tx = (x1+x2)/2 + TB.bx ;
            ty = (y1+y2)/2 - TB.h/2/2;
            ha = 'left';
            va = 'bottom';                                                
        end
        
        s = sprintf('%.0f%%', frac*100);
        
        if frac > 0
        annotation('arrow', [x1 x2], [y1 y2], 'LineWidth', frac*A.w100);
        end
        annotation('textbox', [tx ty TB.w/2 TB.h/2], 'String', s,...
            'HorizontalAlignment', ha, 'VerticalAlignment', 'bottom',...
            'EdgeColor', 'none');
    end
end
