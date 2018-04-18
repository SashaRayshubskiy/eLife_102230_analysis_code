function [rois, wedge_unit_vectors] = get_EPG_wedge_rois(diameter_in, diameter_out, xCenter, yCenter)

% diameter = 4;
% diameter_1 = 22;
% xCenter = 16;
% yCenter = 13;

diameter = diameter_in;
diameter_1 = diameter_out;

radius = diameter / 2;
xMin= xCenter - radius;
yMin = yCenter - radius;
hEllipse = imellipse(gca,[xMin, yMin, diameter, diameter]);

radius_1 = diameter_1 / 2;
xMin_1 = xCenter - radius_1;
yMin_1 = yCenter - radius_1;
hEllipse = imellipse(gca,[xMin_1, yMin_1, diameter_1, diameter_1]);

% Plot center
plot(xCenter, yCenter, 'r+', 'LineWidth', 2, 'MarkerSize', 5);
grid on;

WEDGE_CNT = 16;
wedge_unit_vectors = zeros(WEDGE_CNT, 2);
starting_angle = 90;
WEDGE_INCREMENT_IN_DEG = 360.0/WEDGE_CNT;

for w = 1:WEDGE_CNT    
    
    cur_ang = starting_angle - (WEDGE_INCREMENT_IN_DEG/2.0);
    
    % Unit vector 
    x = cos(deg2rad(cur_ang));
    y = sin(deg2rad(cur_ang));
    
    wedge_unit_vectors(w,1) = x;
    wedge_unit_vectors(w,2) = y;
    
    SCALING_FACT = 10;
    plot([xCenter, xCenter+SCALING_FACT*x], [yCenter yCenter+SCALING_FACT*y], 'w', 'LineWidth', 1.0 );
    
    starting_angle = starting_angle - WEDGE_INCREMENT_IN_DEG;        
end

slopes = [ inf, sin(deg2rad(3*45/2.0))/cos(deg2rad(3*45/2.0)), 1, sin(deg2rad(45/2.0))/cos(deg2rad(45/2.0)), 0, ... 
               -sin(deg2rad(45/2.0))/cos(deg2rad(45/2.0)), -1, -sin(deg2rad(3*45/2.0))/cos(deg2rad(3*45/2.0)) ];

intercepts = zeros(1, length(slopes));

intercepts(1) = xCenter;

for i = 2:length(slopes)
    intercepts(i) = yCenter - slopes(i)*xCenter;
end

% Create first 4 ROIs that will be reflected 3 times
x_outer = cell(1,4);
y_outer = cell(1,4);
x_inner = cell(1,4);
y_inner = cell(1,4);
for i = 1:5
    [ x_outer{i}, y_outer{i} ] = linecirc( slopes(i), intercepts(i), xCenter, yCenter, radius   );
    [ x_inner{i}, y_inner{i} ] = linecirc( slopes(i), intercepts(i), xCenter, yCenter, radius_1 );
end

wedge_rois = cell(1,WEDGE_CNT);
for i = 1:4
    cur_x = [ x_outer{i}(1) x_inner{i}(1) x_inner{i+1}(1) x_outer{i+1}(1) ];
    cur_y = [ y_outer{i}(1) y_inner{i}(1) y_inner{i+1}(1) y_outer{i+1}(1) ];
    wedge_rois{i} = {cur_x, cur_y};
end

SLOPE_ID = 7;

for i = 1:4    
    [ cur_x_1_16(i), cur_y_1_16(i) ] = reflect_in_vert_line([ wedge_rois{1}{1}(i) wedge_rois{1}{2}(i) ], xCenter );
    [ cur_x_1_8(i), cur_y_1_8(i) ] = reflect_in_horz_line([ wedge_rois{1}{1}(i) wedge_rois{1}{2}(i) ], yCenter );
    [ cur_x_1_9(i), cur_y_1_9(i) ] = reflect_in_vert_line([ cur_x_1_8(i), cur_y_1_8(i) ], xCenter );

    [ cur_x_2_15(i), cur_y_2_15(i) ] = reflect_in_vert_line([ wedge_rois{2}{1}(i) wedge_rois{2}{2}(i) ], xCenter );
    [ cur_x_2_7(i), cur_y_2_7(i) ] = reflect_in_horz_line([ wedge_rois{2}{1}(i) wedge_rois{2}{2}(i) ], yCenter );
    [ cur_x_2_10(i), cur_y_2_10(i) ] = reflect_in_vert_line([ cur_x_2_7(i), cur_y_2_7(i) ], xCenter );

    [ cur_x_3_14(i), cur_y_3_14(i) ] = reflect_in_vert_line([ wedge_rois{3}{1}(i) wedge_rois{3}{2}(i) ], xCenter );
    [ cur_x_3_6(i), cur_y_3_6(i) ] = reflect_in_horz_line([ wedge_rois{3}{1}(i) wedge_rois{3}{2}(i) ], yCenter );
    [ cur_x_3_11(i), cur_y_3_11(i) ] = reflect_in_vert_line([ cur_x_3_6(i), cur_y_3_6(i) ], xCenter );

    [ cur_x_4_13(i), cur_y_4_13(i) ] = reflect_in_vert_line([ wedge_rois{4}{1}(i) wedge_rois{4}{2}(i) ], xCenter );
    [ cur_x_4_5(i), cur_y_4_5(i) ] = reflect_in_horz_line([ wedge_rois{4}{1}(i) wedge_rois{4}{2}(i) ], yCenter );
    [ cur_x_4_12(i), cur_y_4_12(i) ] = reflect_in_vert_line([ cur_x_4_5(i), cur_y_4_5(i) ], xCenter );
end

wedge_rois{5} = { cur_x_4_5, cur_y_4_5 };
wedge_rois{6} = { cur_x_3_6, cur_y_3_6 };
wedge_rois{7} = { cur_x_2_7, cur_y_2_7 };
wedge_rois{8} = { cur_x_1_8, cur_y_1_8 };
wedge_rois{9} = { cur_x_1_9, cur_y_1_9 };
wedge_rois{10} = { cur_x_2_10, cur_y_2_10 };
wedge_rois{11} = { cur_x_3_11, cur_y_3_11 };
wedge_rois{12} = { cur_x_4_12, cur_y_4_12 };
wedge_rois{13} = { cur_x_4_13, cur_y_4_13 };
wedge_rois{14} = { cur_x_3_14, cur_y_3_14 };
wedge_rois{15} = { cur_x_2_15, cur_y_2_15 };
wedge_rois{16} = { cur_x_1_16, cur_y_1_16 };

% Close the rectangle 
for w = 1:WEDGE_CNT
    xv = wedge_rois{w}{1};
    yv = wedge_rois{w}{2};

    % Close the rectangle 
    wedge_rois{w}{1} = [ xv xv(1) ];
    wedge_rois{w}{2} = [ yv yv(1) ];
end

rois = wedge_rois;

if 0
colorindex = 0;
ac = get_analysis_constants;
order    = ac.order;

for w = 1:16
    xv = wedge_rois{w}{1};
    yv = wedge_rois{w}{2};
    
    % Close the rectangle 
    xv = [ xv xv(1) ];
    yv = [ yv yv(1) ];
    
    currcolor    = order(1+mod(colorindex,size(order,1)),:);

    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    
    colorindex = colorindex + 1;
end
end

end

% reflection in a vertical line defined by xline
function [xr, yr] = reflect_in_vert_line(P, xline)
    xr = (xline - (P(1) - xline));
    yr = P(2);
end

function [xr, yr] = reflect_in_horz_line(P, yline)
    xr = P(1);
    yr = (yline - (P(2) - yline));
end

function [xr, yr] = reflect_in_x_y_line(P, n)

xr = P(2);
yr = P(1) - n;

end
    
function [xr, yr] = reflect_in_line(P, m, n)

    % line of symmetry: y = m*x + n;
    Md = zeros(2,1); % Middle point between given point and its symmetric
    Md(1) = (P(1) + m*P(2) - m*n)/(m^2 + 1);
    Md(2) = m*Md(1) + n;
    S = 2*Md - P; % symmetric point of P about given line
    xr = S(1);
    yr = S(2);
end
