function [ accel ] = get_accel( t_data, vel_data )

% Assumes input is [velocity events, velocity times]
if( size(vel_data,1) == 1 )
    avg_vel = vel_data;
else
    avg_vel = mean(vel_data);
end

% Fit a slope 
%p = polyfit(t_data, avg_vel, 1);    
%accel = p(1);

[fobject,gof] = fit( t_data', avg_vel', 'poly1');

% RSQ_THRESHOLD = 0.5;
% if( gof.rsquare > RSQ_THRESHOLD )
%     accel = fobject.p1;
% else
%     accel = 0;    
% end

if(mean(avg_vel) < 0)
    accel = 0;
else
    accel = fobject.p1;
end
        
end

