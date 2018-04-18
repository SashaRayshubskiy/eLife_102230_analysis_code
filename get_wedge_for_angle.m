function [wedge_id] = get_wedge_for_angle( angle_vec )

theta_deg = rad2deg(atan2( angle_vec(2) , angle_vec(1) ));

if(theta_deg < 0)
    theta_deg = theta_deg + 360.0;
end

if( ( theta_deg < 90 ) && ( theta_deg >= 67.5 ))
    wedge_id = 1;
elseif( ( theta_deg < 67.5 ) && ( theta_deg >= 45.0 ))
    wedge_id = 2;    
elseif( ( theta_deg < 45.0 ) && ( theta_deg >= 22.5 ))
    wedge_id = 3;    
elseif( ( theta_deg < 22.5 ) && ( theta_deg >= 0 ))
    wedge_id = 4;    
elseif( ( theta_deg < 360 ) && ( theta_deg >= 337.5 ))
    wedge_id = 5;    
elseif( ( theta_deg < 337.5 ) && ( theta_deg >= 315 ))
    wedge_id = 6;    
elseif( ( theta_deg < 315 ) && ( theta_deg >= 292.5 ))
    wedge_id = 7;    
elseif( ( theta_deg < 292.5 ) && ( theta_deg >= 270 ))
    wedge_id = 8;    
elseif( ( theta_deg < 270 ) && ( theta_deg >= 247.5 ))
    wedge_id = 9;    
elseif( ( theta_deg < 247.5 ) && ( theta_deg >= 225 ))
    wedge_id = 10;    
elseif( ( theta_deg < 225 ) && ( theta_deg >= 202.5 ))
    wedge_id = 11;    
elseif( ( theta_deg < 202.5 ) && ( theta_deg >= 180.0 ))
    wedge_id = 12;
elseif( ( theta_deg < 180 ) && ( theta_deg >= 157.5 ))
    wedge_id = 13;    
elseif( ( theta_deg < 157.5 ) && ( theta_deg >= 135 ))
    wedge_id = 14;    
elseif( ( theta_deg < 135 ) && ( theta_deg >= 112.5 ))
    wedge_id = 15;    
elseif( ( theta_deg < 112.5 ) && ( theta_deg >= 90 ))
    wedge_id = 16;    
else
    disp(['ERROR: angle: ' num2str(theta_deg) ' not recognized.']);
    wedge_id = 1;
end


end


% target_wedge = floor(theta_deg / WEDGE_INCREMENT_IN_DEG) + 1;
% 
% % Wedge at [90:67.5] deg is 1, followed by counter-clockwise numbering
% offset = 4;
% 
% wedge_id = target_wedge + offset;
% 
% if(wedge_id > 16)
%     wedge_id = wedge_id - 16;
% end
