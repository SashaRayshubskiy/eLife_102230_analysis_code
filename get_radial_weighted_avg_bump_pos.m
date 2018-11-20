function  bump_motion_ids = get_radial_weighted_avg_bump_pos( bump_gloms )

% assumes input { glomerulus, time }

n_gloms  = size( bump_gloms, 1 );
n_frames = size( bump_gloms, 2 );

bump_motion_ids = zeros(1, n_frames);

WEDGE_CNT = 8;
wedge_unit_vectors = zeros(WEDGE_CNT, 2);
starting_angle = 0;
WEDGE_INCREMENT_IN_DEG = 360.0/WEDGE_CNT;

% Create the unit vectors
for w = 1:WEDGE_CNT    
    
    cur_ang = starting_angle + (WEDGE_INCREMENT_IN_DEG/2.0);
    
    % Unit vector 
    x = cos(deg2rad(cur_ang));
    y = sin(deg2rad(cur_ang));
    
    wedge_unit_vectors(w,1) = x;
    wedge_unit_vectors(w,2) = y;
       
    starting_angle = starting_angle + WEDGE_INCREMENT_IN_DEG;        
end

for ts = 1:n_frames
    
    cur_ts_df_f = squeeze( bump_gloms(:,ts));
    
    cur_wedge_df_f_sum = sum(cur_ts_df_f);
    
    PVA_ts_x_0 = squeeze(wedge_unit_vectors(:,1));
    PVA_ts_x_1 = PVA_ts_x_0 .* cur_ts_df_f;
    PVA_ts_x_sum = sum(PVA_ts_x_1,1);
    PVA_ts_x = PVA_ts_x_sum / cur_wedge_df_f_sum;
    
    PVA_ts_y_0 = squeeze(wedge_unit_vectors(:,2));
    PVA_ts_y_1 = PVA_ts_y_0 .* cur_ts_df_f;
    PVA_ts_y_sum = sum(PVA_ts_y_1,1);
    PVA_ts_y = PVA_ts_y_sum / cur_wedge_df_f_sum;
    
    theta_deg = wrapTo360( rad2deg(atan2( PVA_ts_y, PVA_ts_x )));
    
    %bump_location = theta_deg / WEDGE_INCREMENT_IN_DEG;
    bump_location = theta_deg;
    
    bump_motion_ids(ts) = bump_location;
    
end

end

