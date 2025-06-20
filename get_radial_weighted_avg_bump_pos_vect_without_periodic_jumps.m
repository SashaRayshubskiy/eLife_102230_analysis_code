function  [ smoothed_bump, bump_motion_ids, bump_motion_ids_unwrapped, bump_in_deg_unwrapped, vect_strength_check ] = get_radial_weighted_avg_bump_pos_vect_without_periodic_jumps( bump_gloms, BUMP_QUALITY_THRESHOLD )

% assumes input { glomerulus, time }

n_gloms  = size( bump_gloms, 1 );
n_frames = size( bump_gloms, 2 );

bump_motion_ids = zeros(1, n_frames);

UPSAMPLE_FACTOR = 100;

WEDGE_CNT = n_gloms * UPSAMPLE_FACTOR;
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

tmp_bump_motion_ids = zeros(1,n_frames);
bump_motion_ids = zeros(1,n_frames);

tmp_bump_motion_ids_unwrapped = zeros(1,n_frames);
bump_motion_ids_unwrapped     = zeros(1,n_frames);
tmp_bump_motion_deg_unwrapped = zeros(1,n_frames);

if( nargin == 2 )
    NO_BUMP_THRESHOLD = BUMP_QUALITY_THRESHOLD;
else   
    NO_BUMP_THRESHOLD = 0.2;
end

BUMP_LOCATION_OFFSET = 1;

stop_me = 0;

vect_strength_check = [];

prev_deg = -1;

for ts = 1:n_frames
    
    if(ts == 85)
        stop_me = stop_me + 1;
    end
    
    cur_ts_df_f = squeeze( bump_gloms(:,ts));
    
    cur_ts_df_f_tmp = cur_ts_df_f;
    cur_ts_df_f_tmp(end+1) = cur_ts_df_f(1);
    
    cur_EB_bump_upsample = interp1( cur_ts_df_f_tmp, [ 1 : 1/UPSAMPLE_FACTOR : length(cur_ts_df_f_tmp)-(1/UPSAMPLE_FACTOR) ] )';
    
    % cur_EB_bump_up_han = hanningsmooth(cur_EB_bump_upsample, 200);
    cur_EB_bump_up_han = hanningsmooth_circ(cur_EB_bump_upsample, 200);
    
    smoothed_bump( :, ts ) = cur_EB_bump_up_han;
    
    cur_EB_bump_final = cur_EB_bump_up_han;
    
    if 0
        % hold on;
        hold on;
        % plot( cur_ts_df_f );
        plot( cur_EB_bump_upsample );
        plot( hanningsmooth(cur_EB_bump_upsample, 200), '--' );
        ylim([-0.5 4]);
        waitforbuttonpress;
        cla();
    end
    
    cur_wedge_df_f_sum = sum( cur_EB_bump_final );
    
    PVA_ts_x_0 = squeeze(wedge_unit_vectors(:,1));
    PVA_ts_x_1 = PVA_ts_x_0 .* cur_EB_bump_final;
    PVA_ts_x_sum = sum(PVA_ts_x_1,1);
    PVA_ts_x = PVA_ts_x_sum / cur_wedge_df_f_sum;
    
    PVA_ts_y_0 = squeeze(wedge_unit_vectors(:,2));
    PVA_ts_y_1 = PVA_ts_y_0 .* cur_EB_bump_final;
    PVA_ts_y_sum = sum(PVA_ts_y_1,1);
    PVA_ts_y = PVA_ts_y_sum / cur_wedge_df_f_sum;
    
    vect_strength = sqrt( PVA_ts_x*PVA_ts_x + PVA_ts_y*PVA_ts_y );

    vect_strength_check(ts) = vect_strength;
    
    if( vect_strength < NO_BUMP_THRESHOLD )
        tmp_bump_motion_ids(ts) = NaN;
        tmp_bump_motion_ids_unwrapped(ts) = NaN;
        cur_deg = prev_deg;
    else

        theta_deg = wrapTo360( rad2deg(atan2( PVA_ts_y, PVA_ts_x )));
        % theta_deg = rad2deg(atan2( PVA_ts_y, PVA_ts_x ));        
        bump_location = ((theta_deg / WEDGE_INCREMENT_IN_DEG) / UPSAMPLE_FACTOR) + BUMP_LOCATION_OFFSET;
        %bump_location = theta_deg;        
        tmp_bump_motion_ids(ts) = bump_location;

        theta_deg_unwrapped = rad2deg(atan2( PVA_ts_y, PVA_ts_x ));                
        bump_location_unwrapped = ((theta_deg_unwrapped / WEDGE_INCREMENT_IN_DEG) / UPSAMPLE_FACTOR) + BUMP_LOCATION_OFFSET;
        tmp_bump_motion_ids_unwrapped(ts) = bump_location_unwrapped;
        
        if( prev_deg ~= -1 )
            delta_deg = shortest_distance_in_deg( prev_deg, theta_deg_unwrapped );
            cur_deg = cur_deg + delta_deg;
        else
            % This should be set once at the start
            cur_deg = theta_deg_unwrapped;
        end        
    end
    
    tmp_bump_motion_deg_unwrapped(ts) = cur_deg;
    prev_deg = cur_deg;
end

% Check distribution of vector strengths
% figure;
% histogram( vector_strengths_test );

bump_motion_ids           = tmp_bump_motion_ids;
bump_motion_ids_unwrapped = tmp_bump_motion_ids_unwrapped;
bump_in_deg_unwrapped     = tmp_bump_motion_deg_unwrapped;

% Fix any reasonable discountinueties in bump motion
if 0
    %
    if( tmp_bump_motion_ids(1) == NaN )
        % bump_motion_ids = [];
        % return;
        bump_motion_ids(1) = 0;
    end
    
    MAX_NAN_TIME_POINTS = 18;
    last_bump_idx = 1;
    NAN_cnt = 0;
    nan_buffer = [];
    bump_motion_ids(1) = tmp_bump_motion_ids(1);
    
    for ts = 2:n_frames
        
        cur_bump_loc = tmp_bump_motion_ids(ts);
        
        if( isnan(cur_bump_loc) == 1 )
            nan_buffer(end+1) = ts;
            NAN_cnt = NAN_cnt + 1;
            
            if( NAN_cnt == MAX_NAN_TIME_POINTS )
                bump_motion_ids = [];
                return;
            end
        else
            
            bump_motion_ids(ts) = tmp_bump_motion_ids(ts);
            if( length( nan_buffer ) > 0 )
                % Interpolate between last bump and current
                nan_interp = interp1( [tmp_bump_motion_ids(last_bump_idx) cur_bump_loc], [1:1/(length(nan_buffer)+1):2] );
                bump_motion_ids(nan_buffer) = nan_interp(2:end-1);
            end
            
            nan_buffer = [];
            NAN_cnt = 0;
            last_bump_idx = ts;
        end
    end
    
    if 0
        figure;
        hold on;
        plot(tmp_bump_motion_ids);
        plot(bump_motion_ids+20);
        waitforbuttonpress();
    end
end

end

