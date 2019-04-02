function [ bump_motion_ids ] = assess_and_fix_bump_quality( cur_bump_in, VPS )

% If there are sufficient consecutive gaps in bump quality then disqualify this trial,
% Otherwise, fill in the gaps.

TIME_GAP_WITH_NO_BUMP    = 1.0;
SAMPLES_GAP_WITH_NO_BUMP = floor( TIME_GAP_WITH_NO_BUMP * VPS );

last_bump_idx = 1;
NAN_cnt = 0;
nan_buffer = [];
bump_motion_ids = NaN(1, length( cur_bump_in ) );
bump_motion_ids(1) = cur_bump_in(1);

n_frames = length( cur_bump_in );

for ts = 2:n_frames
    
    cur_bump_loc = cur_bump_in(ts);
    
    if( isnan(cur_bump_loc) == 1 )
        nan_buffer(end+1) = ts;
        NAN_cnt = NAN_cnt + 1;
        
        if( NAN_cnt == SAMPLES_GAP_WITH_NO_BUMP )
            bump_motion_ids = [];
            return;
        end
    else        
        bump_motion_ids(ts) = cur_bump_loc;
        
        if( length( nan_buffer ) > 0 )
            % Interpolate between last bump and current
            nan_interp = interp1( [cur_bump_in(last_bump_idx) cur_bump_loc], [1:1/(length(nan_buffer)+1):2] );
            bump_motion_ids(nan_buffer) = nan_interp(2:end-1);
        end
        
        nan_buffer = [];
        NAN_cnt = 0;
        last_bump_idx = ts;
    end
end

end

