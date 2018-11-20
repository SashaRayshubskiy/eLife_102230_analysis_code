function  display_bump_dynamics_turning_ephys( roi_ts_per_trial, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid, pre_stim, stim_dur )

ac = get_analysis_constants;
post_stim = pre_stim + stim_dur;

nframes = size(roi_ts_per_trial,3);

t = [0:nframes-1]./VPS;
stim_t_vec = find( (t >= (pre_stim-0.02)) & (t <=  post_stim) );
stim_t = stim_t_vec(1);

% Calculate stability of bump dynamics before stim and separate trials 
% by bump jump to the left or right of the currently stable bump.
left_bump_jump = [];
left_bump_jump_yaw = [];
left_bump_jump_ephys = [];

right_bump_jump = [];
right_bump_jump_yaw = [];
right_bump_jump_ephys = [];

STABILITY_THRESHOLD = 1;
bump_jump_t = find( (t >= pre_stim) & (t <=  (pre_stim+1.0)) );

NUM_PB_GLOMERULI_HALF = 8;
PB_GLOMERULI_MID = 5;

PB_half = [ 1:NUM_PB_GLOMERULI_HALF ];

stability_markers = [];
for tr = 1:size(roi_ts_per_trial,1)
    
    cur_trial = squeeze(roi_ts_per_trial(tr,:,:));
    
    pre_stim_bump_dyn = cur_trial( PB_half, 1:stim_t );
    
    [ val, ind ] = max( pre_stim_bump_dyn );
    
    stability_marker = std(ind,1);
    stability_markers(end+1) = stability_marker;
    
    % Check if ind is stable
    if( stability_marker <  STABILITY_THRESHOLD )
        % Check if the bump jump is to the left or the right of the sBump
        % (stable bump)
        % classify trials
        
        sBump_inx = ceil(mean(ind));
        
        [bump_jump_val, bump_jump_idx] = max(squeeze(mean(cur_trial(:,bump_jump_t),2)));
        bump_jump_idx = ceil( bump_jump_idx );
        
        sb_offset = PB_GLOMERULI_MID - sBump_inx;
        bjump_offset = mod((bump_jump_idx + sb_offset), NUM_PB_GLOMERULI_HALF );
        jump_delta = PB_GLOMERULI_MID - bjump_offset;
                
        [dummy, bump_locations] = max( cur_trial );
        bump_locations_shifted = mod((bump_locations + sb_offset), NUM_PB_GLOMERULI_HALF );
        
        if( ( jump_delta >= -4 ) && ( jump_delta <= -2 ) )
            % Add to left
            left_bump_jump(end+1,:)       = bump_locations_shifted;
            left_bump_jump_yaw(end+1,:)   = squeeze(bdata_vel{ 1 }( tr, ac.VEL_YAW, : ));            
            left_bump_jump_ephys(end+1,:) = squeeze(ephys_data{ 1 }( tr, : ));
        elseif( ( jump_delta >= 2 ) && ( jump_delta <= 4 ) )
            % Add to right
            right_bump_jump(end+1,:)       = bump_locations_shifted;            
            right_bump_jump_yaw(end+1,:)   = squeeze(bdata_vel{ 1 }( tr, ac.VEL_YAW, : ));
            right_bump_jump_ephys(end+1,:) = squeeze(ephys_data{ 1 }( tr, : ));
        end        
    end    
end

figure; 
histogram( stability_markers, 20 );

% Plot left/right bump jump and yaw
f = figure;

ax(1) = subplot(3,1,1);
hold on;
plot( t, squeeze(mean( left_bump_jump )), 'color', rgb('FireBrick' ), 'LineWidth', 2.0, 'DisplayName', ['Num trials: ' num2str(size(left_bump_jump,1))] );
plot( t, squeeze(mean( right_bump_jump )), 'color', rgb('SeaGreen' ), 'LineWidth', 2.0, 'DisplayName', ['Num trials: ' num2str(size(right_bump_jump,1))] );
ylabel('Bump rep');

legend();

ax(2) = subplot(3,1,2);
hold on;
plot( bdata_vel_time, squeeze(mean( left_bump_jump_yaw )),  'color', rgb('FireBrick' ), 'LineWidth', 2.0 );
plot( bdata_vel_time, squeeze(mean( right_bump_jump_yaw )), 'color', rgb('SeaGreen' ), 'LineWidth', 2.0 );
ylabel('Yaw (au)');

ax(3) = subplot(3,1,3);
hold on;
plot( ephys_time, squeeze(mean( left_bump_jump_ephys )),  'color', rgb('FireBrick' ), 'LineWidth', 2.0 );
plot( ephys_time, squeeze(mean( right_bump_jump_ephys )), 'color', rgb('SeaGreen' ), 'LineWidth', 2.0 );

ylabel('Vm (mV)');
xlabel('Time (s)');

linkaxes(ax,'x');

xlim([0 t(end)]);
axis tight;

saveas(f,[analysis_path '/bump_dynamics_vs_behaviour_vs_ephys_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/bump_dynamics_vs_behaviour_vs_ephys_' num2str(sid) '.png']);

end

