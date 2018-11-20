function  display_bump_dynamics_turning_ephys_var_stim( roi_ts_per_trial, pico_stim_data, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid )

ac = get_analysis_constants;
settings = sensor_settings;

nframes = size(roi_ts_per_trial,3);
t = [0:nframes-1]./VPS;

num_trials = size(pico_stim_data{1},1);

bump_around_stim  = [];
yaw_around_stim   = [];
ephys_around_stim = [];

before_stim_t = 2.0; % seconds
after_stim_t  = 7.0;

before_frames_bump = before_stim_t * VPS;
after_frames_bump = after_stim_t * VPS;

before_frames_yaw = before_stim_t * settings.sensorPollFreq;
after_frames_yaw = after_stim_t * settings.sensorPollFreq;

before_frames_ephys = before_stim_t * settings.sampRate;
after_frames_ephys = after_stim_t * settings.sampRate;

dt_bump  = 2*(1.0/VPS);
dt_yaw   = 2*(1.0/settings.sensorPollFreq);
dt_ephys = 2*(1.0/settings.sampRate);

for tr = 1:num_trials
       
    cur_stims = pico_stim_data{1}(tr,:);

    cur_bump_dynamics  = squeeze(roi_ts_per_trial( tr, :, :));
    cur_yaw_dynamics   = squeeze(bdata_vel{ 1 }( tr, ac.VEL_YAW, : ));            
    cur_ephys_dynamics = squeeze(ephys_data{ 1 }( tr, : ));     
    
    stim_ids = find(diff(cur_stims) > 4);
    stim_times = ephys_time( stim_ids );
    
    % Process each stim
    for st = 1:length(stim_times)
        cur_stim_time = stim_times(st);
                
        % Extract bump dynamics 
        bump_stim_time_idx = find( (t> (cur_stim_time-dt_bump)) & ((t<(cur_stim_time+dt_bump)) ));
        if( length(bump_stim_time_idx) == 0 )
            error(['Could not find a bump time index for time: ' num2str(cur_stim_time)]);
        end        
        
        cur_bump_stim_time_idx = bump_stim_time_idx(1);
        
        bump_start_idx = cur_bump_stim_time_idx - before_frames_bump;
        bump_end_idx = cur_bump_stim_time_idx + after_frames_bump;
        
        bump_around_stim(end+1,:,:) = cur_bump_dynamics(:, bump_start_idx:bump_end_idx );
        
        % Extract yaw dynamics
        yaw_stim_time_idx = find( (bdata_vel_time> (cur_stim_time-dt_yaw)) & ((bdata_vel_time < (cur_stim_time+dt_yaw)) ));
        if( length(yaw_stim_time_idx) == 0 )
            error(['Could not find a yaw time index for time: ' num2str(cur_stim_time)]);
        end        
        
        cur_yaw_stim_time_idx = yaw_stim_time_idx(1);
        
        yaw_start_idx = cur_yaw_stim_time_idx - before_frames_yaw;
        yaw_end_idx = cur_yaw_stim_time_idx + after_frames_yaw;
        
        yaw_around_stim(end+1,:) = cur_yaw_dynamics( yaw_start_idx:yaw_end_idx );
        
        % Extract ephys dynamics                
        ephys_stim_time_idx = find( (ephys_time > (cur_stim_time-dt_ephys)) & ((ephys_time < (cur_stim_time+dt_ephys)) ));
        if( length(ephys_stim_time_idx) == 0 )
            error(['Could not find a ephys time index for time: ' num2str(cur_stim_time)]);
        end        
        
        cur_ephys_stim_time_idx = ephys_stim_time_idx(1);
        
        ephys_start_idx = cur_ephys_stim_time_idx - before_frames_ephys;
        ephys_end_idx = cur_ephys_stim_time_idx + after_frames_ephys;
        
        ephys_around_stim(end+1,:) = cur_ephys_dynamics( ephys_start_idx:ephys_end_idx );
    end
end

% Plot left/right bump jump and yaw
bump_t_base  = ([0:size(bump_around_stim,3)-1] ./ VPS)                    - before_stim_t;
yaw_t_base   = ([0:size(yaw_around_stim,2)-1] ./ settings.sensorPollFreq) - before_stim_t;
ephys_t_base = ([0:size(ephys_around_stim,2)-1] ./ settings.sampRate)     - before_stim_t;

f = figure;

num_stims = size(ephys_around_stim,1);
cm = colormap(jet(num_stims));

active_gloms = [1:8];

for st = 1:num_stims
    
    cur_bump_gloms = squeeze(bump_around_stim(st,active_gloms,:));

    % cur_bump_gloms(find(cur_bump_gloms < 0.0 )) = 0.0;
    
    % [dummy, bump_motion_ids] = max(cur_bump_gloms);
    % bump_motion_ids = get_weigh_avg_bump_pos( cur_bump_gloms );
    bump_motion_ids = get_radial_weighted_avg_bump_pos( cur_bump_gloms );
    
    ax(1) = subplot(3,1,1);
    hold on;
    plot( bump_t_base, squeeze(bump_motion_ids) ,'color', cm(st,:), 'LineWidth', 1.0, 'DisplayName', ['Num stims: ' num2str(num_stims)] );
    %ylim([1 10]);
    ylabel('Bump rep');
        
    ax(2) = subplot(3,1,2);
    hold on;
    plot( yaw_t_base, squeeze(yaw_around_stim(st,:)), 'color', cm(st,:), 'LineWidth', 1.0 );
    ylabel('Yaw (au)');
    
    ax(3) = subplot(3,1,3);
    hold on;
    plot( ephys_t_base, squeeze(ephys_around_stim(st,:)), 'color', cm(st,:), 'LineWidth', 1.0 );
    
    ylabel('Vm (mV)');
    xlabel('Time (s)');
    xlim([0 t(end)]);
    axis tight;
end

ax(1) = subplot(3,1,1);
legend();

linkaxes(ax,'x');

saveas(f,[analysis_path '/var_stim_bump_dynamics_vs_behaviour_vs_ephys_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/var_stim_bump_dynamics_vs_behaviour_vs_ephys_' num2str(sid) '.png']);



end

