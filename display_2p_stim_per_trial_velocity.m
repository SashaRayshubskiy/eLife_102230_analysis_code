function display_2p_stim_per_trial_velocity( sid, bdata_raw, bdata_vel, bdata_vel_time, analysis_path, with_single_trials )

ac = get_analysis_constants;
settings = sensor_settings;

f = figure;

one_trial_bdata = squeeze(bdata_raw{ ac.LEFT}(1,:,:));
left_odor_stim  = squeeze(one_trial_bdata(:,6));

first_stim = find((left_odor_stim > 2.0), 1, 'first') ./ settings.sampRate;
last_stim = find((left_odor_stim > 2.0), 1, 'last') ./ settings.sampRate;
prestim = 3.0;
stim = 1.0;

pre_stim_t = find(bdata_vel_time < prestim);
stim_t = find((bdata_vel_time >= prestim) & (bdata_vel_time <= (prestim+stim)));

for trial_type = 1:length(bdata_vel)
    subplot(length(bdata_vel),1,trial_type);

    cur_trial_cnt = size(bdata_vel{trial_type}, 1); 
    avg_pre_stim_fwd_vel = zeros(1,cur_trial_cnt);
    avg_stim_yaw_vel = zeros(1,cur_trial_cnt);

    for tid = 1:cur_trial_cnt
        avg_pre_stim_fwd_vel( tid ) = mean( squeeze(bdata_vel{ trial_type }( tid, ac.VEL_FWD, pre_stim_t )) ); 
        avg_stim_yaw_vel( tid ) = mean( squeeze(bdata_vel{ trial_type }( tid, ac.VEL_YAW, stim_t )) ); 
    end

    colormap jet;
    bar([avg_pre_stim_fwd_vel', avg_stim_yaw_vel'], 2.0, 'EdgeColor', 'none');
    
    if( trial_type == 1 )
        legend('Avg pre-stim fwd vel','Avg stim yaw vel');
    end

    xlim([0 cur_trial_cnt]);
    xlabel('Trial #', 'FontSize', 14);
    ylabel('Velocity (au/s)', 'FontSize', 14);
    title([ac.task_str{trial_type} ': ' num2str(cur_trial_cnt)], 'FontSize', 14);
end

saveas(f, [analysis_path '/trial_by_trial_vel_' num2str( sid ) '.fig']);
saveas(f, [analysis_path '/trial_by_trial_vel_' num2str( sid ) '.png']);

end

