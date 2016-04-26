function display_single_trial_trajectories( sid, bdata_vel_time, traj, analysis_path )

ac = get_analysis_constants;
settings = sensor_settings;

prestim = settings.pre_stim;
stim    = settings.stim;

stim_t = find((bdata_vel_time > prestim) & (bdata_vel_time < (prestim+stim)));

for trial_type = 1:length(traj)

    f = figure; 
    
    for trial = 1:size( traj{ trial_type }, 1 )
        disp_x = squeeze(traj{ trial_type }( trial, 1, :));
        disp_y = squeeze(traj{ trial_type }( trial, 2, :));
        
        hold on;
        plot(disp_x, disp_y, 'color', rgb('Purple'), 'LineWidth', 1.0 );
        plot(disp_x(stim_t), disp_y(stim_t), 'color', rgb('Crimson'), 'LineWidth', 2.0 );
    end
    
    xlim([-0.1 0.3]);
    ylim([-0.05 1.5]);
    title(ac.task_str{ trial_type });

    fname = [ analysis_path '/sid_' num2str( sid ) '_single_trial_traj_' ac.task_str{trial_type}];
    saveas(f, [ fname '.fig'] );
    saveas(f, [ fname '.png'] );
end

end

