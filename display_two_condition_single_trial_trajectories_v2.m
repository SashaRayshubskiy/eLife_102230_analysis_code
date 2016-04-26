function display_two_condition_single_trial_trajectories_v2( sid, condition_trials, condition_trials_str, bdata_vel_time, traj, analysis_path )

ac = get_analysis_constants;
settings = sensor_settings;

prestim = settings.pre_stim;
stim    = settings.stim;

stim_t = find((bdata_vel_time > prestim) & (bdata_vel_time < (prestim+stim)));

for trial_type = 1:size( condition_trials, 1 )
    
    f = figure; 

    for cond = 1:size( condition_trials, 2 )    
        
        subplot(1,2,cond);
        
        for trial = 1:length(condition_trials{trial_type, cond})
            
            cur_trial = condition_trials{trial_type,cond}(trial);
            
            disp_x = squeeze(traj{trial_type}( cur_trial, 1, :));
            disp_y = squeeze(traj{trial_type}( cur_trial, 2, :));
            
            hold on;
            plot(disp_x-disp_x(stim_t(1)), disp_y-disp_y(stim_t(1)), 'color', rgb('Purple'), 'LineWidth', 1.0 );
            plot(disp_x(stim_t)-disp_x(stim_t(1)), disp_y(stim_t)-disp_y(stim_t(1)), 'color', rgb('Crimson'), 'LineWidth', 2.0 );
        
        end
        
        xlim([ -0.1  0.3 ]);
        ylim([ -0.05 1.5 ]);
    end
    
    title( [ ac.task_str{ trial_type } ':' condition_trials_str{ cond } ]);
    fname = [ analysis_path '/sid_' num2str( sid ) '_traj_' ac.task_str{ trial_type } '_' condition_trials_str{ cond } ];
    saveas(f, [ fname '.fig'] );
    saveas(f, [ fname '.png'] );
end

end

