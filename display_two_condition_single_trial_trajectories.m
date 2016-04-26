function display_two_condition_single_trial_trajectories( sid, condition_trials, condition_trials_str, bdata_vel_time, traj, analysis_path )

ac = get_analysis_constants;
settings = sensor_settings;

prestim = settings.pre_stim;
stim    = settings.stim;

stim_t = find((bdata_vel_time > prestim) & (bdata_vel_time < (prestim+stim)));

tbt_analysis_path = [analysis_path '/sid_' num2str(sid) '_trial_by_trial_behavior/'];

if(~exist(tbt_analysis_path, 'dir'))
    mkdir( tbt_analysis_path );
end

for trial_type = 1:size( condition_trials, 1 )
    
    for cond = 1:size( condition_trials, 2 )
        
        for trial = 1:length(condition_trials{trial_type, cond})
            
            f = figure; 
            cur_trial = condition_trials{trial_type,cond}(trial);
            
            disp_x = squeeze(traj{trial_type}( cur_trial, 1, :));
            disp_y = squeeze(traj{trial_type}( cur_trial, 2, :));
            
            hold on;
            plot(disp_x, disp_y, 'color', rgb('Purple'), 'LineWidth', 2.0 );
            plot(disp_x(stim_t), disp_y(stim_t), 'color', rgb('Crimson'), 'LineWidth', 3.0 );
        
            xlim([ -0.1  0.3 ]);
            ylim([ -0.05 1.5 ]);
            title( [ ac.task_str{ trial_type } ':' condition_trials_str{ cond } ]);
            fname = [ tbt_analysis_path '/sid_' num2str(sid) '_traj_' ac.task_str{trial_type} '_' condition_trials_str{ cond } '_tid_' num2str(trial)];
            saveas(f, [ fname '.fig'] );
            saveas(f, [ fname '.png'] );
            close( f );
         end

    end    
end

end

