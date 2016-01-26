function display_two_condition_trials( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel )

ac = get_analysis_constants;
settings = sensor_settings;

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

total_time = prestim + stim + poststim;
first_stim_t = prestim;
last_stim_t = stim + prestim;

for trial_type = 1:size(condition_trials,1)
    for cond_ord = 1:2        
        
        if(cond_ord == 1)
            cur_cond_symbol = '-';
        else
            cur_cond_symbol = '--';
        end
        
        subaxis(1,3, trial_type, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
        
        hold on;
        
        yaw_traces = [];
        fwd_traces = [];
        for trial_ord = 1:size(condition_trials{trial_type}, 2)           
            
            cur_trial = condition_trials{ trial_type }( cond_ord );
            yaw_trace = bdata_vel{trial_type}( cur_trial, ac.VEL_YAW, :);
            fwd_trace = bdata_vel{trial_type}( cur_trial, ac.VEL_FWD, :);

            yaw_traces(end+1,:) = yaw_trace;
            fwd_traces(end+1,:) = fwd_trace;
            
            plot( bdata_vel_time, yaw_trace, 'color', rgb('LightSalmon'), 'LineSpec', cur_cond_symbol, 'LineWidth', 0.5 );
            plot( bdata_vel_time, fwd_trace, 'color', rgb('PaleGreen'), 'LineSpec', cur_cond_symbol, 'LineWidth', 0.5 );      
        
        end
        
        % Plot average trial
        avg_yaw_trace = mean( yaw_traces );
        avg_fwd_trace = mean( fwd_traces );
        plot( bdata_vel_time, avg_yaw_trace, 'color', rgb('FireBrick'), 'LineSpec', cur_cond_symbol, 'LineWidth', 2.0 );
        plot( bdata_vel_time, avg_fwd_trace, 'color', rgb('SeaGreen'), 'LineSpec', cur_cond_symbol, 'LineWidth', 2.0 );              
        
        if( ( trial_type == 1 ) && ( cond_ord == 2 ))
            legend( [ phdl(1,1), phdl(2,1), phdl(1,2), phdl(2,2) ], ...
                ['Vel fwd - ' condition_trials_str(1)], ['Vel yaw - ' condition_trials_str(1)], ...
                ['Vel fwd - ' condition_trials_str(2)], ['Vel yaw - ' condition_trials_str(2)] );
        end
        
        yy = ylim;
        y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');
        
        xlim([0, total_time]);
        xlabel('Time (s)');
        ylabel('Velocity (au/s)');
        
        drawnow;
    end
end

end

