function display_two_condition_trials( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel )

ac = get_analysis_constants;
settings = sensor_settings;
order = ac.order;

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

for trial_type = 1:size(condition_trials,1)
    for cond_ord = 1:2        
        
        if(cond_ord == 1)
            cur_cond_symbol = '-';
        else
            cur_cond_symbol = '--';
        end
        
        subaxis(1,3, trial_type, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
        
        hold on;
        
        for trial_ord = 1:size(
        
        avg_trace_fwd = mean(squeeze(btraces_per_condition( cond_ord, trial_type, :, ac.VEL_FWD, : )));
        avg_trace_yaw = mean(squeeze(btraces_per_condition( cond_ord, trial_type, :, ac.VEL_YAW, : )));
        
        phdl(cond_ord, 1) = plot( bdata_vel_time, avg_trace_fwd, 'color', rgb('FireBrick'), 'LineSpec', cur_cond_symbol );
        phdl(cond_ord, 2) = plot( bdata_vel_time, avg_trace_yaw, 'color', rgb('SeaGreen'), 'LineSpec', cur_cond_symbol );
        
        if( ( c == 1 ) & ( cond_ord == 2 ))
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

