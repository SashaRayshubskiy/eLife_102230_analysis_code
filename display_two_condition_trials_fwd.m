function display_two_condition_trials_fwd( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, filepath_prefix, with_single_trials )

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

ac = get_analysis_constants;
settings = sensor_settings;

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

total_time = prestim + stim + poststim;
first_stim_t = prestim;
last_stim_t = stim + prestim;

f = figure('units','normalized','outerposition',[0 0 1 1]);

%for trial_type = 1:size(condition_trials,1)
for trial_type = [1:2]
    for cond_ord = [2 1]        

        cur_cond_symbol = '-';              

        if(cond_ord == 1)
            if( trial_type == ac.BOTH )
                cur_color_avg = rgb('DimGray');
                cur_color_single = rgb('DarkGray');                                        
            elseif( trial_type == ac.LEFT )
                cur_color_avg = rgb('FireBrick');
                cur_color_single = rgb('LightSalmon');                            
            elseif( trial_type == ac.RIGHT )
                cur_color_avg = rgb('SeaGreen');
                cur_color_single = rgb('PaleGreen');              
            end            
        else
            cur_color_avg = rgb('DimGray');
            cur_color_single = rgb('DarkGray');                                                    
        end
        
        %subaxis( 3, 1, trial_type, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN );
       subplot( 2, 1, trial_type );
        
        hold on;
        
        yaw_traces = [];
        fwd_traces = [];
        for trial_ord = 1:length(condition_trials{trial_type, cond_ord})           
            
            cur_trial = condition_trials{ trial_type, cond_ord }( trial_ord );
            yaw_trace = squeeze(bdata_vel{trial_type}( cur_trial, ac.VEL_YAW, :));
            fwd_trace = squeeze(bdata_vel{trial_type}( cur_trial, ac.VEL_FWD, :));

            yaw_traces(end+1,:) = yaw_trace;
            fwd_traces(end+1,:) = fwd_trace;
            
            if( with_single_trials )
                plot( bdata_vel_time, fwd_trace', 'color', cur_color_single, 'LineStyle', cur_cond_symbol, 'LineWidth', 0.5 );
            end
            %plot( bdata_vel_time, fwd_trace', 'color', rgb('PaleGreen'), 'LineStyle', cur_cond_symbol, 'LineWidth', 0.5 );      
        
        end
        
        % Plot average trial
        avg_yaw_trace = mean( yaw_traces );
        avg_fwd_trace = mean( fwd_traces );
        phdl(cond_ord,1) = plot( bdata_vel_time, avg_fwd_trace, 'color', rgb('SeaGreen'), 'LineStyle', cur_cond_symbol, 'LineWidth', 2.0 );              
        %phdl(cond_ord,1) = plot( bdata_vel_time, avg_yaw_trace, 'color', cur_color_avg, 'LineStyle', cur_cond_symbol, 'LineWidth', 2.0 );
        
        traces_cnt(cond_ord) = size(yaw_traces,1);
        
        if( cond_ord == 1 )
            ll = legend([ phdl(1,1), phdl(2,1) ], ...
                ['Vel fwd - ' condition_trials_str{ 1 } '( ' num2str( traces_cnt( 1 ) ) ' )'], ...
                ['Vel fwd - ' condition_trials_str{ 2 } '( ' num2str( traces_cnt( 2 ) ) ' )'] );
            set(ll, 'Interpreter', 'none');
        end
        
        ylim([-0.5 0.5]);
        yy = ylim;
        y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');
        
        xlim([0, total_time]);
        if( trial_type == 3 )
            xlabel('Time (s)');
        else
            set(gca(), 'XTickLabels', '');
        end
        title([ac.task_str(trial_type)]);
        
        ylabel('Fwd velocity (au/s)');
        
        drawnow;
    end
end

if( with_single_trials )
    saveas(f, [ filepath_prefix '_single_t.fig']);
    saveas(f, [ filepath_prefix '_single_t.png']);
else
    saveas(f, [ filepath_prefix '.fig']);
    saveas(f, [ filepath_prefix '.png']);
end

end

