function display_two_behavioral_condition_traces( condition_trials_str, btraces_per_condition, ctraces_in_roi_per_condition, bdata_vel_time, frame_start_offsets, VPS, filename_prefix )

ac = get_analysis_constants;
settings = sensor_settings;
order = ac.order;

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;


prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

total_time = prestim + stim + poststim;

first_stim_t = prestim;
last_stim_t = stim + prestim;
nframes = size(ctraces_in_roi_per_condition{1,1,1}, 3);

PLANES = size(ctraces_in_roi_per_condition,3);
IMAGE_ROWS = floor(sqrt(PLANES));
IMAGE_COLS = IMAGE_ROWS;

t = zeros(PLANES,nframes,'double');
for p=1:PLANES
    t(p,:) = (([0:nframes-1]))./VPS + frame_start_offsets(p);
end

for trial_type = 1:size( btraces_per_condition, 2 )
        
    f = figure('units','normalized','outerposition',[0 0 1 1]);
        
    for cond_ord = 1:size( btraces_per_condition, 1 )
        
        if(cond_ord == 1)
            cur_cond_symbol = '-';
        else
            cur_cond_symbol = '--';
        end
        
        for p=1:PLANES
            subaxis( IMAGE_ROWS+1, IMAGE_COLS, p, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN );
            
            colorindex = 0;
            
            cur_plane_data = ctraces_in_roi_per_condition{ cond_ord, trial_type, p };
            
            for roi_id = 1:size(cur_plane_data, 2);
                hold on;
                currcolor = order(1+mod(colorindex,size(order,1)),:);

                avg_trace = mean(squeeze(cur_plane_data(:,roi_id,:)));
                
                plot( squeeze(t(p,:)), avg_trace, 'color', currcolor, 'LineStyle', cur_cond_symbol );
                colorindex = colorindex + 1;
            end
            
            if(p>=1 & p<=8 )
                ylim([ -0.2 1.2 ]);
            elseif (p>=9 & p<=12 )
                ylim([ -0.2 0.8 ]);
            elseif (p>=13 & p<=16 )
                ylim([ -0.2 0.4 ]);
            end
            
            yy = ylim;
            y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
            hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
            set(gca,'children',circshift(get(gca,'children'),-1));
            set(hh, 'EdgeColor', 'None');
            
            xlim([0, total_time]);
            if(mod((p-1),IMAGE_COLS) == 0 )
                ylabel('dF/F');
            else
                set(gca, 'YTickLabel', '');
            end
            
            set(gca, 'XTickLabel', '');
            
            if( p == 2 )
                tt = title(ac.task_str(trial_type));
                set(tt, 'Interpreter', 'none')
            end
            drawnow;
        end
    
        for c = 1:IMAGE_COLS
            
            % Axis for behavioral data
            subaxis(IMAGE_ROWS+1, IMAGE_COLS, PLANES + c, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
            
            hold on;
            %avg_trace_fwd = mean(squeeze(btraces_per_condition{ cond_ord, trial_type }( :, ac.VEL_FWD, : )));
            avg_trace_yaw = mean(squeeze(btraces_per_condition{ cond_ord, trial_type }( :, ac.VEL_YAW, : )));
            
            %phdl(cond_ord, 1) = plot( bdata_vel_time, avg_trace_fwd, 'color', rgb('FireBrick'), 'LineStyle', cur_cond_symbol );
            phdl(cond_ord) = plot( bdata_vel_time, avg_trace_yaw, 'color', rgb('SeaGreen'), 'LineStyle', cur_cond_symbol );
            
            cond_num_trials(cond_ord) = size( btraces_per_condition{ cond_ord, trial_type }( :, ac.VEL_YAW, : ), 1 );
            
            if( ( c == 1 ) & ( cond_ord == 2 ))
                ll = legend( [ phdl(1), phdl(2) ], ... 
                        [ condition_trials_str{ 1 } '(' num2str( cond_num_trials( 1 ) ) ')'], ...
                        [ condition_trials_str{ 2 } '(' num2str( cond_num_trials( 2 ) ) ')'] );
                set(ll, 'Interpreter', 'none')
            end
            
            yy = ylim;
            y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
            hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
            set(gca,'children',circshift(get(gca,'children'),-1));
            set(hh, 'EdgeColor', 'None');
            
            xlim([0, total_time]);
            xlabel('Time (s)');
            if( c == 1 )
                ylabel('Yaw velocity (au/s)');
            else
                set(gca, 'YTickLabel', '');
            end
            
            drawnow;
        end
    end
        
    saveas(f, [ filename_prefix '_' ac.task_str{trial_type} '.fig']);
    saveas(f, [ filename_prefix '_' ac.task_str{trial_type} '.png']);
    % close(f);
end

end

