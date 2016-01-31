function display_two_behavioral_condition_diff_traces( condition_trials_str, btraces_per_condition, ctraces_in_roi_per_condition, bdata_vel_time, VPS, filename_prefix )

ac = get_analysis_constants;
settings = sensor_settings;
order = ac.order;

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

IMAGE_ROWS = 4;
IMAGE_COLS = 4;
PLANES = IMAGE_ROWS * IMAGE_COLS;

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

base_begin = 1;
base_end = floor(prestim*VPS);

total_time = prestim + stim + poststim;

first_stim_t = prestim;
last_stim_t = stim + prestim;
nframes = size(ctraces_in_roi_per_condition{1,1,1}, 3);

t = [1:nframes]./VPS;

for trial_type = 1:size( btraces_per_condition, 2 )
        
    f = figure('units','normalized','outerposition',[0 0 1 1]);
       
    for p=1:PLANES
        subaxis( IMAGE_ROWS+1, IMAGE_COLS, p, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN );
        
        colorindex = 0;
        
        cur_plane_data_cond_1 = ctraces_in_roi_per_condition{ 1, trial_type, p };
        cur_plane_data_cond_2 = ctraces_in_roi_per_condition{ 2, trial_type, p };
        
        for roi_id = 1:size(cur_plane_data_cond_1, 2);
            hold on;
            currcolor = order(1+mod(colorindex,size(order,1)),:);
            
            avg_trace_cond_1 = mean(squeeze(cur_plane_data_cond_1(:,roi_id,:)));
            avg_trace_cond_2 = mean(squeeze(cur_plane_data_cond_2(:,roi_id,:)));
            
            plot( t, (avg_trace_cond_1-avg_trace_cond_2), 'color', currcolor );
            colorindex = colorindex + 1;
        end
        
        min_ylim = -0.3;
        max_ylim = 0.35;
        
        if(p>=1 & p<=8 )
            ylim([ min_ylim max_ylim ]);
        elseif (p>=9 & p<=12 )
            ylim([ min_ylim max_ylim ]);
        elseif (p>=13 & p<=16 )
            ylim([ min_ylim max_ylim ]);
        end            
        
        yy = ylim;
        y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');
        
        xlim([0, total_time]);
        if(mod((p-1),4) == 0 )
            ylabel('diff dF/F');
        else
            set(gca, 'YTickLabel', '');
        end
        
        set(gca, 'XTickLabel', '');
        
        if( p == 2 )
            tt = title(ac.task_str(trial_type));
            set(tt, 'Interpreter', 'none');
        end
        drawnow;
    end
    
    for c = 1:IMAGE_COLS
        
        % Axis for behavioral data
        subaxis(IMAGE_ROWS+1, IMAGE_COLS, PLANES + c, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
        
        hold on;

        avg_trace_yaw_cond_1 = mean(squeeze(btraces_per_condition{ 1, trial_type }( :, ac.VEL_YAW, : )));
        avg_trace_yaw_cond_2 = mean(squeeze(btraces_per_condition{ 2, trial_type }( :, ac.VEL_YAW, : )));
        
        %phdl(cond_ord, 1) = plot( bdata_vel_time, avg_trace_fwd, 'color', rgb('FireBrick'), 'LineStyle', cur_cond_symbol );
        phdl(1) = plot( bdata_vel_time, avg_trace_yaw_cond_1, 'color', rgb('SeaGreen'), 'LineStyle', '-' );
        phdl(2) = plot( bdata_vel_time, avg_trace_yaw_cond_2, 'color', rgb('SeaGreen'), 'LineStyle', '--' );
        
        cond_num_trials( 1 ) = size( btraces_per_condition{ 1, trial_type }( :, ac.VEL_YAW, : ), 1 );
        cond_num_trials( 2 ) = size( btraces_per_condition{ 2, trial_type }( :, ac.VEL_YAW, : ), 1 );
        
        if( c == 1 )
            ll = legend( [ phdl(1), phdl(2) ], ...
                [ condition_trials_str{ 1 } '(' num2str( cond_num_trials( 1 ) ) ')'], ...
                [ condition_trials_str{ 2 } '(' num2str( cond_num_trials( 2 ) ) ')'] );
            set(ll, 'Interpreter', 'none');
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

    saveas(f, [ filename_prefix '_' ac.task_str{trial_type} '.fig']);
    saveas(f, [ filename_prefix '_' ac.task_str{trial_type} '.png']);
    % close(f);
end

end

