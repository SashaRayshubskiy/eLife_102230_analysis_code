function display_two_condition_difference_image( condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets, VPS, filename_prefix )

ac = get_analysis_constants;
settings = sensor_settings;

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

IMAGE_ROWS = 4;
IMAGE_COLS = 4;
PLANES = IMAGE_ROWS * IMAGE_COLS;

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

total_time = prestim + stim + poststim;

first_stim_t = prestim;
last_stim_t = stim + prestim;

x_size = size(avg_df_f_per_condition_per_plane, 3);
y_size = size(avg_df_f_per_condition_per_plane, 4);
nframes = size(avg_df_f_per_condition_per_plane, 5);

t = zeros(PLANES,nframes,'double');
for p=1:PLANES
    t(p,:) = (([0:nframes-1]))./VPS + frame_start_offsets(p);
end

%for trial_type = 1:size( btraces_per_condition, 2 )
for trial_type = 2
        
    f = figure('units','normalized','outerposition',[0 0 1 1]);
       
    % for p=1:PLANES
    for p = 10
        subaxis( IMAGE_ROWS+1, IMAGE_COLS, p, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN );
                        
        cur_plane_avg_df_f_cond_1 = squeeze(avg_df_f_per_condition_per_plane(1,p,:,:,:));
        cur_plane_avg_df_f_cond_1_filt_r = smoothts(reshape(cur_plane_avg_df_f_cond_1,[x_size*y_size, nframes]));
        cur_plane_avg_df_f_cond_1_filt = reshape(cur_plane_avg_df_f_cond_1_filt_r, [x_size, y_size, nframes]);

        cur_plane_avg_df_f_cond_2 = squeeze(avg_df_f_per_condition_per_plane(2,p,:,:,:));
        cur_plane_avg_df_f_cond_2_filt_r = smoothts(reshape(cur_plane_avg_df_f_cond_2,[x_size*y_size, nframes]));
        cur_plane_avg_df_f_cond_2_filt = reshape(cur_plane_avg_df_f_cond_2_filt_r, [x_size, y_size, nframes]);
        
        cur_t = t(p,:);

        % Extract frames during stim only for now.
        cur_frames = find((cur_t >= prestim) & (cur_t<=(prestim+stim)));

        avg_df_f_img_cond_1 = squeeze(mean(cur_plane_avg_df_f_cond_1_filt(:,:,cur_frames),3));
        avg_df_f_img_cond_2 = squeeze(mean(cur_plane_avg_df_f_cond_2_filt(:,:,cur_frames),3));

        diff_img =  avg_df_f_img_cond_1 - avg_df_f_img_cond_2;
        
        %TOL = 0.3;
        %diff_img( find( diff_img >  TOL ) ) = 0.0;
        %diff_img( find( diff_img <  -1.0*TOL ) ) = 0.0;        

        imagesc(diff_img);
        axis off;
        axis image;
        colormap jet;
        caxis([-0.05 0.15]);
        %caxis([-0.05 0.15]);

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

