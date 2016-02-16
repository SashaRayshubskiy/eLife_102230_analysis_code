function display_two_condition_difference_image( down_params, ref_imgs, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets, VPS, filename_prefix )

ac = get_analysis_constants;
settings = sensor_settings;

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

PLANES = size(avg_df_f_per_condition_per_plane,3);
IMAGE_ROWS = floor(sqrt(PLANES));
IMAGE_COLS = IMAGE_ROWS;

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

total_time = prestim + stim + poststim;

first_stim_t = prestim;
last_stim_t = stim + prestim;

x_size = size(avg_df_f_per_condition_per_plane, 4);
y_size = size(avg_df_f_per_condition_per_plane, 5);
nframes = size(avg_df_f_per_condition_per_plane, 6);

t = zeros(PLANES,nframes,'double');
for p=1:PLANES
    t(p,:) = (([0:nframes-1]))./VPS + frame_start_offsets(p);
end

for trial_type = 1:size( btraces_per_condition, 2 )
        
    f = figure('units','normalized','outerposition',[0 0 1 1]);
       
     for p=1:PLANES
        subaxis( IMAGE_ROWS+1, IMAGE_COLS, p, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN );
                        
        cur_plane_avg_df_f_cond_1 = squeeze(avg_df_f_per_condition_per_plane(trial_type,1,p,:,:,:));
        cur_plane_avg_df_f_cond_1(~isfinite(cur_plane_avg_df_f_cond_1)) = 0.0;

        cur_plane_avg_df_f_cond_2 = squeeze(avg_df_f_per_condition_per_plane(trial_type,2,p,:,:,:));
        cur_plane_avg_df_f_cond_2(~isfinite(cur_plane_avg_df_f_cond_2)) = 0.0;
        
        cur_t = t(p,:);

        % Extract frames during stim only for now.
        cur_frames = find((cur_t >= prestim) & (cur_t<=(prestim+stim)));
        ref_img = squeeze(ref_imgs(p,:,:));
        ref_img_mask = get_dead_pixel_mask( ref_img );

        dx = 16;
        dy = 16;        
                
        cur_plane_avg_df_f_cond_1_down = downsample_with_mask(cur_plane_avg_df_f_cond_1, ref_img_mask, dx, dy);
        cur_plane_avg_df_f_cond_2_down = downsample_with_mask(cur_plane_avg_df_f_cond_2, ref_img_mask, dx, dy);            
        
        frames_of_interest_cond_1 = cur_plane_avg_df_f_cond_1_down(:,:,cur_frames);
        frames_of_interest_cond_2 = cur_plane_avg_df_f_cond_2_down(:,:,cur_frames);
        
        diff_img_down = trapz(frames_of_interest_cond_1,3) - trapz(frames_of_interest_cond_2,3);
        diff_img_down_scaled = expand_img(diff_img_down, dx, dy);
        ref_img_filt = ref_img.*ref_img_mask;

        C = imfuse(ref_img_filt, diff_img_down_scaled, 'ColorChannels', [2 1 0]);
        imagesc([1:y_size*down_params(2)], [1:x_size*down_params(1)], C);
        axis image;
        axis off;
        
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

