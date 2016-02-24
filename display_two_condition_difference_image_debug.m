function display_two_condition_difference_image_debug( ref_img, PLANE_OF_INTEREST, TRIAL_TYPE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets, VPS, filename_prefix )

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

x_size = size(avg_df_f_per_condition_per_plane, 4);
y_size = size(avg_df_f_per_condition_per_plane, 5);
nframes = size(avg_df_f_per_condition_per_plane, 6);

t = zeros(PLANES,nframes,'double');
for p=1:PLANES
    t(p,:) = (([0:nframes-1]))./VPS + frame_start_offsets(p);
end

npts = 1;
colorindex = 0;

order    = [ rgb('Blue'); rgb('Green'); rgb('Red'); rgb('Black'); rgb('Purple'); rgb('Brown'); rgb('Indigo'); rgb('DarkRed') ];
nroi = 1;
intens = [];
[x, y] = meshgrid(1:y_size, 1:x_size);
baseline_start = 0;
baseline_end = 2.8;

%for trial_type = 1:size( btraces_per_condition, 2 )
for trial_type = TRIAL_TYPE_OF_INTEREST
        
%    f1 = figure('units','normalized','outerposition',[0 0 1 1]);
%    f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    f1 = figure();
    f2 = figure();
       
    % for p=1:PLANES
    for p = PLANE_OF_INTEREST
        %subaxis( IMAGE_ROWS+1, IMAGE_COLS, p, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN );
        %subplot(1,3,1); 
               
        cur_plane_avg_df_f_cond_1 = squeeze(avg_df_f_per_condition_per_plane(trial_type,1,p,:,:,:));
        cur_plane_avg_df_f_cond_1(~isfinite(cur_plane_avg_df_f_cond_1)) = 0.0;
        %cur_plane_avg_df_f_cond_1_filt_r = smoothts(reshape(cur_plane_avg_df_f_cond_1,[x_size*y_size, nframes]));
        %cur_plane_avg_df_f_cond_1_filt = reshape(cur_plane_avg_df_f_cond_1_filt_r, [x_size, y_size, nframes]);

        cur_plane_avg_df_f_cond_2 = squeeze(avg_df_f_per_condition_per_plane(trial_type,2,p,:,:,:));
        cur_plane_avg_df_f_cond_2(~isfinite(cur_plane_avg_df_f_cond_2)) = 0.0;
        %cur_plane_avg_df_f_cond_2_filt_r = smoothts(reshape(cur_plane_avg_df_f_cond_2,[x_size*y_size, nframes]));
        %cur_plane_avg_df_f_cond_2_filt = reshape(cur_plane_avg_df_f_cond_2_filt_r, [x_size, y_size, nframes]);
        
        cur_t = squeeze(t(p,:));

        % Extract frames during stim only for now.
        PAD = 1.0;
        cur_frames = find((cur_t >= (prestim-PAD)) & (cur_t<=(prestim+stim+PAD)));

        avg_df_f_img_cond_1 = squeeze(mean(cur_plane_avg_df_f_cond_1(:,:,cur_frames),3));
        avg_df_f_img_cond_2 = squeeze(mean(cur_plane_avg_df_f_cond_2(:,:,cur_frames),3));

        figure(f1);

        ax1 = subplot(2,2,1); 
        ref_img_mask = get_dead_pixel_mask(ref_img);
        save('/tmp/ref_img.mat', 'ref_img');

        
        [xsize, ysize] = size(ref_img);
        %imagesc(imresize(ref_img, [xsize 2*ysize]));
        imagesc( ref_img.*ref_img_mask );
        colormap(ax1, 'gray');
        axis image;
        caxis([0 3000]);
        title([ac.task_str{trial_type}]);

        ax2 = subplot(2,2,2); 
        imagesc(avg_df_f_img_cond_1.*ref_img_mask);
        axis image;
        colormap(ax2, 'jet');
        caxis([-0 0.5]);
        tt = title(['Condition 1: ' condition_trials_str{1}]);
        set(tt, 'Interpreter', 'none');

        ax3 = subplot(2,2,4); 
        imagesc(avg_df_f_img_cond_2.*ref_img_mask);
        axis image;
        colormap(ax3, jet);
        caxis([-0 0.5]);       
        tt = title(['Condition 2: ' condition_trials_str{2}]);
        set(tt, 'Interpreter', 'none');

        ax4 = subplot(2,2,3);
       
        dx = 16;
        dy = 16;        
        
        %cur_plane_avg_df_f_cond_1_down = squeeze(mean(mean(reshape(cur_plane_avg_df_f_cond_1, [dx, xsize/dx, dy, ysize/dy, nframes ]),3),1));        
        %cur_plane_avg_df_f_cond_2_down = squeeze(mean(mean(reshape(cur_plane_avg_df_f_cond_2, [dx, xsize/dx, dy, ysize/dy, nframes ]),3),1));        
        
        UPSAMPLE_FACTOR = 100;

        cur_t_up = squeeze(t(p,:));

        cur_plane_avg_df_f_cond_1_down = downsample_with_mask(cur_plane_avg_df_f_cond_1, ref_img_mask, dx, dy);
        cur_plane_avg_df_f_cond_2_down = downsample_with_mask(cur_plane_avg_df_f_cond_2, ref_img_mask, dx, dy);

        cur_plane_avg_df_f_cond_1_up = upsample_in_t( cur_plane_avg_df_f_cond_1_down, UPSAMPLE_FACTOR );
        cur_plane_avg_df_f_cond_2_up = upsample_in_t( cur_plane_avg_df_f_cond_2_down, UPSAMPLE_FACTOR );        
                    
        VPS_UPSAMPLE = VPS * UPSAMPLE_FACTOR;
        cutoff_freq = 0.5;
        cur_plane_avg_df_f_cond_1_up_filt = fft_filter_3D( cur_plane_avg_df_f_cond_1_up, cutoff_freq, VPS_UPSAMPLE );
        cur_plane_avg_df_f_cond_2_up_filt = fft_filter_3D( cur_plane_avg_df_f_cond_2_up, cutoff_freq, VPS_UPSAMPLE );                
        
        %diff_frames = find((cur_t >= (prestim)) & (cur_t<=(prestim+stim)));
        diff_frames = [1:length(cur_plane_avg_df_f_cond_1_up_filt)];
        frames_of_interest_cond_1 = cur_plane_avg_df_f_cond_1_up_filt( :, :, diff_frames );
        frames_of_interest_cond_2 = cur_plane_avg_df_f_cond_2_up_filt( :, :, diff_frames );
        
        %diff_img_down = trapz(frames_of_interest_cond_1,3) - trapz(frames_of_interest_cond_2,3);
        dx_size = size(frames_of_interest_cond_1,1);
        dy_size = size(frames_of_interest_cond_1,2);
        
        diff_img_down = ones(dx_size, dy_size);
        kept_trials = zeros(dx_size * dy_size, 3);
        
        kept_trials_index = 1;
        
        if 0
        P_VALUE_THRESHOLD = 0.01;
        for ii = 1:dx_size
            for jj = 1:dy_size
                cur_p = signrank(squeeze(frames_of_interest_cond_1(ii,jj,:)), squeeze(frames_of_interest_cond_2(ii,jj,:)));
                if(cur_p < P_VALUE_THRESHOLD )
                    diff_img_down(ii,jj) = cur_p;
                    
                    % add debugging info
                    kept_trials( kept_trials_index, : ) = [ cur_p, ii, jj ];
                    
                    kept_trials_index = kept_trials_index + 1;
                end
            end
        end
        
        kept_trials_p_value_sorted = sortrows( kept_trials(1:kept_trials_index-1,:), [1 2 3]);
        
        figure;
        active_roi_cnt = length(kept_trials_p_value_sorted);
        NUM_ROI_PER_AXES = 2;
        num_axes = floor(active_roi_cnt/NUM_ROI_PER_AXES);
        AXES_ROWS = floor(sqrt(num_axes));
        AXES_COLS = ceil(num_axes/AXES_ROWS);

        cur_t = squeeze(t(p,:));
        roi_idx = 1;
        for a = 1:num_axes
            subaxis( AXES_ROWS, AXES_COLS, a, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN );
            for kk = 0:NUM_ROI_PER_AXES-1

                currcolor    = order(1+mod(kk,size(order,1)),:);

                pv = kept_trials_p_value_sorted( roi_idx, 1 );
                xx = kept_trials_p_value_sorted( roi_idx, 2 );
                yy = kept_trials_p_value_sorted( roi_idx, 3 );
                
                itrace_1 = squeeze(cur_plane_avg_df_f_cond_1_down(xx,yy,:));
                itrace_2 = squeeze(cur_plane_avg_df_f_cond_2_down(xx,yy,:));

                hold on;
                plt_1(kk+1) = plot( cur_t, itrace_1, 'Color', currcolor, 'LineWidth', 2);
                plt_2(kk+1) = pv;
                plot( cur_t, itrace_2, 'Color', currcolor, 'LineWidth', 2, 'LineStyle', '--');
                
                roi_idx = roi_idx + 1;
            end
            
            legend([plt_1(1), plt_1(2)], ['pv: ' num2str(plt_2(1), '%10.1e')], ['pv: ' num2str(plt_2(2), '%10.1e')], 'location', 'southeast');            
            
            xlim([0 6.5]);
            ylim([-0.5 0.5]);
            yy = ylim;
            y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
            hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
            set(gca,'children',circshift(get(gca,'children'),-1));
            set(hh, 'EdgeColor', 'None');
        end
        end
        
        figure(f1);
        imagesc( [1:ysize], [1:xsize], diff_img_down );
        xlim([1 ysize]);
        ylim([1 xsize]);
        axis image;
        colormap(ax4, 'jet');
        caxis(ax4,[0.0 0.05]);       
        %colorbar;
        title('Diff img');
        
        return;
        
        a_data_1 = cur_plane_avg_df_f_cond_1;
        a_data_2 = cur_plane_avg_df_f_cond_2;

        plt_cond_1 = [];
        plt_cond_2 = [];
        
        clicky_plane = 3;
        while(npts > 0)
            
            figure(f1)
            subplot(2,2,clicky_plane);
            % subplot(1,3,1)
            [xv, yv] = (getline(gca, 'closed'));
            if size(xv,1) < 3  % exit loop if only a line is drawn
                break
            end
            
            inpoly = inpolygon(x,y,xv,yv);
                                   
            %draw the bounding polygons and label them
            currcolor    = order(1+mod(colorindex,size(order,1)),:);
            hold on;
            plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
            text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
            
            cur_t = squeeze(t(p,:));
            
            if 0
                bline_s = 1;
                bline_e = floor(baseline_end*VPS);
                
                itrace_1 = squeeze(sum(sum(double(a_data_1).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
                itrace_2 = squeeze(sum(sum(double(a_data_2).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
            end
            
            xx = xv(1);
            xx_down = ceil( abs(xx)/dx );
            yy = yv(1);
            yy_down = ceil( abs(yy)/dy );
            
            itrace_1 = squeeze(cur_plane_avg_df_f_cond_1_down(yy_down,xx_down,:));
            itrace_2 = squeeze(cur_plane_avg_df_f_cond_2_down(yy_down,xx_down,:));
            
            figure(f2);
            subplot(1,1,1)
            hold on;
            plt_cond_1(end+1) = plot( cur_t, itrace_1, 'Color', currcolor, 'LineWidth', 2);
            plt_cond_2(end+1) = plot( cur_t, itrace_2, 'Color', currcolor, 'LineWidth', 2, 'LineStyle', '--');
            
            xlim([0 max(cur_t)]);
            xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
            ylabel('dF/F');
            set(gca, 'FontSize', 14 );
            set(gca, 'FontWeight', 'bold');
            
            colorindex = colorindex+1;
          
            roi_points{nroi} = [xv, yv];
            nroi = nroi + 1;
        end                
        
        figure( f2 );
        %ax1 = subplot(1,3,2:3); % plot the trace
        subplot(1,1,1)

        yy = ylim;
        y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');

        cond_1_num_trials = size( btraces_per_condition{ 1, trial_type }( :, ac.VEL_YAW, : ), 1 );
        cond_2_num_trials = size( btraces_per_condition{ 2, trial_type }( :, ac.VEL_YAW, : ), 1 );
        
        ll = legend( [ plt_cond_1(1), plt_cond_2(1) ], ...
                     [ condition_trials_str{ 1 } '(' num2str( cond_1_num_trials ) ')'], ...
                     [ condition_trials_str{ 2 } '(' num2str( cond_2_num_trials ) ')'], 'Location', 'southeast');
        set(ll, 'Interpreter', 'none');
        
        drawnow;
    end
    
if 0
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
end

    saveas(f1, [ filename_prefix '_' ac.task_str{trial_type} '_rois.fig']);
    saveas(f1, [ filename_prefix '_' ac.task_str{trial_type} '_rois.png']);
    saveas(f2, [ filename_prefix '_' ac.task_str{trial_type} '_tc.fig']);
    saveas(f2, [ filename_prefix '_' ac.task_str{trial_type} '_tc.png']);
    % close(f);
end

end

