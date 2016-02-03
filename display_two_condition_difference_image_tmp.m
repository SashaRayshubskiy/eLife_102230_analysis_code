function display_two_condition_difference_image_tmp( condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets, VPS, filename_prefix )

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
for trial_type = 2
        
    f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    f2 = figure('units','normalized','outerposition',[0 0 1 1]);
       
    % for p=1:PLANES
    for p = 10
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
        cur_frames = find((cur_t >= prestim) & (cur_t<=(prestim+stim)));

        avg_df_f_img_cond_1 = squeeze(max(cur_plane_avg_df_f_cond_1(:,:,cur_frames),[],3));
        avg_df_f_img_cond_2 = squeeze(max(cur_plane_avg_df_f_cond_2(:,:,cur_frames),[],3));

        diff_img =  avg_df_f_img_cond_1 - avg_df_f_img_cond_2;
        
        %TOL = 0.3;
        %diff_img( find( diff_img >  TOL ) ) = 0.0;
        %diff_img( find( diff_img <  -1.0*TOL ) ) = 0.0;        

        figure(f1);

        subplot(3,1,1); 
        imagesc(avg_df_f_img_cond_1);
        axis image;
        colormap jet;
        caxis([-0 0.5]);

        subplot(3,1,2); 
        imagesc(avg_df_f_img_cond_2);
        axis image;
        colormap jet;
        caxis([-0 0.5]);       

        subplot(3,1,3);
        save('/tmp/diff_img.mat', 'diff_img');
        diff_img_filt = filter_dithered_image( diff_img );
        imagesc( diff_img_filt );
        axis image;
        colormap jet;
        caxis([-0.1 0.8]);       
        
        a_data_1 = cur_plane_avg_df_f_cond_1;
        a_data_2 = cur_plane_avg_df_f_cond_2;

        while(npts > 0)
            
            figure(f1)
            subplot(3,1,3);
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
            
            %bline_s = floor(baseline_start*VPS);
            bline_s = 1;
            bline_e = floor(baseline_end*VPS);
            
            itrace_1 = squeeze(sum(sum(double(a_data_1).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
            itrace_2 = squeeze(sum(sum(double(a_data_2).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
            
            figure(f2);
            %ax1 = subplot(1,3,2:3); % plot the trace
            subplot(1,1,1)
            hold on;
            cur_t = squeeze(t(p,:));
            plot( cur_t, itrace_1, 'Color', currcolor, 'LineWidth', 2);
            plot( cur_t, itrace_2, 'Color', currcolor, 'LineWidth', 2, 'LineStyle', '--');
            
            xlim([0 max(cur_t)]);
            %ylim([-0.2 0.75]);
            xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
            ylabel('dF/F');
            set(gca, 'FontSize', 14 );
            set(gca, 'FontWeight', 'bold');
            
            colorindex = colorindex+1;
          
            roi_points{nroi} = [xv, yv];
            nroi = nroi + 1;
        end
        
        if( p == 2 )
            tt = title(ac.task_str(trial_type));
            set(tt, 'Interpreter', 'none');
        end

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

    %saveas(f, [ filename_prefix '_' ac.task_str{trial_type} '.fig']);
    %saveas(f, [ filename_prefix '_' ac.task_str{trial_type} '.png']);
    % close(f);
end

end

