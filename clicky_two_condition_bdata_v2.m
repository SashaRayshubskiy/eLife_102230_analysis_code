function clicky_two_condition_bdata_v2( ref_img, PLANE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets, VPS, filename_prefix )
% This plots both left and right time courses

ac = get_analysis_constants;
settings = sensor_settings;

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

total_time = prestim + stim + poststim;

first_stim_t = prestim;
last_stim_t = stim + prestim;

x_size = size(avg_df_f_per_condition_per_plane, 4);
y_size = size(avg_df_f_per_condition_per_plane, 5);
nframes = size(avg_df_f_per_condition_per_plane, 6);

t = zeros(1,nframes,'double');
t = (([0:nframes-1]))./VPS + frame_start_offsets(PLANE_OF_INTEREST);

npts = 1;
colorindex = 0;

order    = ac.order;
nroi = 1;
intens = [];
[x, y] = meshgrid(1:y_size, 1:x_size);
baseline_start = 0;
baseline_end = 2.8;

p = PLANE_OF_INTEREST;

f1 = figure();
               
ax1 = subplot(2,2,1);
ref_img_mask = get_dead_pixel_mask(ref_img);

[xsize, ysize] = size(ref_img);
%imagesc(imresize(ref_img, [xsize 2*ysize]));
%imagesc( ref_img.*ref_img_mask );
imagesc( ref_img );
colormap(ax1, 'gray');
axis image;
caxis([0 3500]);
% title([ac.task_str{trial_type}]);

% Get ROI
clicky_plane = 1;
clear roi_points;
while(npts > 0)
    
    subplot(2,2,clicky_plane);
    % subplot(1,3,1)
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
        break
    end
        
    subplot(2,2,clicky_plane);
    %draw the bounding polygons and label them
    currcolor    = order(1+mod(colorindex,size(order,1)),:);
    hold on;
    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
            
    roi_points{nroi} = [xv, yv];
    nroi = nroi + 1;
    colorindex = colorindex +1;
end

for trial_type = 1:2
    
    cur_plane_avg_df_f_cond_1 = squeeze(avg_df_f_per_condition_per_plane(trial_type,1,p,:,:,:));
    cur_plane_avg_df_f_cond_1(~isfinite(cur_plane_avg_df_f_cond_1)) = 0.0;
    
    cur_plane_avg_df_f_cond_2 = squeeze(avg_df_f_per_condition_per_plane(trial_type,2,p,:,:,:));
    cur_plane_avg_df_f_cond_2(~isfinite(cur_plane_avg_df_f_cond_2)) = 0.0;
    
    cur_t = t;
            
    a_data_1 = cur_plane_avg_df_f_cond_1;
    a_data_2 = cur_plane_avg_df_f_cond_2;
    
    plt_cond_1 = [];
    plt_cond_2 = [];
    
    colorindex = 0;
    
    t_test_per_roi = [];
    
    for r = 1:length(roi_points)
    
        xv = roi_points{ r }(:,1);
        yv = roi_points{ r }(:,2);
        
        inpoly = inpolygon(x,y,xv,yv);
        
        subplot(2,2,trial_type*2);
        %draw the bounding polygons and label them
        currcolor    = order(1+mod(colorindex,size(order,1)),:);
        hold on;
        plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
        text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
                
        itrace_1 = squeeze(sum(sum(double(a_data_1).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
        itrace_2 = squeeze(sum(sum(double(a_data_2).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
        
        [t_test, pval] = ttest( itrace_1, itrace_2 );
        
        t_test_per_roi( r )  = pval;
        
        hold on;
        cur_t = t;
        plt_cond_1(end+1) = plot( cur_t, itrace_1, 'Color', currcolor, 'LineWidth', 2);
        plt_cond_2(end+1) = plot( cur_t, itrace_2, 'Color', currcolor, 'LineWidth', 2, 'LineStyle', '--');
        
        xlim([0 max(cur_t)]);
        ylim([-0.1 0.4]);
        xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('dF/F');
        set(gca, 'FontSize', 14 );
        set(gca, 'FontWeight', 'bold');
        
        colorindex = colorindex+1;        
    end
    
    t_test_str = '';
    for r = 1:length( t_test_per_roi )
        t_test_str = strcat(t_test_str, ['_r_' num2str(r) '__p_' num2str(sprintf('%03f',t_test_per_roi(r)))]);
    end
    
    ax1 = subplot(2,2,trial_type*2); % plot the trace    
    
    yy = ylim;
    y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
    hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
    set(gca,'children',circshift(get(gca,'children'),-1));
    set(hh, 'EdgeColor', 'None');
    tt = title([ac.task_str{trial_type} '_' t_test_str]);
    set( tt, 'Interpreter', 'none' );
    
    cond_1_num_trials = size( btraces_per_condition{ 1, trial_type }( :, ac.VEL_YAW, : ), 1 );
    cond_2_num_trials = size( btraces_per_condition{ 2, trial_type }( :, ac.VEL_YAW, : ), 1 );
    
    ll = legend( [ plt_cond_1(1), plt_cond_2(1) ], ...
        [ condition_trials_str{ 1 } '(' num2str( cond_1_num_trials ) ')'], ...
        [ condition_trials_str{ 2 } '(' num2str( cond_2_num_trials ) ')'], 'Location', 'southeast');
    set(ll, 'Interpreter', 'none');
    
    drawnow;
end

saveas(f1, [ filename_prefix '_' ac.task_str{trial_type} '_rois.fig']);
saveas(f1, [ filename_prefix '_' ac.task_str{trial_type} '_rois.png']);

end

