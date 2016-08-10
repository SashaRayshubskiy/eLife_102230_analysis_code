function clicky_two_condition_bdata_8roi_LR( ref_img, PLANE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets, VPS, filename_prefix )

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

colorindex = 0;

order    = ac.order;
intens = [];
[x, y] = meshgrid(1:y_size, 1:x_size);
baseline_start = 0;
baseline_end = 2.8;

p = PLANE_OF_INTEREST;

cur_t = t;
% Extract frames during stim only for now.
cur_frames = find((cur_t >= prestim) & (cur_t<=(prestim+stim)));

f1 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
f2 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
               
l_cur_plane_avg_df_f_cond_1 = squeeze(avg_df_f_per_condition_per_plane(ac.LEFT,1,p,:,:,:));
l_cur_plane_avg_df_f_cond_1(~isfinite(l_cur_plane_avg_df_f_cond_1)) = 0.0;

l_cur_plane_avg_df_f_cond_2 = squeeze(avg_df_f_per_condition_per_plane(ac.LEFT,2,p,:,:,:));
l_cur_plane_avg_df_f_cond_2(~isfinite(l_cur_plane_avg_df_f_cond_2)) = 0.0;

r_cur_plane_avg_df_f_cond_1 = squeeze(avg_df_f_per_condition_per_plane(ac.RIGHT,1,p,:,:,:));
r_cur_plane_avg_df_f_cond_1(~isfinite(r_cur_plane_avg_df_f_cond_1)) = 0.0;

r_cur_plane_avg_df_f_cond_2 = squeeze(avg_df_f_per_condition_per_plane(ac.RIGHT,2,p,:,:,:));
r_cur_plane_avg_df_f_cond_2(~isfinite(r_cur_plane_avg_df_f_cond_2)) = 0.0;
        
l_avg_df_f_img_cond_1 = squeeze(mean(l_cur_plane_avg_df_f_cond_1(:,:,cur_frames),3));
l_avg_df_f_img_cond_2 = squeeze(mean(l_cur_plane_avg_df_f_cond_2(:,:,cur_frames),3));

figure(f1);
ax1 = subplot(5,4,1);
ref_img_mask = get_dead_pixel_mask(ref_img);
[xsize, ysize] = size(ref_img);
%imagesc(imresize(ref_img, [xsize 2*ysize]));
%imagesc( ref_img.*ref_img_mask );
imagesc( ref_img );
colormap(ax1, 'gray');
axis image;
axis off;
caxis([0 900]);
title([ac.task_str{ac.LEFT} ' plane: ' num2str(p)]);

figure(f2);
ax2 = subplot(1,1,1);
imagesc(l_avg_df_f_img_cond_1.*ref_img_mask);
axis image;
axis off;
colormap(ax2, 'jet');
caxis([-0 0.5]);
tt = title(['Condition 1: ' condition_trials_str{1}]);
set(tt, 'Interpreter', 'none');

figure(f1);
ax2 = subplot(5,4,2);
imagesc(l_avg_df_f_img_cond_1.*ref_img_mask);
axis image;
axis off;
colormap(ax2, 'jet');
caxis([-0 0.5]);
tt = title(['Condition 1: ' condition_trials_str{1}]);
set(tt, 'Interpreter', 'none');

ax3 = subplot(5,4,3);
imagesc(l_avg_df_f_img_cond_2.*ref_img_mask);
axis image;
axis off;
colormap(ax3, jet);
caxis([-0 0.5]);
tt = title(['Condition 2: ' condition_trials_str{2}]);
set(tt, 'Interpreter', 'none');
    
l_a_data_1 = l_cur_plane_avg_df_f_cond_1;
l_a_data_2 = l_cur_plane_avg_df_f_cond_2;

r_a_data_1 = r_cur_plane_avg_df_f_cond_1;
r_a_data_2 = r_cur_plane_avg_df_f_cond_2;

plt_cond_1 = [];
plt_cond_2 = [];

roi_map = { [5 6], [9 10], [13 14], [17 18], [7 8], [11 12], [15 16], [19 20] };

clicky_plane = 2;

for r = 1:8
    
    %subplot(5,4,clicky_plane);
    % subplot(1,3,1)
    figure(f2);    
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
        break
    end
    inpoly = inpolygon(x,y,xv,yv);
        
    %subplot(5,4,2);
    %draw the bounding polygons and label them
    currcolor    = order(1+mod(colorindex,size(order,1)),:);
    hold on;
    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    
    figure(f1);
    subplot(5,4,2);
    %draw the bounding polygons and label them
    currcolor    = order(1+mod(colorindex,size(order,1)),:);
    hold on;
    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    
    pause(1.0);
    
    %bline_s = floor(baseline_start*VPS);
    bline_s = 1;
    bline_e = floor(baseline_end*VPS);
    
    l_itrace_1 = squeeze(sum(sum(double(l_a_data_1).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
    l_itrace_2 = squeeze(sum(sum(double(l_a_data_2).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));

    r_itrace_1 = squeeze(sum(sum(double(r_a_data_1).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
    r_itrace_2 = squeeze(sum(sum(double(r_a_data_2).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
        
    ax1 = subplot(5,4,roi_map{r}); % plot the trace
    hold on;
    cur_t = t;
    plt_cond_1(end+1) = plot( cur_t, l_itrace_1, 'Color', 'r', 'LineWidth', 2 );
    plt_cond_2(end+1) = plot( cur_t, l_itrace_2, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--' );
    plot( cur_t, r_itrace_1, 'Color', 'g', 'LineWidth', 2 );
    plot( cur_t, r_itrace_2, 'Color', 'g', 'LineWidth', 2, 'LineStyle', '--' );
    
    xlim([0 max(cur_t)]);
    ylim([-0.2 0.75]);
    xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('dF/F');
    set(gca, 'FontSize', 14 );
    set(gca, 'FontWeight', 'bold');
    
    yy = ylim;
    y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
    hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
    set(gca,'children',circshift(get(gca,'children'),-1));
    set(hh, 'EdgeColor', 'None');
    
    colorindex = colorindex+1;   
end

figure(f1);
ax1 = subplot(5,4,roi_map{1}); % plot the trace

l_cond_1_num_trials = size( btraces_per_condition{ 1, ac.LEFT }( :, ac.VEL_YAW, : ), 1 );
l_cond_2_num_trials = size( btraces_per_condition{ 2, ac.LEFT }( :, ac.VEL_YAW, : ), 1 );
r_cond_1_num_trials = size( btraces_per_condition{ 1, ac.RIGHT }( :, ac.VEL_YAW, : ), 1 );
r_cond_2_num_trials = size( btraces_per_condition{ 2, ac.RIGHT }( :, ac.VEL_YAW, : ), 1 );

ll = legend( [ plt_cond_1(1), plt_cond_2(1) ], ...
    [ condition_trials_str{ 1 } '(l: ' num2str( l_cond_1_num_trials ) ' r: ' num2str( r_cond_1_num_trials ) ')'], ...
    [ condition_trials_str{ 2 } '(l: ' num2str( l_cond_2_num_trials ) ' r: ' num2str( r_cond_2_num_trials ) ')'], 'Location', 'southeast');
set(ll, 'Interpreter', 'none');

drawnow;

saveas(f1, [ filename_prefix '_left_right_tc.fig']);
saveas(f1, [ filename_prefix '_left_right_tc.png']);
saveas(f2, [ filename_prefix '_left_right_rois.fig']);
saveas(f2, [ filename_prefix '_left_right_rois.png']);

end

