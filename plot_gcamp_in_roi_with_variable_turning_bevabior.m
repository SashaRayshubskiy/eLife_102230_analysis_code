function [ output_args ] = plot_gcamp_in_roi_with_variable_turning_behav(condition_trials, ref_img, PLANE_OF_INTEREST, TRIAL_TYPE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

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

trial_type = TRIAL_TYPE_OF_INTEREST;
p = PLANE_OF_INTEREST;

f1 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
f2 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
               
cur_plane_avg_df_f_cond_1 = squeeze(avg_df_f_per_condition_per_plane(trial_type,1,p,:,:,:));
cur_plane_avg_df_f_cond_1(~isfinite(cur_plane_avg_df_f_cond_1)) = 0.0;

cur_plane_avg_df_f_cond_2 = squeeze(avg_df_f_per_condition_per_plane(trial_type,2,p,:,:,:));
cur_plane_avg_df_f_cond_2(~isfinite(cur_plane_avg_df_f_cond_2)) = 0.0;
        
cur_t = t;

% Extract frames during stim only for now.
cur_frames = find((cur_t >= prestim) & (cur_t<=(prestim+stim)));

avg_df_f_img_cond_1 = squeeze(mean(cur_plane_avg_df_f_cond_1(:,:,cur_frames),3));
avg_df_f_img_cond_2 = squeeze(mean(cur_plane_avg_df_f_cond_2(:,:,cur_frames),3));

[xsize, ysize] = size(ref_img);
figure(f1)
imagesc( ref_img );
colormap(ax1, 'gray');
axis image;
caxis([0 300]);
title([ac.task_str{trial_type}]);
    
a_data_1 = cur_plane_avg_df_f_cond_1;
a_data_2 = cur_plane_avg_df_f_cond_2;

plt_cond_1 = [];
plt_cond_2 = [];

clicky_plane = 2;
while(npts > 0)
    
    % subplot(1,3,1)
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
        break
    end
    inpoly = inpolygon(x,y,xv,yv);
        
    subplot(2,2,clicky_plane);
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
    subplot(2,1,1);
    hold on;
    cur_t = t;
    plt_cond_1(end+1) = plot( cur_t, itrace_1, 'Color', currcolor, 'LineWidth', 2);
    plt_cond_2(end+1) = plot( cur_t, itrace_2, 'Color', currcolor, 'LineWidth', 2, 'LineStyle', '--');
    
    xlim([0 max(cur_t)]);
    ylim([-0.2 0.75]);
    xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('dF/F');
    set(gca, 'FontSize', 14 );
    set(gca, 'FontWeight', 'bold');
    
    colorindex = colorindex+1;
    
    roi_points{nroi} = [xv, yv];
    nroi = nroi + 1;
end

figure(f2);
subplot(2,1,1);
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

saveas(f1, [ filename_prefix '_' ac.task_str{trial_type} '_rois_with_behav.fig']);
saveas(f1, [ filename_prefix '_' ac.task_str{trial_type} '_rois_with_behav.png']);

trial_type = 1;
yaw_traces = [];
fwd_traces = [];
for trial_ord = 1:length(condition_trials{trial_type, cond_ord})
    
    cur_trial = condition_trials{ trial_type, cond_ord }( trial_ord );
    yaw_trace = squeeze(bdata_vel{trial_type}( cur_trial, ac.VEL_YAW, :));
    fwd_trace = squeeze(bdata_vel{trial_type}( cur_trial, ac.VEL_FWD, :));
    
    yaw_traces(end+1,:) = yaw_trace;
    fwd_traces(end+1,:) = fwd_trace;
end

% Plot average trial
avg_yaw_trace = mean( yaw_traces );
sem_yaw_trace = std( yaw_traces, 1 ) ./ sqrt( size(yaw_traces,1) );

figure(f2);
subplot(2,1,2);
fh = fill( [bdata_vel_time, fliplr(bdata_vel_time)], ...
    [(avg_yaw_trace+sem_yaw_trace) fliplr((avg_yaw_trace-sem_yaw_trace))], ...
    cur_color_single);
set(fh, 'EdgeColor', 'None');

plot( bdata_vel_time, avg_yaw_trace, 'color', cur_color_avg, 'LineStyle', cur_cond_symbol, 'LineWidth', 2.0 );
        
saveas(f2, [ filename_prefix '_' ac.task_str{trial_type} '_tc_with_behav.fig']);
saveas(f2, [ filename_prefix '_' ac.task_str{trial_type} '_tc_with_behav.png']);

end

