function display_PB_roi_dynamics( cdata_raw, bdata_vel_time, bdata_vel, VPS, analysis_path )

% trials, x, y, color, plane, time
EPG_data       = squeeze(cdata_raw{1}( :, :, :, 1, :, : ));
PEN_tdTom_data = squeeze(cdata_raw{1}( :, :, :, 2, :, : ));

% Get max for each plane
EPG_data_max = max(EPG_data, [], 4);
PEN_data_max = max(PEN_tdTom_data, [], 4);

% Get max for each trial
EPG_data_max_tr = max( EPG_data_max, [], 1 );
PEN_data_max_tr = max( PEN_data_max, [], 1 );


rois = get_rois(squeeze(mean(EPG_data_max_tr,3)), analysis_path);

ysize = size(EPG_data_max,2);
xsize = size(EPG_data_max,3);
nframes = size(EPG_data,5);
[x, y] = meshgrid(1:xsize, 1:ysize);

num_trials = size(EPG_data, 1);

% Compute F 
EPG_roi_min = cell(1,length(rois));
PEN_roi_min = cell(1,length(rois));
for r = 1:length(rois)
    EPG_roi_min{r} = [];
    PEN_roi_min{r} = [];
end

for tr = 1:num_trials

    % get dF/f for an roi
    for r = 1:length(rois)
        cur_roi = rois{r};
        
        xv = cur_roi(:,1);
        yv = cur_roi(:,2);
    
        inpoly = inpolygon(x,y,xv,yv);
        
        tmp = squeeze(sum(sum(squeeze(EPG_data_max(tr,:,:,:)).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));        
        EPG_roi_min{r} = horzcat(EPG_roi_min{r}, tmp');

        tmp = squeeze(sum(sum(squeeze(PEN_data_max(tr,:,:,:)).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
        PEN_roi_min{r} = horzcat(PEN_roi_min{r}, tmp');
    end 
end

EPG_roi_all = zeros(length(rois), length(EPG_roi_min{1}));
PEN_roi_all = zeros(length(rois), length(EPG_roi_min{1}));
for r = 1:length(rois)
    EPG_roi_all(r,:) = EPG_roi_min{r};
    PEN_roi_all(r,:) = PEN_roi_min{r};
end


f = figure;
subplot(2,1,1);
imagesc( EPG_roi_all );

subplot(2,1,2);
imagesc( PEN_roi_all );

N_min = ceil(0.05 * length(EPG_roi_min{1}));
for r = 1:length(rois)
    xx = sort(EPG_roi_min{r});
    EPG_min(r) = mean(xx(1:N_min));

    xx = sort(PEN_roi_min{r});
    PEN_min(r) = mean(xx(1:N_min));
end

for tr = 1:num_trials

    % get dF/f for an roi
    for r = 1:length(rois)
        cur_roi = rois{r};
        
        xv = cur_roi(:,1);
        yv = cur_roi(:,2);
    
        inpoly = inpolygon(x,y,xv,yv);
        
        tmp = squeeze(sum(sum(squeeze(EPG_data_max(tr,:,:,:)).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));        
        baseline = repmat( EPG_min(r), [1 1 size(tmp,2)]);
        EPG_roi_trace(tr,r,:) = (tmp-baseline) ./ baseline;

        tmp = squeeze(sum(sum(squeeze(PEN_data_max(tr,:,:,:)).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
        baseline = repmat( PEN_min(r), [1 1 size(tmp,2)]);
        PEN_roi_trace(tr,r,:) = (tmp-baseline) ./ baseline;
    end 
end

t = [0:nframes-1]./VPS;

% Identify the stim ROI using the PEN tdTom data
stim_t = find( (t < 4.0) & (t>3.0) );

[dummy, stim_roi] = max(mean(PEN_data_max(:,stim_t),2));

pre_stim_t = find( (t < 3.0) & (t>2.0) );

left_trials = [];
right_trials = [];

for tr = 1:num_trials
    cur_trial = squeeze(EPG_roi_trace(tr,:,:));

    pre_stim_per_roi = squeeze(mean(cur_trial(:,pre_stim_t),2));
    [y,i] = max(pre_stim_per_roi);
    
    if( i < stim_roi )
        left_trials(end+1) = tr;
    else
        right_trials(end+1) = tr;
    end
    
end

ac = get_analysis_constants;

left_yaw = [];
left_fwd = [];
for tr = 1:length( left_trials )    
    cur_trial_idx = left_trials(tr);
    
    left_fwd(tr, :) = bdata_vel{ 1 }( cur_trial_idx, ac.VEL_FWD, : );
    left_yaw(tr, :) = bdata_vel{ 1 }( cur_trial_idx, ac.VEL_YAW, : );
end

right_yaw = [];
right_fwd = [];
for tr = 1:length( right_trials )    
    cur_trial_idx = right_trials(tr);
    
    right_fwd(tr, :) = bdata_vel{ 1 }( cur_trial_idx, ac.VEL_FWD, : );
    right_yaw(tr, :) = bdata_vel{ 1 }( cur_trial_idx, ac.VEL_YAW, : );
end

f = figure;
subplot(2,1,1);
hold on;
plot(bdata_vel_time, mean(left_fwd), 'DisplayName', 'EPG left of stim');
plot(bdata_vel_time, mean(right_fwd), 'DisplayName', 'EPG right of stim');
legend();
ylabel('Fwd vel (mm/s)');
xlabel('Time (s)');

subplot(2,1,2);
hold on;
plot(bdata_vel_time, mean(left_yaw), 'DisplayName', 'EPG left of stim');
plot(bdata_vel_time, mean(right_yaw), 'DisplayName', 'EPG right of stim');
ylabel('Yaw vel (deg/s)');
xlabel('Time (s)');

saveas(f,[analysis_path '/behavior_by_epg_phase.fig']);
saveas(f,[analysis_path '/behavior_by_epg_phase.png']);

