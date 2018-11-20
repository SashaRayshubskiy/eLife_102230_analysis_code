function display_PB_roi_dynamics_w_ephys( EPG_data, PEN_data,  ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid )

SKIP_PEN = 0;
if(length(PEN_data) == 0)
    SKIP_PEN = 1;
end

ac = get_analysis_constants;

rois = get_rois(squeeze(mean(squeeze(mean(EPG_data(:,:,:,:),4)),1)), analysis_path);

ysize = size(EPG_data,2);
xsize = size(EPG_data,3);
nframes = size(EPG_data,4);
[x, y] = meshgrid(1:xsize, 1:ysize);        

num_trials = size(EPG_data, 1);
%num_trials = 3;

EPG_smallest_F_tmp = [];
PEN_smallest_F_tmp = [];
EPG_smallest_F = [];
PEN_smallest_F = [];

% Filter out 'broken scanimage' trials
avg_signal = [];

for r = 1:length(rois)
    cur_roi = rois{r};
    
    xv = cur_roi(:,1);
    yv = cur_roi(:,2);
            
    inpoly = inpolygon(x,y,xv,yv);
    
    for tr = 1:num_trials                
        cur_F_tc_in_roi_EPG = squeeze(sum(sum(squeeze(EPG_data(tr,:,:,:)).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));              
        
        avg_signal(r, tr) = mean( cur_F_tc_in_roi_EPG );
    end    
end

avg_signal_all_roi = squeeze(mean(avg_signal));
%figure;
%plot(avg_signal_all_roi);

BAD_SCANIMAGE_TRIAL_THRESHOLD = 450;

trials_to_include = [];

for tr = 1:num_trials
    %%if(avg_signal_all_roi(tr) > BAD_SCANIMAGE_TRIAL_THRESHOLD )
        trials_to_include(end+1) = tr;        
    %%end
end

% Use trials_to_include from this point forward
roi_ts_per_trial = zeros( length(rois), length( trials_to_include ), nframes );

for r = 1:length(rois)
    cur_roi = rois{r};
    
    xv = cur_roi(:,1);
    yv = cur_roi(:,2);
            
    inpoly = inpolygon(x,y,xv,yv);

    
    tmp_list_EPG = [];
    tmp_list_PEN = [];
    
    for j = 1:length(trials_to_include)
    
        tr = trials_to_include(j);
                
        cur_F_tc_in_roi_EPG = squeeze(sum(sum(squeeze(EPG_data(tr,:,:,:)).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));  
        tmp_list_EPG = horzcat( tmp_list_EPG, cur_F_tc_in_roi_EPG' );
        
        if( SKIP_PEN == 0 )
            cur_F_tc_in_roi_PEN = squeeze(sum(sum(squeeze(PEN_data(tr,:,:,:)).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));  
            tmp_list_PEN = horzcat( tmp_list_PEN, cur_F_tc_in_roi_PEN' );
        end
                   
    end
    
    EPG_smallest_F(r,:) = tmp_list_EPG;
    
    if( SKIP_PEN == 0 )
        PEN_smallest_F(r,:) = tmp_list_PEN;
    end
end

% Get the smallest 5% of the time course data for each ROI
baseline_per_roi_PEN = [];
baseline_per_roi_EPG = [];
for r = 1:length(rois)
    baseline_per_roi_EPG( r ) = mean( EPG_smallest_F(r) ) - 2*std( EPG_smallest_F(r), 1 );
    
    if( SKIP_PEN == 0 )
        baseline_per_roi_PEN( r ) = mean( PEN_smallest_F(r) ) - 2*std( PEN_smallest_F(r), 1 );
    end
end

yaw_all = [];
yaw_time_all = [];

ephys_data_all = [];
ephys_time_all = [];

dt_bvel = bdata_vel_time(2) - bdata_vel_time(1);
dt_ephys_t = ephys_time(2) - ephys_time(1);

pre_stim = 5.0;
post_stim = pre_stim + 0.05;

for j = 1:length(trials_to_include)
    
    tr = trials_to_include(j);
    
    cur_yaw = squeeze(bdata_vel{ 1 }( tr, ac.VEL_YAW, : ));       
    yaw_all = horzcat( yaw_all, cur_yaw' );
    
    cur_ephys = squeeze(ephys_data{ 1 }( tr, : ));       
    ephys_data_all = horzcat( ephys_data_all, cur_ephys' );

    if( length(yaw_time_all) == 0 )
        yaw_time_all = horzcat( yaw_time_all, bdata_vel_time );
        ephys_time_all = horzcat( ephys_time_all, ephys_time );
    else
        yaw_time_all = horzcat( yaw_time_all, bdata_vel_time + (yaw_time_all(end)+dt_bvel) );        
        ephys_time_all = horzcat( ephys_time_all, (ephys_time + ephys_time_all(end)+dt_ephys_t) );        
    end
end

stim_all = [];
t = [0:nframes-1]./VPS;
stim_t = find( (t >= pre_stim) & (t <=  post_stim) );

stim = zeros(1, nframes);

STIM_VAL = 50;
if( length(stim_t) == 0 )
    stim_t = find( (t >= (pre_stim-0.2)) & (t <=  post_stim) );   
    stim(stim_t(1)) = STIM_VAL;
else
    stim(stim_t) = STIM_VAL;
end

EPG_roi_min = {};
PEN_roi_min = {};
for r = 1:length(rois)
    EPG_roi_min{r} = [];
    PEN_roi_min{r} = [];    
end
    
time_all = [];
dt = t(2) - t(1);
for j = 1:length(trials_to_include)
    
    tr = trials_to_include(j);
    
    % get dF/f for an roi
    for r = 1:length(rois)
        cur_roi = rois{r};
        
        xv = cur_roi(:,1);
        yv = cur_roi(:,2);
    
        inpoly = inpolygon(x,y,xv,yv);
        
        tmp = squeeze(sum(sum(squeeze(EPG_data(tr,:,:,:)).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));  
        F = baseline_per_roi_EPG(r);
        df_f = (tmp - F) ./ F;
        EPG_roi_min{r} = horzcat(EPG_roi_min{r}, df_f');
        
        roi_ts_per_trial( r, j, : ) = df_f;

        if( SKIP_PEN == 0 )
            tmp = squeeze(sum(sum(squeeze(PEN_data(tr,:,:,:)).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
            F = baseline_per_roi_PEN(r);
            df_f = (tmp - F) ./ F;
            PEN_roi_min{r} = horzcat(PEN_roi_min{r}, df_f');
        end
    end
    
    stim_all = horzcat( stim_all, stim );
    if( length(time_all) == 0 )
        time_all = horzcat(time_all, t);
    else
        time_all = horzcat(time_all, time_all(end)+dt+t);
    end
end

EPG_roi_all = zeros(length(rois), length(EPG_roi_min{1}));
PEN_roi_all = zeros(length(rois), length(EPG_roi_min{1}));
for r = 1:length(rois)
    EPG_roi_all(r,:) = EPG_roi_min{r};

    if( SKIP_PEN == 0 )
        PEN_roi_all(r,:) = PEN_roi_min{r};
    end
end

% Add stimulus markers, as the last ROI
EPG_roi_all(end+1,:) = stim_all;
PEN_roi_all(end+1,:) = stim_all;

f = figure;

if( SKIP_PEN == 0 )    
    ax(1) = subplot(4,1,1);
    imagesc( time_all, [1:size(EPG_roi_all,1)], EPG_roi_all );
    colormap(flipud(gray));
    caxis([-0.5 1.5]);
    %colorbar;
    
    ax(2) = subplot(4,1,2);
    imagesc( time_all, [1:size(EPG_roi_all,1)], PEN_roi_all );
    colormap(flipud(gray));
    caxis([-0.5 1.5]);
    %colorbar;
    
    ax(3) = subplot(4,1,3);
    plot( yaw_time_all, yaw_all );

    ax(4) = subplot(4,1,4);
    plot( ephys_time_all, ephys_data_all );        
else
    ax(1) = subplot(2,1,1);
    colormap(flipud(gray));
    imagesc( time_all, [1:size(EPG_roi_all,1)], EPG_roi_all );
    caxis([-0.5 1.5]);
    %colorbar;
        
    ax(2) = subplot(2,1,2);
    plot( yaw_time_all, yaw_all);
end    

linkaxes(ax,'x');

saveas(f,[analysis_path '/roi_traces_all_data_sid_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/roi_traces_all_data_sid_' num2str(sid) '.png']);

f = figure;

avg_yaw = squeeze(mean(squeeze(bdata_vel{ 1 }( :, ac.VEL_YAW, : ))));
avg_fwd = squeeze(mean(squeeze(bdata_vel{ 1 }( :, ac.VEL_FWD, : ))));

ax_1(1) = subplot(2,1,1); 
plot( bdata_vel_time, avg_fwd );
xlim([0 bdata_vel_time(end)]);
ylabel('Fwd vel (au/s)');

ax_1(2) = subplot(2,1,2);
plot( bdata_vel_time, avg_yaw );
xlim([0 bdata_vel_time(end)]);
xlabel('Time (s)');
ylabel('Yaw (au/s)');
axis tight;

linkaxes(ax_1,'x');

saveas(f,[analysis_path '/fwd_yaw_tc_sid_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/fwd_yaw_tc_sid_' num2str(sid) '.png']);

f = figure;
hold on;

for r = 1:length(rois)
    avg_in_roi = mean(squeeze(roi_ts_per_trial( r, :, : )));
    
    plot( t, avg_in_roi, 'DisplayName', [ 'roi: ' num2str(r)] );
end

stim_start = pre_stim;
stim_stop = post_stim;
yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ stim_start stim_start stim_stop stim_stop ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

legend();
saveas(f,[analysis_path '/avg_tc_in_roi_sid_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/avg_tc_in_roi_sid_' num2str(sid) '.png']);


end

