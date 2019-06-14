function [ PFL_LAL_dF_F_per_trial, PFL_FB_dF_F_per_trial ] = display_PFLm_dynamics( imaging_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid )

% imaging_data and PEN_data is 
% { trials, x, y, planes, time }

ac = get_analysis_constants;

% rois is a cell array of { planes, roi }, there will always be 10 rois
% { left PFL.LAL, right PFL.LAL, 8 FB glomeruli starting left to right }

rois = get_PFLm_rois( imaging_data, analysis_path );

ysize = size(imaging_data,2);
xsize = size(imaging_data,3);
nframes = size(imaging_data,5);     

num_trials = size(imaging_data, 1);

% Use trials_to_include from this point forward
PFL_LAL_F_per_trial = get_PB_F_per_trial( rois(1:4, :), imaging_data );
PFL_LAL_dF_F_per_trial = get_PB_dF_F_per_trial( PFL_LAL_F_per_trial );

PFL_FB_F_per_trial = get_PB_F_per_trial( rois(5:end, :), imaging_data );
PFL_FB_dF_F_per_trial = get_PB_dF_F_per_trial( PFL_FB_F_per_trial );

yaw_all = [];
yaw_time_all = [];

dt_bvel = bdata_vel_time(2) - bdata_vel_time(1);

stim_all = [];
t = [0:nframes-1]./VPS;
    
time_all = [];
dt = t(2) - t(1);

for tr = 1:num_trials     
    
    cur_yaw = squeeze(bdata_vel{ 1 }( tr, ac.VEL_YAW, : ));       
    yaw_all = horzcat( yaw_all, cur_yaw' );
    
    if( length(yaw_time_all) == 0 )
        yaw_time_all = horzcat( yaw_time_all, bdata_vel_time );
    else
        yaw_time_all = horzcat( yaw_time_all, bdata_vel_time + (yaw_time_all(end)+dt_bvel) );        
    end
       
    if( length(time_all) == 0 )
        time_all = horzcat(time_all, t);
    else
        time_all = horzcat(time_all, time_all(end)+dt+t);
    end
end

% right/left PFL.LAL tuft and axon
NUM_LAL_ROIS = 4;

% Fan-shaped body glomeruli
NUM_FB_ROIS = 8;

PFL_LAL_roi_all_c = cell( 1, length(NUM_LAL_ROIS));
PFL_FB_roi_all_c = cell( 1, length(NUM_FB_ROIS));

for g = 1:NUM_LAL_ROIS    
    PFL_LAL_roi_all_c{ g } = [];
    
    for tr = 1:num_trials
        PFL_LAL_roi_all_c{ g } = horzcat( PFL_LAL_roi_all_c{ g }, squeeze(PFL_LAL_dF_F_per_trial( tr, g, : ))' );
    end
end

for g = 1:NUM_FB_ROIS
    PFL_FB_roi_all_c{ g } = [];
    
    for tr = 1:num_trials
        PFL_FB_roi_all_c{ g } = horzcat( PFL_FB_roi_all_c{ g }, squeeze( PFL_FB_dF_F_per_trial( tr, g, : ))' );
    end
end

PFL_LAL_roi_all = zeros( NUM_LAL_ROIS, length(PFL_LAL_roi_all_c{ 1 }) );
for g = 1:NUM_LAL_ROIS
    PFL_LAL_roi_all(g,:) = PFL_LAL_roi_all_c{ g };
end

PFL_FB_roi_all = zeros( NUM_FB_ROIS, length(PFL_FB_roi_all_c{ 1 }) );
for g = 1:NUM_FB_ROIS
    PFL_FB_roi_all(g,:) = PFL_FB_roi_all_c{ g };
end

% Add stimulus markers, as the last ROI
% downsample stim_all
imaging_data_len = size( PFL_LAL_roi_all, 2 );

f = figure;

ax(1) = subplot(3,1,1);
imagesc( time_all, [1:size(PFL_LAL_roi_all,1)], PFL_LAL_roi_all );
colormap(flipud(gray));
caxis([-0.5 4]);

ax(2) = subplot(3,1,2);
imagesc( time_all, [1:size(PFL_FB_roi_all,1)], PFL_FB_roi_all );
colormap(flipud(gray));
caxis([-0.5 4]);

ax(3) = subplot(3,1,3);
plot( yaw_time_all, yaw_all );

linkaxes(ax,'x');

saveas(f,[analysis_path '/roi_traces_all_data_sid_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/roi_traces_all_data_sid_' num2str(sid) '.png']);

for tr = 1:num_trials
    f = figure;
        
    ax( 1 ) = subplot( 3, 1, 1 );
    cur_PFL_LAL = squeeze(PFL_LAL_dF_F_per_trial(tr, :,:));        
    imagesc( t, [1:size(PFL_LAL_dF_F_per_trial,2)], cur_PFL_LAL );
    colormap(flipud(gray));
    caxis([-0.5 2]);
    
    ax( 2 ) = subplot( 3, 1, 2 );
    cur_PFL_FB = squeeze(PFL_FB_dF_F_per_trial(tr, :,:));        
    imagesc( t, [1:size(PFL_FB_dF_F_per_trial,2)], cur_PFL_FB );
    colormap(flipud(gray));
    caxis([-0.5 2]);
    
    ax( 3 ) = subplot( 3, 1, 3 );
    plot( bdata_vel_time,  squeeze(bdata_vel{ 1 }( tr, ac.VEL_YAW, : )) );
    
    xlabel('Time (s)');    
    linkaxes(ax,'x');
        
    saveas(f,[analysis_path '/roi_traces_sid_' num2str(sid) '_tid_' num2str(tr) '.fig']);
    saveas(f,[analysis_path '/roi_traces_sid_' num2str(sid) '_tid_' num2str(tr) '.png']);
    close(f);
end
end

