function [ glom_EPG_dF_F_per_trial ] = display_EB_dynamics_var_stim( EPG_data, pico_stim_data, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid )

% EPG_data and PEN_data is 
% { trials, x, y, planes, time }

ac = get_analysis_constants;

% glom is a cell array of { planes, roi }, there will always be 8 glomeruli
% glom = get_left_half_PN_glomeruli( EPG_data, analysis_path );
% glom = get_left_half_PN_glomeruli_10202018( EPG_data, analysis_path );
glom = get_EB_glomeruli_16_planes( EPG_data, analysis_path );

ysize = size(EPG_data,2);
xsize = size(EPG_data,3);
nframes = size(EPG_data,5);     

num_trials = size(EPG_data, 1);

% Use trials_to_include from this point forward
glom_EPG_F_per_trial = get_PB_F_per_trial( glom, EPG_data );
glom_EPG_dF_F_per_trial = get_PB_dF_F_per_trial( glom_EPG_F_per_trial );

yaw_all = [];
yaw_time_all = [];

ephys_data_all = [];
ephys_time_all = [];

dt_bvel = bdata_vel_time(2) - bdata_vel_time(1);
dt_ephys_t = ephys_time(2) - ephys_time(1);

stim_all = [];
t = [0:nframes-1]./VPS;
    
time_all = [];
dt = t(2) - t(1);

for tr = 1:num_trials     
    
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
   
    stim_all = horzcat( stim_all, 10.0*pico_stim_data{1}( tr, : ) );
    
    if( length(time_all) == 0 )
        time_all = horzcat(time_all, t);
    else
        time_all = horzcat(time_all, time_all(end)+dt+t);
    end
end

EPG_roi_all_c = cell( 1, length(glom) );

for g = 1:length(glom)
    
    EPG_roi_all_c{ g } = [];
    
    for tr = 1:num_trials
        EPG_roi_all_c{ g } = horzcat( EPG_roi_all_c{ g }, squeeze(glom_EPG_dF_F_per_trial( tr, g, : ))' );
    end
end

EPG_roi_all = zeros( length(glom)+1, length(EPG_roi_all_c{ 1 }) );
for g = 1:length(glom)
    EPG_roi_all(g,:) = EPG_roi_all_c{ g };
end


% Add stimulus markers, as the last ROI
% downsample stim_all
imaging_data_len = length(EPG_roi_all_c{ 1 });
stim_all_len     = length(stim_all);

dstim = floor(stim_all_len / imaging_data_len);

stim_all_d = squeeze(mean(reshape(stim_all(1:dstim*imaging_data_len), [dstim, imaging_data_len])));

EPG_roi_all(end,:) = stim_all_d;

f = figure;

ax(1) = subplot(3,1,1);
imagesc( time_all, [1:size(EPG_roi_all,1)], EPG_roi_all );
colormap(flipud(gray));
caxis([-0.5 4]);
%colorbar;

ax(2) = subplot(3,1,2);
plot( yaw_time_all, yaw_all );

ax(3) = subplot(3,1,3);
plot( ephys_time_all, ephys_data_all );

linkaxes(ax,'x');

saveas(f,[analysis_path '/roi_traces_all_data_sid_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/roi_traces_all_data_sid_' num2str(sid) '.png']);

for tr = 1:num_trials
    f = figure;

    cur_stim = 10.0*squeeze(pico_stim_data{1}( tr, : ));
    dstim = floor(length(cur_stim) / nframes);
    stim_d = squeeze(mean(reshape(cur_stim(1:dstim*nframes), [dstim, nframes])));
        
    ax(1) = subplot( 3, 1, 1 );
    cur_EPG = squeeze(glom_EPG_dF_F_per_trial(tr, :,:));
    
    cur_EPG(end+1,:) = stim_d;
    
    imagesc( t, [1:size(glom_EPG_dF_F_per_trial,2)+1], cur_EPG );
    colormap(flipud(gray));
    caxis([-0.5 2]);
    %colorbar;
    
    ax(2) = subplot( 3, 1, 2 );
    plot( bdata_vel_time,  squeeze(bdata_vel{ 1 }( tr, ac.VEL_YAW, : )) );
    
    ax(3) = subplot( 3, 1, 3 );
    plot( ephys_time,  squeeze(ephys_data{ 1 }( tr, : )) );
    xlabel('Time (s)');
    
    linkaxes(ax,'x');
        
    saveas(f,[analysis_path '/roi_traces_sid_' num2str(sid) '_tid_' num2str(tr) '.fig']);
    saveas(f,[analysis_path '/roi_traces_sid_' num2str(sid) '_tid_' num2str(tr) '.png']);
    close(f);
end
end

