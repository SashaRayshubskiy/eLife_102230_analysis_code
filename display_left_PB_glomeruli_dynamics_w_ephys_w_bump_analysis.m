function [ glom_EPG_dF_F_per_trial ] = display_left_PB_glomeruli_dynamics_w_ephys_w_bump_analysis( EPG_data, PEN_data,  ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid, pre_stim, stim_dur )

% EPG_data and PEN_data is 
% { trials, x, y, planes, time }

ac = get_analysis_constants;

% glom is a cell array of { planes, roi }, there will always be 8 glomeruli
glom = get_left_half_PN_glomeruli( EPG_data, analysis_path );

ysize = size(EPG_data,2);
xsize = size(EPG_data,3);
nframes = size(EPG_data,5);     

num_trials = size(EPG_data, 1);

% Use trials_to_include from this point forward
glom_EPG_F_per_trial = get_PB_F_per_trial( glom, EPG_data );
glom_EPG_dF_F_per_trial = get_PB_dF_F_per_trial( glom_EPG_F_per_trial );

glom_PEN_F_per_trial = get_PB_F_per_trial( glom, PEN_data );
glom_PEN_dF_F_per_trial = get_PB_dF_F_per_trial( glom_PEN_F_per_trial );

yaw_all = [];
yaw_time_all = [];

ephys_data_all = [];
ephys_time_all = [];

dt_bvel = bdata_vel_time(2) - bdata_vel_time(1);
dt_ephys_t = ephys_time(2) - ephys_time(1);

post_stim = pre_stim + stim_dur;

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
   
    stim_all = horzcat( stim_all, stim );
    if( length(time_all) == 0 )
        time_all = horzcat(time_all, t);
    else
        time_all = horzcat(time_all, time_all(end)+dt+t);
    end
end

EPG_roi_all_c = cell( 1, length(glom) );
PEN_roi_all_c = cell( 1, length(glom) );

for g = 1:length(glom)
    
    EPG_roi_all_c{ g } = [];
    PEN_roi_all_c{ g } = [];
    
    for tr = 1:num_trials
        EPG_roi_all_c{ g } = horzcat( EPG_roi_all_c{ g }, squeeze(glom_EPG_dF_F_per_trial( tr, g, : ))' );
        PEN_roi_all_c{ g } = horzcat( PEN_roi_all_c{ g }, squeeze(glom_PEN_dF_F_per_trial( tr, g, : ))' );
    end
end

EPG_roi_all = zeros( length(glom)+1, length(EPG_roi_all_c{ 1 }) );
PEN_roi_all = zeros( length(glom)+1, length(EPG_roi_all_c{ 1 }) );
for g = 1:length(glom)
    EPG_roi_all(g,:) = EPG_roi_all_c{ g };
    PEN_roi_all(g,:) = PEN_roi_all_c{ g };
end


% Add stimulus markers, as the last ROI
EPG_roi_all(end,:) = stim_all;
PEN_roi_all(end,:) = stim_all;

f = figure;

ax(1) = subplot(4,1,1);
imagesc( time_all, [1:size(EPG_roi_all,1)], EPG_roi_all );
colormap(flipud(gray));
caxis([-0.5 4]);
%colorbar;

ax(2) = subplot(4,1,2);
imagesc( time_all, [1:size(EPG_roi_all,1)], PEN_roi_all );
colormap(flipud(gray));
caxis([-0.5 10]);
%colorbar;

ax(3) = subplot(4,1,3);
plot( yaw_time_all, yaw_all );

ax(4) = subplot(4,1,4);
plot( ephys_time_all, ephys_data_all );

linkaxes(ax,'x');

saveas(f,[analysis_path '/roi_traces_all_data_sid_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/roi_traces_all_data_sid_' num2str(sid) '.png']);

stim_start = pre_stim;
stim_stop = post_stim;

f = figure;

avg_yaw = squeeze(mean(squeeze(bdata_vel{ 1 }( :, ac.VEL_YAW, : ))));
avg_fwd = squeeze(mean(squeeze(bdata_vel{ 1 }( :, ac.VEL_FWD, : ))));
avg_ephys = squeeze(mean(ephys_data{ 1 }));

ax_1(1) = subplot(3,1,1); 
hold on;
plot( bdata_vel_time, avg_fwd );
xlim([0 bdata_vel_time(end)]);
ylabel('Fwd vel (au/s)');

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ stim_start stim_start stim_stop stim_stop ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

ax_1(2) = subplot(3,1,2);
hold on;
plot( bdata_vel_time, avg_yaw );
xlim([0 bdata_vel_time(end)]);
ylabel('Yaw (au/s)');
axis tight;

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ stim_start stim_start stim_stop stim_stop ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

ax_1(3) = subplot(3,1,3);
hold on;
plot( ephys_time, avg_ephys );
xlabel('Time (s)');
ylabel('Vm (mV)');
xlim([0 ephys_time(end)]);
axis tight;

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ stim_start stim_start stim_stop stim_stop ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

linkaxes(ax_1,'x');

saveas(f,[analysis_path '/fwd_yaw_tc_sid_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/fwd_yaw_tc_sid_' num2str(sid) '.png']);

f = figure;
hold on;

for g = 1:length(glom)
    avg_in_glom = mean(squeeze(glom_EPG_dF_F_per_trial( :, g, : )));
    
    plot( t, avg_in_glom, 'DisplayName', [ 'glom: ' num2str(g)] );
end

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ stim_start stim_start stim_stop stim_stop ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

legend();
saveas(f,[analysis_path '/avg_tc_in_roi_sid_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/avg_tc_in_roi_sid_' num2str(sid) '.png']);

end

