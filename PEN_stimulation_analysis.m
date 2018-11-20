%% Load imaging and behavioral data
global slash;

clear all;

if isunix() == 1
    slash = '/';
else
    slash = '\';
end

%trial_exclusion_list = nsyb_83blexA_01_blank_trials;
trial_exclusion_list = {[],[],[]};

%datapath = '/data/drive2/sasha/180706_R12D09_P2X2_op6s_60D05_09/';
%datapath = '/data/drive2/sasha/180707_R12D09_lexA_LexAOp_GFP_01/';
datapath = '/data/drive2/sasha/181022_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15/';
%datapath = '/data/drive2/sasha/181019_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_13_first_trials/';

analysis_path = [datapath slash 'analysis'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

sid = [ 0 ];

aconstants = get_analysis_constants;
trial_type_cnt = 1; 

% Load behavioral data
bdata_path = [datapath  slash 'ball' slash ];
tic; [ b_rawdata, b_time, btrial_meta ] = load_behavioral_data(sid, bdata_path, trial_type_cnt ); toc

% Get behavioral data that is usable for analysis
%[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data( b_time, b_rawdata );
[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data_berg1( b_time, b_rawdata );

% Get ephys data
[ephys_time, ephys_data] = reformat_ephys_data_berg1( b_time, b_rawdata );

% Get pico monitor data
[ pico_stim_data ] = reformat_pico_stim_data_berg1( b_rawdata );

% Load imaging data
cdata_path = [datapath  slash '2p' slash ];
dx = 1;
dy = 1;
dt = 1;
down_params(1) = dx;
down_params(2) = dy;

LOAD_IMAGE_DATA = 1;
if(LOAD_IMAGE_DATA == 1)
    if( dt == 1 )
        % tic; [ cdata_raw, cdata_meta, ctrial_meta ] = load_imaging_data_2(sid, cdata_path, trial_type_cnt, dx, dy, dt ); toc
        tic; [ cdata_raw, cdata_meta, ctrial_meta ] = load_imaging_data_hack(sid, cdata_path, trial_type_cnt, dx, dy, dt ); toc
    else
        tic; [ cdata_raw, cdata_meta, ctrial_meta ] = load_imaging_data_single_plane(sid, cdata_path, trial_type_cnt, dx, dy, dt ); toc
    end
    
    
    settings = sensor_settings;
    global file_writer_cnt;
    file_writer_cnt = 1;
    
    % Create the reverse mapping from
    external_trial_id_to_internal_ordinal_map = get_external_trial_id_to_internal_ordinal_map(btrial_meta);
    
    % Create a map of frame start offsets per plane. (This is key because
    % there's about 150 ms interval from first plane to first plane in a
    % volume.
    if(dt == 1)
        planes = size( cdata_raw{ 1 }, 4 );
        VPS = cdata_meta.volume_rate;
        frame_start_offsets_per_plane = generate_frame_start_offsets_per_plane( planes, b_rawdata, b_time );
        
        ref_imgs = generate_ref_imgs(cdata_raw);
    else
        VPS = cdata_meta.frame_rate;
        %frame_start_offsets_per_plane = generate_frame_start_offsets_per_plane( 1, b_rawdata, b_time );
        frame_start_offsets_per_plane = generate_frame_start_offsets_per_plane( 1, b_rawdata, b_time );
        
        % frame_start_offsets_per_plane = [1:16] .* 9.8/1000;
        
        ref_imgs = squeeze(mean(cdata_raw{ 1 }(1,:,:,:),4));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finished data load.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% USE FOR 2 CHANNEL DATA: Pre-process time consuming parts of the analysis
% Drop the first 4 planes
% trials, x, y, color, plane, time
START_PLANE = 4;
EPG_data_1       = squeeze(cdata_raw{1}( :, :, :, 1, START_PLANE:end, : ));
% Get max for each plane
EPG_data = squeeze(max(EPG_data_1, [], 4));
% Dont run if only one channel saved
PEN_Alexa_data = squeeze(cdata_raw{1}( :, :, :, 2, START_PLANE:end, : ));
PEN_data = squeeze(max(PEN_Alexa_data, [], 4));

%% USE FOR 1 CHANNEL DATA: Pre-process time consuming parts of the analysis
% Drop the first 4 planes
% trials, x, y, color, plane, time
EPG_data_1       = squeeze(cdata_raw{1}( :, :, :, 5:end, : ));
% Get max for each plane
EPG_data = squeeze(max(EPG_data_1, [], 4));

%% 2 channel 
display_PB_roi_dynamics( EPG_data, PEN_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid );

%% 2 Channel with ephys
display_PB_roi_dynamics_w_ephys( EPG_data, PEN_data, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid );

%% 2 Channel with ephys and analysis of bump dynamics with left/right trial classification.

pre_stim = 15.0;
stim_dur = 0.05;        

[df_f_in_roi_per_trial] = display_PB_roi_dynamics_w_ephys_w_bump_analysis( EPG_data, PEN_data, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid, pre_stim, stim_dur );

%% Special left half PB glomeruli code, 2 channel with ephys and analysis of bump dynamics with left/right trial classification.

pre_stim = 15.0;
stim_dur = 0.05;        

[df_f_in_roi_per_trial] = display_left_PB_glomeruli_dynamics_w_ephys_w_bump_analysis( EPG_data_1, PEN_Alexa_data, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid, pre_stim, stim_dur );

%% Code for manual pico stim

[df_f_in_roi_per_trial] = display_left_PB_glomeruli_dynamics_var_stim( EPG_data_1, PEN_Alexa_data, pico_stim_data, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid );

%% 
display_bump_dynamics_turning_ephys( df_f_in_roi_per_trial, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid, pre_stim, stim_dur );

%% Dump dynamics with yaw and ephys with variable stim

display_bump_dynamics_turning_ephys_var_stim( df_f_in_roi_per_trial, pico_stim_data, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid );


%% 1 channel
display_PB_roi_dynamics( EPG_data, [], bdata_vel_time, bdata_vel, VPS, analysis_path, sid );

%%  Analyze Vm vs. yaw
FILT_FACTOR = 0.04;
ac = get_analysis_constants;

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;
BIN_SIZE = 0.050; % s
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

% { {trial, {t_start, t_finish} } } in seconds
% 181019_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_13_first_trials
% counter_turn_epochs =  { { 1, { 6.02, 7.025 } }, ...
%                          { 1, { 24.75, 25.4 } }, ...
%                          { 1, { 47.2, 48.2 }  }, ...
%                          { 2, { 8.6, 9.6 } } };

% 181022_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15
counter_turn_epochs =  { { 2, { 11.03, 13.32 } }, ...
                         { 3, { 16.07, 17.32 } }, ...
                         { 4, { 24.26, 27.45 }  }, ...
                         { 8, { 46.64, 48.28 } }, ...
                         { 9, { 18.3, 19.52 } }, ...
                         { 10, { 19.11, 20.83 } }};

f = figure;
for tr = 1:size(ephys_data{1},1)
    cur_ephys = squeeze(ephys_data{1}(tr,:));
    cur_yaw   = squeeze(bdata_vel{ 1 }( tr, ac.VEL_YAW, : ));            
    
    VmFilt_A2 = medfilt1( cur_ephys, FILT_FACTOR * ephys_SR, 'truncate' );
    VmFilt_A2_corr = VmFilt_A2 - mean(VmFilt_A2);

    t_down = squeeze(mean(reshape(ephys_time, [DT_EPHYS, length(ephys_time)/DT_EPHYS]), 1));
    A2_Vm_down = squeeze(mean(reshape(VmFilt_A2_corr, [ DT_EPHYS, length(VmFilt_A2)/DT_EPHYS ] ),1));
    
    yaw_t_down = squeeze(mean(reshape(bdata_vel_time, [DT_YAW, length(bdata_vel_time)/DT_YAW]),1));
    yaw_down = squeeze(mean(reshape(cur_yaw, [DT_YAW, length(bdata_vel_time)/DT_YAW]),1));    

    hold on;
    
    SHIFT_FACTOR = 3;
    for ii = 1:(length(A2_Vm_down)-SHIFT_FACTOR)
        cur_yaw_1 = yaw_down( ii+SHIFT_FACTOR );
        plot(A2_Vm_down(ii), cur_yaw_1, 'o', 'MarkerSize', 3, 'color', 'b' );
    end
    
    for jj = 1:length( counter_turn_epochs )
        cur_tr = counter_turn_epochs{jj}{1};
        if( cur_tr == tr )
            cur_start_t = counter_turn_epochs{jj}{2}{1};
            cur_end_t = counter_turn_epochs{jj}{2}{2};
            
            cur_t_down = find( ( t_down >= cur_start_t ) & ( t_down <= cur_end_t ) );
            cur_yaw_t_down = find( ( yaw_t_down >= cur_start_t ) & ( yaw_t_down <= cur_end_t ) );
            
            for ii = 1:(length(cur_t_down)-SHIFT_FACTOR)
                cur_index_yaw   = cur_yaw_t_down(ii);
                cur_index_ephys = cur_t_down(ii);
                cur_yaw_2 = yaw_down( cur_index_yaw +SHIFT_FACTOR );
                plot(A2_Vm_down(cur_index_ephys), cur_yaw_2, 'o', 'MarkerSize', 10, 'color', 'r' );
            end            
        end    
    end
end

xlabel('A2 Vm (mV)');
ylabel('Yaw (au)');

saveas(f,[analysis_path '/Vm_vs_yaw_spontaneous_CX_turning_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/Vm_vs_yaw_spontaneous_CX_turning_' num2str(sid) '.png']);

%%
traj = get_single_trial_trajectories( sid, bdata_vel_time, bdata_vel );

num_trials = size(traj{1},1);

%TRIAL_MAX = 6;
TRIAL_MAX = num_trials;

f = figure;
cm = colormap(jet(TRIAL_MAX));
hold on;
for tr = 1:TRIAL_MAX
    disp_x = squeeze(traj{1}(tr,1,:));
    disp_y = squeeze(traj{1}(tr,2,:));
    plot(disp_x, disp_y, 'color', cm(tr,:), 'DisplayName', ['trial: ' num2str(tr)]);
end

legend();
xlabel('X displacement (au)');
ylabel('Y displacement (au)');

saveas(f,[analysis_path '/running_trajectories_' num2str(sid) '_max_trials_' num2str( TRIAL_MAX ) '.fig']);
saveas(f,[analysis_path '/running_trajectories_' num2str(sid) '_max_trials_' num2str( TRIAL_MAX ) '.png']);

%% Exploring large right turns

f = figure;
tr = 11;
disp_x = squeeze(traj{1}(tr,1,:));
disp_y = squeeze(traj{1}(tr,2,:));
plot3( disp_x, disp_y, bdata_vel_time, 'color', cm(tr,:), 'DisplayName', ['trial: ' num2str(tr)]);
xlabel('Dist X');
ylabel('Dist Y');
zlabel('Time (s)');
legend();

saveas(f,[analysis_path '/running_trajectories_' num2str(sid) '_trial_' num2str( tr ) '.fig']);
saveas(f,[analysis_path '/running_trajectories_' num2str(sid) '_trial_' num2str( tr ) '.png']);

%% 
figure;
for f=1:size( EPG_data, 1)

    imagesc(squeeze(mean(squeeze(EPG_data(f,:,:,:)),3)));
    
    waitforbuttonpress;
end


%% Plot avg behavior
ac = get_analysis_constants;

first_stim = 3.0;
last_stim = 3.5;

f = figure;

avg_yaw = squeeze(mean(squeeze(bdata_vel{ 1 }( :, ac.VEL_YAW, : ))));
avg_fwd = squeeze(mean(squeeze(bdata_vel{ 1 }( :, ac.VEL_FWD, : ))));

ax_1(1) = subplot(2,1,1); 
hold on;

p = plot( bdata_vel_time, avg_fwd );
xlim([0 bdata_vel_time(end)]);
ylabel('Fwd vel (au/s)');
axis tight;
legend(p,['n = (' num2str( size(bdata_vel{1},1) ) ')']);

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ first_stim first_stim last_stim last_stim ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

ax_1(2) = subplot(2,1,2);
hold on;
plot( bdata_vel_time, avg_yaw );
xlim([0 bdata_vel_time(end)]);
xlabel('Time (s)');
ylabel('Yaw (au/s)');
axis tight;

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ first_stim first_stim last_stim last_stim ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');


linkaxes(ax_1,'x');


saveas(f,[analysis_path '/fwd_yaw_tc_sid_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/fwd_yaw_tc_sid_' num2str(sid) '.png']);








