%% Load imaging and behavioral data
global slash;

clear all;

if isunix() == 1
    slash = '/';
else
    slash = '\';
end

% Experiment

% sid 0, 1 -- very low number of trials, no clear bump returns
% datapath = '/data/drive2/sasha/181022_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15/';

% sid 0, 2
% datapath = '/data/drive2/sasha/181203_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16/'; 

% sid 0
datapath = '/data/drive2/sasha/181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17/'; 

% sid 1
% datapath = '/data/drive2/sasha/181211_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18/'; 

% sid 0
% datapath = '/data/drive2/sasha/190131_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19/'; 

% sid 0
% datapath = '/data/drive2/sasha/190204_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22/'; 

% sid 1 
% datapath = '/data/drive2/sasha/190205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23/'; 


% Control

% sid 0
datapath = '/data/drive2/sasha/181206_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_01/';

% sid 0
% datapath = '/data/drive2/sasha/190208_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_02/';

% sid 1
% datapath = '/data/drive2/sasha/190211_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_03/';

% sid 0
% datapath = '/data/drive2/sasha/190212_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_04/';

% sid 0 -- poor bump data
% datapath = '/data/drive2/sasha/190213_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_05/';

% sid 1 -- bad recording
% datapath =  '/data/drive2/sasha/181211_Lex_6f_60D05_Gal4_P2X2_control_02';

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
% Drop the first 4 planes..
% trials, x, y, color, plane, time
% START_PLANE = 4; % For 9 planes
START_PLANE = 5; % For 16 planes

EPG_data_1       = squeeze(cdata_raw{1}( :, :, :, 1, START_PLANE:end, : ));
% Get max for each plane
EPG_data = squeeze(max(EPG_data_1, [], 4));

ephys_data_delta_Vm = zeros( size(ephys_data{1}) );

for tr = 1:size(ephys_data{1},1)
    
    cur_ephys = ephys_data{1}(tr,:);
    
    ephys_data_delta_Vm(tr,:) = cur_ephys - mean(cur_ephys);
end

%% Display EPG bump, yaw, and ephys on the same plot for quick view of experiment
[df_f_in_roi_per_trial] = display_EB_dynamics_var_stim( EPG_data_1, pico_stim_data, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid );

%% Plot avg yaw and ephys, triggered on stim onset
EBYE = EB_yaw_ephys_data; 

% Experiment
% stims_to_include = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15_181022_sid_1.stims_to_include;

% stims_to_include = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16_181203.stims_to_include;
% stims_to_include = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17_181205.stims_to_include;
% stims_to_include = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18_181211.stims_to_include;

stims_to_include = [];

% Control 
% stims_to_include = EBYE.Lex_6f_60D05_Gal4_P2X2_control_01_181206.stims_to_include;

[ trial_stim_id_map ] = display_avg_yaw_ephys_triggered_on_stim_onset( df_f_in_roi_per_trial, pico_stim_data, ephys_time, ephys_data_delta_Vm, bdata_vel_time, bdata_vel, VPS, analysis_path, sid, stims_to_include );
%[stim_events, trial_stim_id_map] = display_avg_yaw_ephys_triggered_on_stim_onset_left_right_bump( df_f_in_roi_per_trial, pico_stim_data, ephys_time, ephys_data_delta_Vm, bdata_vel_time, bdata_vel, VPS, analysis_path, sid, stims_to_include );



%% Show bump analysis 
pad = '/data/drive2/sasha/';
EBYE = EB_yaw_ephys_data; 

% Path to bump_yaw_ephys_in_window_data_sid_<>.mat file
% exp_dirs = { {[ pad '181022_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15_181022_sid_1.stims_to_include } };
%              {[ pad '181022_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15/analysis/'], 1, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15_181022_sid_1.stims_to_include }, ...
%              {[ pad '181203_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16_181203.stims_to_include }, ...
%              {[ pad '181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17_181205.stims_to_include }, ...
%              {[ pad '181211_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18/analysis/'], 1, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18_181211.stims_to_include } };

% exp_dirs = { {[ pad '181203_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16_181203.stims_to_include } };
% exp_dirs = { {[ pad '181022_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15_181022.stims_to_include } };
% exp_dirs = { {[ pad '181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17_181205.stims_to_include } };
% exp_dirs = { {[ pad '181211_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18/analysis/'], 1, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18_181211.stims_to_include } };
% exp_dirs = { {[ pad '190131_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19_190131.stims_to_include } };

% exp_dirs = { {[ pad '190204_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22_190204.stims_to_include } };
exp_dirs = { {[ pad '190205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23/analysis/'], 1, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23_190205.stims_to_include } };

% analysis_all_path = [pad '/EB_bump_yaw_ephys_all/'];
% if(~exist(analysis_all_path, 'dir'))
%     mkdir(analysis_all_path);
% end   

%% Run this to obtain bump velocity plots, used for bump return alignment
analysis_of_bump_speed_vs_yaw_v4( exp_dirs, VPS, analysis_path );


%% Show bump return analysis
clear EBYE;
EBYE = EB_yaw_ephys_data; 
% cur_bump_return = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16_181203_cbr_down;
% cur_bump_return = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17_181205_cbr_down;
% cur_bump_return = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18_181211_cbr_down;

% cur_bump_return = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19_190131_cbr_up; direction_of_return = 'up';
% cur_bump_return = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19_190131_cbr_down; direction_of_return = 'down';

% cur_bump_return = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22_190204_cbr_up; direction_of_return = 'up';
% cur_bump_return = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22_190204_cbr_down; direction_of_return = 'down';

% cur_bump_return = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23_190205_cbr_up; direction_of_return = 'up';
cur_bump_return = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23_190205_cbr_down; direction_of_return = 'down';

analysis_of_bump_return( exp_dirs, cur_bump_return, VPS, analysis_path, direction_of_return );

%% Analysis of bump return for all flies
pad = '/data/drive2/sasha/';
EBYE = EB_yaw_ephys_data; 

VPS = 12.4;

direction_of_return = 'down';
cur_all_fly_dirs_bump_down = { {[ pad '181203_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16_181203.stims_to_include, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16_181203_cbr_down }, ...
                           {[ pad '181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17_181205.stims_to_include, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17_181205_cbr_down }, ...
                           {[ pad '181211_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18/analysis/'], 1, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18_181211.stims_to_include, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18_181211_cbr_down }, ...
                           {[ pad '190131_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19_190131.stims_to_include, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19_190131_cbr_down }, ...
                           {[ pad '190204_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22_190204.stims_to_include, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22_190204_cbr_down }, ...
                           {[ pad '190205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23/analysis/'], 1, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23_190205.stims_to_include, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23_190205_cbr_down } };

[t_win_down, cx_vars_down] = analysis_of_bump_return_v2( cur_all_fly_dirs_bump_down, VPS, direction_of_return );

direction_of_return = 'up';
cur_all_fly_dirs_bump_up = { {[ pad '181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17_181205.stims_to_include, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17_181205_cbr_up }, ...
                             {[ pad '190131_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19_190131.stims_to_include, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19_190131_cbr_up }, ...
                             {[ pad '190204_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22/analysis/'], 0, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22_190204.stims_to_include, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22_190204_cbr_up }, ...
                             {[ pad '190205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23/analysis/'], 1, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23_190205.stims_to_include, EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23_190205_cbr_up } };

[t_win_up, cx_vars_up] = analysis_of_bump_return_v2( cur_all_fly_dirs_bump_up, VPS, direction_of_return );


%% Plots both up and down CX vars
cx_file_path = [pad 'CX_summary/'];

f = figure;

down_clr_avg = rgb('SeaGreen');
up_clr_avg = rgb('DarkMagenta');

ax1(1) = subplot(3,1,1);
hold on;
avg_bump_down = mean(cx_vars_down{1});
sem_bump_down = get_sem(cx_vars_down{1}, 1);

avg_bump_up = mean(cx_vars_up{1});
sem_bump_up = get_sem(cx_vars_up{1},1);

fh = fill( [t_win_down{1}, fliplr(t_win_down{1})], ...
    [(avg_bump_down-sem_bump_down) fliplr((avg_bump_down+sem_bump_down))], ...
    rgb('PaleGreen'));
set(fh, 'EdgeColor', 'None');

fh = fill( [t_win_up{1}, fliplr(t_win_up{1})], ...
    [(avg_bump_up-sem_bump_up) fliplr((avg_bump_up+sem_bump_up))], ...
    rgb('Violet'));
set(fh, 'EdgeColor', 'None');

plot( t_win_down{1}, avg_bump_down, 'color', down_clr_avg, 'LineWidth', 2 );
plot( t_win_up{1}, avg_bump_up, 'color', up_clr_avg, 'LineWidth', 2 );

ylabel('EB vel(au/s)');

ax1(2) = subplot(3,1,2);
hold on;
avg_yaw_down = mean(cx_vars_down{2});
sem_yaw_down = get_sem(cx_vars_down{2}, 1);

avg_yaw_up = mean(cx_vars_up{2});
sem_yaw_up = get_sem(cx_vars_up{2},1);

fh = fill( [t_win_down{2}, fliplr(t_win_down{2})], ...
    [(avg_yaw_down-sem_yaw_down) fliplr((avg_yaw_down+sem_yaw_down))], ...
    rgb('PaleGreen'));
set(fh, 'EdgeColor', 'None');

fh = fill( [t_win_up{2}, fliplr(t_win_up{2})], ...
    [(avg_yaw_up-sem_yaw_up) fliplr((avg_yaw_up+sem_yaw_up))], ...
    rgb('Violet'));
set(fh, 'EdgeColor', 'None');

plot( t_win_down{2}, avg_yaw_down, 'color', down_clr_avg, 'LineWidth', 2 );
plot( t_win_up{2}, avg_yaw_up, 'color', up_clr_avg, 'LineWidth', 2 );

ylabel('Yaw vel (au/s)');

ax1(3) = subplot(3,1,3);
hold on;
avg_ephys_down = mean(cx_vars_down{3});
sem_ephys_down = get_sem(cx_vars_down{3}, 1);

avg_ephys_up = mean(cx_vars_up{3});
sem_ephys_up = get_sem(cx_vars_down{3},1);

fh = fill( [t_win_down{3}, fliplr(t_win_down{3})], ...
    [(avg_ephys_down-sem_ephys_down) fliplr((avg_ephys_down+sem_ephys_down))], ...
    rgb('PaleGreen'));
set(fh, 'EdgeColor', 'None');

fh = fill( [t_win_up{3}, fliplr(t_win_up{3})], ...
    [(avg_ephys_up-sem_ephys_up) fliplr((avg_ephys_up+sem_ephys_up))], ...
    rgb('Violet'));
set(fh, 'EdgeColor', 'None');

plot( t_win_down{3}, avg_ephys_down, 'color', down_clr_avg, 'LineWidth', 2 );
plot( t_win_up{3}, avg_ephys_up, 'color', up_clr_avg, 'LineWidth', 2 );

ylabel('Vm (mV)');
xlabel('Time (s)');

linkaxes(ax1, 'x');
saveas(f, [ cx_file_path 'bump_return_yaw_ephys_all_flies.fig' ]);
saveas(f, [ cx_file_path 'bump_return_yaw_ephys_all_flies.png' ]);

%%

trial_stim_id_map = [];
show_running_trajectories_for_good_bump_jumps(  sid, bdata_vel_time, bdata_vel, trial_stim_id_map, pico_stim_data, analysis_path );

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


%%%% ATTIC 
%% 
display_bump_dynamics_turning_ephys( df_f_in_roi_per_trial, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid, pre_stim, stim_dur );

%% Bump dynamics with yaw and ephys with variable stim
display_bump_dynamics_turning_ephys_var_stim( df_f_in_roi_per_trial, pico_stim_data, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid );


%% Analyze the data saved for time warping
warp_EB_data_upsample( stim_events, VPS, analysis_path );


