%% Load imaging and behavioral data
global slash;

clear all;

if isunix() == 1
    slash = '/';
else
    slash = '\';
end

% Must end with a slash
%datapath = '/data/drive_fast/sasha/160118_R84C10_83blexA_02/';
datapath = '/data/drive_fast/sasha/160122_nsyb_83blexA_01/';

analysis_path = [datapath slash 'analysis'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

sid = 0;

trial_type_cnt = 3;

% Load behavioral data
bdata_path = [datapath  slash 'ball' slash ];
tic; [ b_rawdata, b_time, btrial_meta ] = load_behavioral_data(sid, bdata_path, trial_type_cnt ); toc

% Load imaging data
cdata_path = [datapath  slash '2p' slash ];
tic; [ cdata_raw, cdata_meta, ctrial_meta ] = load_imaging_data(sid, cdata_path, trial_type_cnt ); toc

% Check that the behavioral and imaging trials match up
check_bdata_and_cdata_trial_integrity( btrial_meta, ctrial_meta );

% Get behavioral data that is usable for analysis

[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data( b_time, b_rawdata );

settings = sensor_settings;
global file_writer_cnt;
file_writer_cnt = 1;

%% Save imaging data for faster load
save

%% Display individual trials. Behavior along with gcamp.
a_const = get_analysis_constants;

cur_trial_type = a_const.RIGHT;
cur_trial_type_str = a_const.task_str{cur_trial_type};
trial_id = 4;

cur_cdata     = squeeze(cdata_raw{ cur_trial_type }(trial_id,:,:,:,:,:));
cur_bdata_vel = squeeze(bdata_vel{ cur_trial_type }(trial_id,:,:));

figsave_prefix = [analysis_path '/clicky_with_behaviour_' cur_trial_type_str '_tid_' num2str(trial_id) ];
display_avg_volume(cur_cdata, figsave_prefix);

cur_plane = 8;  

cur_plane_cdata = squeeze(cur_cdata(:,:,cur_plane,:));

figsave_prefix = [analysis_path '/clicky_with_behaviour_' cur_trial_type_str '_tid_' num2str(trial_id) '_plane_' num2str(cur_plane) ];

VPS = cdata_meta.volume_rate;
clicky_with_behaviour( cur_plane_cdata, cur_bdata_vel, bdata_vel_time, VPS, settings, figsave_prefix );

%% Display behavioral data
display_avg_velocity(b_rawdata, bdata_vel, bdata_vel_time, analysis_path);

%% Play movie on a trial
a_const = get_analysis_constants;

cur_trial_type = a_const.RIGHT;
cur_trial_type_str = a_const.task_str{cur_trial_type};
trial_id = 4;

cur_cdata     = squeeze(cdata_raw{ cur_trial_type }(trial_id,:,:,:,:,:));
VPS = cdata_meta.volume_rate;

figsave_prefix = [analysis_path '/mov_' cur_trial_type_str '_tid_' num2str(trial_id) ];
play_volumetric_movie(cur_cdata, VPS, figsave_prefix);

%% Play with stim params
one_trial_bdata = squeeze(b_rawdata{ BOTH }(1,:,:));

right_odor_stim = squeeze(one_trial_bdata(:,6));
left_odor_stim  = squeeze(one_trial_bdata(:,7));

figure;
hold on;
plot(b_time, left_odor_stim);
first_stim = find((left_odor_stim > 2.0), 1, 'first');
last_stim = find((left_odor_stim > 2.0), 1, 'last');
plot(b_time(first_stim), 5.0, 'xr');
plot(b_time(last_stim), 5.0, 'xr');

%% Play with frame clock

one_trial_bdata = squeeze(b_rawdata{ BOTH }(1,:,:));
frame_clock     = squeeze(one_trial_bdata(:,5));

figure;
hold on;
plot(b_time, frame_clock);

first_stim = find((frame_clock > 2.0), 1, 'first');
last_stim = find((frame_clock > 2.0), 1, 'last');

plot(b_time(first_stim), 5.0, 'xr');
plot(b_time(last_stim), 5.0, 'xr');

%% Analyze frame start and duration times

frame_clocks_for_all = vertcat( squeeze(b_rawdata{ BOTH }(:,:,5)), squeeze(b_rawdata{ LEFT }(:,:,5)), squeeze(b_rawdata{ RIGHT }(:,:,5)) );

for i=1:size(frame_clocks_for_all,1)
    cur_frame_clock = squeeze(frame_clocks_for_all(i,:));
    frame_start = find((cur_frame_clock > 2.0), 1, 'first');
    frame_end = find((cur_frame_clock > 2.0), 1, 'last');
    frame_duration = frame_end - frame_start;
    
    frame_data(i,:) = [frame_start, frame_duration];
end

f = figure;
subplot(2,1,1);
hist(frame_data(:,1)./settings.sampRate);
xlabel('Time (s)');
ylabel('Counts (trials)');
title('Frame start variability');

subplot(2,1,2);
hist(frame_data(:,2)./settings.sampRate);
xlabel('Time (s)');
ylabel('Counts (trials)');
title('Imaging run duration variability');

saveas(f, [analysis_path '/image_acq_variability.fig']);

%% 

