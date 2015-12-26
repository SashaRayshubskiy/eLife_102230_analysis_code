%% Load imaging and behavioral data
global slash;

clear all;

if isunix() == 1
    slash = '/';
else
    slash = '\';
end

% Must end with a slash
datapath = '/data/drive_fast/sasha/151224_bfly_2/';

analysis_path = [datapath slash 'analysis'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

sid = 7;

% Load behavioral data
bdata_path = [datapath  slash 'ball' slash ];
tic; [ b_rawdata, b_time, btrial_meta ] = load_natural_odor_daq_data(sid, bdata_path); toc

% Load imaging data
%cdata_path = [datapath  slash '2p' slash ];
%tic; [ cdata_raw, cdata_meta, ctrial_meta ] = load_imaging_data(sid, cdata_path); toc

% Check that the behavioral and imaging trials match up
%check_bdata_and_cdata_trial_integrity( btrial_meta, ctrial_meta );

% Get behavioral data that is usable for analysis

%[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data( b_time, b_rawdata );

settings = sensor_settings;

%% Display behavioral data
display_avg_velocity(b_rawdata, bdata_vel, bdata_vel_time, analysis_path);

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

frame_clocks_for_all = vertcat( squeeze(b_rawdata(:,:,5)) );

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

saveas(f, [analysis_path '/image_acq_variability_externally_triggered.fig']);




