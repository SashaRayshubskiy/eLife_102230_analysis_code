%% Load imaging and behavioral data
global slash;

clear all;

if isunix() == 1
    slash = '/';
else
    slash = '\';
end

% Must end with a slash
datapath = '/data/drive_fast/sasha/151226_z_drift_1/';

analysis_path = [datapath slash 'analysis'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

sid = 2;

trial_type_cnt = 1;

% Load behavioral data
bdata_path = [datapath  slash 'ball' slash ];
tic; [ b_rawdata, b_time, btrial_meta ] = load_behavioral_data(sid, bdata_path, trial_type_cnt ); toc

% Load imaging data
cdata_path = [datapath  slash '2p' slash ];
tic; [ cdata_raw, cdata_meta, ctrial_meta ] = load_imaging_data(sid, cdata_path, trial_type_cnt ); toc

% Check that the behavioral and imaging trials match up
check_bdata_and_cdata_trial_integrity( btrial_meta, ctrial_meta );

settings = sensor_settings;

%% Display behavioral data
display_avg_velocity(b_rawdata, bdata_vel, bdata_vel_time, analysis_path);

%% Play with stim params
one_trial_bdata = squeeze(b_rawdata(1,:,:));

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

one_trial_bdata = squeeze(b_rawdata{1}(1,:,:));
frame_clock     = squeeze(one_trial_bdata(:,5));

figure;
hold on;
plot(b_time, frame_clock);

first_stim = find((frame_clock > 2.0), 1, 'first');
last_stim = find((frame_clock > 2.0), 1, 'last');

plot(b_time(first_stim), 5.0, 'xr');
plot(b_time(last_stim), 5.0, 'xr');

%% Analyze frame start and duration times

frame_clocks_for_all = vertcat( squeeze(b_rawdata{1}(:,:,5)) );

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

saveas(f, [analysis_path '/sid_' num2str(sid) '_image_acq_variability_externally_triggered.fig']);

%% Determine z-level for a single image
single_img = squeeze(cdata_raw{1}(1,:,:,1,2,10));

figure;
%imagesc(mean(single_img));
%axis image;
%colormap gray;
hold on;
x = [1:size(single_img,2)];
y = squeeze(mean(single_img,2));

% plot(x,y);
f = fit(x',y,'gauss2');
plot(f,x,y);

[max_val,coord] = max(y);

%% Try parfor
pobj = parpool(20);

%% Characterize z-drift
PLANES = size(cdata_raw{1}, 5);
VOLUMES = size(cdata_raw{1}, 6);
NUM_TRIALS = size(cdata_raw{1}, 1);

z_level = zeros( NUM_TRIALS, PLANES, VOLUMES );

parfor_progress(NUM_TRIALS);

for tr = 1:NUM_TRIALS
    
    cur_trial_data = squeeze( cdata_raw{1}(tr,:,:,1,:,:));
    
    
    z_levels_p_v = get_z_levels(cur_trial_data, size(cdata_raw{1},3), PLANES, VOLUMES );
        
    z_level(tr,:,:) = z_levels_p_v;
    
    parfor_progress;
end
parfor_progress(0);

%% 
f = figure;
subplot(1,1,1);
surf(squeeze(z_level(:,:,1)));
ylabel('Trial #');
xlabel('Plane #');
zlabel('Z level');

saveas(f, [analysis_path '/sid_' num2str(sid) '_z_level_trials_vs_plane.fig' ]);

f = figure;
for p=1:PLANES
    subplot(4,4,p);
    surf(squeeze(z_level(:,p,:)));
    ylabel('Trial #');
    xlabel('Volume #');
    zlabel('Z level');
    zlim([0 512]);
    xlim([0 VOLUMES]);
    ylim([0 NUM_TRIALS]);
end

saveas(f, [analysis_path '/sid_' num2str(sid) '_z_level_trials_vs_volumes.fig' ]);

%% Z delta histogram

avg_plane_z = squeeze(mean(mean(z_level(:,2:end, :), 3)));

f = figure;
subplot(2,1,1);
plot([2:16], avg_plane_z);

xlabel('Plane #');
ylabel('Z level');

subplot(2,1,2);
%plot([3:16], diff(avg_plane_z));
hist(-1.0*diff(avg_plane_z))
saveas(f, [analysis_path '/sid_' num2str(sid) '_z_level_plane_hist.fig' ]);

%% Z delta histogram

all_z_planes = reshape(permute(z_level, [2 3 1]), [1, size(z_level,1)*size(z_level,2)*size(z_level,3)]);

f = figure;
subplot(2,1,1);
plot(all_z_planes);
xlabel('Frame #');
ylabel('Z level');

subplot(2,1,2);
histogram(-1.0*diff(all_z_planes), 'BinWidth', 0.1);
xlabel('Delta z level');
ylabel('Count (frames)');
saveas(f, [analysis_path '/sid_' num2str(sid) '_z_level_plane_hist_all.fig' ]);