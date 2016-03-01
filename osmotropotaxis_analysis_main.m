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
%datapath = '/data/drive_fast/sasha/160125_nsyb_83blexA_02/';
%datapath = '/data/drive_fast/sasha/160122_nsyb_83blexA_01/';
%disp('CAUTION: Using blank trials from nsyb_83blexA_01');
%nsyb_83blexA_01_blank_trials = { [165], [163, 164, 167, 359], [166, 360] };
%nsyb_83blexA_06_blank_trials = { [44,342], [], [491] };

%trial_exclusion_list = nsyb_83blexA_01_blank_trials;
trial_exclusion_list = {[],[],[]};

datapath = '/data/drive0/sasha/160229_nsyb_83blexA_27/';

analysis_path = [datapath slash 'analysis'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

sid = [0];

aconstants = get_analysis_constants;
trial_type_cnt = 2;

% Load behavioral data
bdata_path = [datapath  slash 'ball' slash ];
tic; [ b_rawdata, b_time, btrial_meta ] = load_behavioral_data(sid, bdata_path, trial_type_cnt ); toc

% Load imaging data
cdata_path = [datapath  slash '2p' slash ];
dx = 1;
dy = 1;
dt = 2;
down_params(1) = dx;
down_params(2) = dy;
tic; [ cdata_raw, cdata_meta, ctrial_meta ] = load_imaging_data(sid, cdata_path, trial_type_cnt, dx, dy ); toc

% Check that the behavioral and imaging trials match up
check_bdata_and_cdata_trial_integrity( btrial_meta, ctrial_meta );

% Get behavioral data that is usable for analysis

[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data( b_time, b_rawdata );

settings = sensor_settings;
global file_writer_cnt;
file_writer_cnt = 1;

% Create the reverse mapping from 
external_trial_id_to_internal_ordinal_map = get_external_trial_id_to_internal_ordinal_map(btrial_meta);

% Create a map of frame start offsets per plane. (This is key because
% there's about 150 ms interval from first plane to first plane in a
% volume.
planes = size( cdata_raw{ 1 }, 5 );
VPS = cdata_meta.volume_rate;
frame_start_offsets_per_plane = generate_frame_start_offsets_per_plane( planes, b_rawdata, b_time );

ref_imgs = generate_ref_imgs(cdata_raw);

%% Display behavioral data
display_avg_velocity(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);

%% Display behavioral data
%datapath_tmp = '/data/drive0/sasha/160211_nsyb_83blexA_11/';
datapath_tmp = '/data/drive0/sasha/160229_nsyb_83blexA_27/';

analysis_path_tmp = [datapath_tmp slash 'analysis'];

if(~exist(analysis_path_tmp, 'dir'))
    mkdir(analysis_path_tmp);
end

sid_tmp = [0];
trial_type_cnt_tmp = 2;

bdata_path_tmp = [datapath_tmp  slash 'ball' slash ];
tic; [ b_rawdata_tmp, b_time_tmp, btrial_meta_tmp ] = load_behavioral_data(sid_tmp, bdata_path_tmp, trial_type_cnt_tmp ); toc

[bdata_vel_time_tmp, bdata_vel_tmp] = reformat_raw_behavioral_data( b_time_tmp, b_rawdata_tmp );

% display_avg_velocity(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);
% display_avg_velocity_exclude_zero_vel(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);

display_avg_velocity_exclude_zero_vel_RL_only(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);
display_avg_velocity_left_right_only(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);
display_per_trial_velocity(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);

%%
display_avg_velocity_exclude_zero_vel_RL_only(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);
display_avg_velocity_left_right_only(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);
display_per_trial_velocity(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);

%% Generate turning clusters using standard unsupervised clustering tools.
MAX_CLUSTERS = 5;
clust = generate_turning_clusters( sid, bdata_vel_time, bdata_vel, analysis_path, MAX_CLUSTERS );


%% Test the idea that activation of 'LAL' turning region between the 2 sides corresponds to turn magnitude.

generate_turning_magnitude_vs_bilateral_calcium_delta_response_plot( bdata_vel_time, bdata_vel, cdata_raw, frame_start_offsets_per_plane, VPS, analysis_path );


%% Generate expected vs. ignored

turn_metadata = generate_turning_metadata( sid, bdata_vel_time, bdata_vel, analysis_path );

% [ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list( bdata_vel_time, bdata_vel );
[ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list_v2( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path );

%% Generate large vs. small counter turn

turn_metadata = generate_turning_metadata( sid, bdata_vel_time, bdata_vel, analysis_path );

[ condition_trials, condition_trials_str, condition_str ] = generate_large_vs_small_turn_trial_list( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path )

%% Generate large vs. small counter turn
[ condition_trials, condition_trials_str, condition_str ] = generate_large_vs_small_counter_turn_trial_list( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path )

%% Generate counter turn vs. no counter turn
[ condition_trials, condition_trials_str, condition_str ] = generate_counter_turn_vs_no_counter_turn_trial_list( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path )

%% Generate early vs. late turns
[ condition_trials, condition_trials_str, condition_str ] = generate_early_vs_late_turn_trial_list( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path )

%% Generate quick vs. delayed counter turns
[ condition_trials, condition_trials_str, condition_str ] = generate_quick_vs_delayed_counter_turn_trial_list( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path )

%% Generate stationary vs. motion

[ condition_trials, condition_trials_str, condition_str ] = gen_stm_vs_ft_motion_trials_160122_nsyb_83blexA_01( bdata_vel_time, bdata_vel, external_trial_id_to_internal_ordinal_map );

%% Display behavioral 2 condition trials.

%condition_trials{2}(2) = 14;

avg_cond_btrace_trace_filepath = [ analysis_path '/' condition_str '_asid_' num2str( asid ) '_sid_' num2str(sid) ];
with_single_trials = 1;
display_two_condition_trials( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, with_single_trials );

with_single_trials = 0;
display_two_condition_trials( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, with_single_trials );

%% Get covariance plots from 2 symmetric ROIs
ac = get_analysis_constants();

PLANE_OF_INTEREST = 9;

roi_session = 0;
rois_path = [analysis_path '/roi_session_' num2str(roi_session) '.mat' ];
rois_1 = get_rois_from_volume( PLANE_OF_INTEREST, squeeze(cdata_raw{ 1 }(1,:,:,:,:,:)), rois_path );

generate_two_rois_comparison_with_turning_metadata( sid, rois_1, PLANE_OF_INTEREST, cdata_raw, turn_metadata, frame_start_offsets_per_plane, VPS, btrial_meta, analysis_path );

%% Create a differece image for each plane
tic; [ btraces_per_condition, avg_df_f_per_condition_per_plane ] = collect_two_behavioral_condition_and_df_f_per_cond( condition_trials, cdata_raw, bdata_vel, VPS, trial_exclusion_list, btrial_meta ); toc;

%%
a_const = get_analysis_constants;

cur_trial_type = a_const.RIGHT;
cur_trial_type_str = a_const.task_str{cur_trial_type};
trial_id = 4;

cur_cdata     = squeeze(cdata_raw{ cur_trial_type }(trial_id,:,:,:,:,:));
figsave_prefix = [analysis_path '/avg_volume_' cur_trial_type_str '_tid_' num2str(trial_id) ];
display_avg_volume( down_params, cur_cdata, figsave_prefix );

%%
roi_session = 30;

diff_image_path = [ analysis_path '/' condition_str '_diff_image_asid_' num2str( asid ) '_sid_' num2str(sid) '_roi_session_' num2str(roi_session) ];
display_two_condition_difference_image(down_params, ref_imgs, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% Clicky showing avg data for 2 conditionsquick

ac = get_analysis_constants();

roi_session = 19;
PLANE_OF_INTEREST = 9;
TRIAL_TYPE_OF_INTEREST = ac.RIGHT;
ref_img = mean(squeeze(cdata_raw{ 1 }(1,:,:,PLANE_OF_INTEREST,:)),3);

diff_image_path = [ analysis_path '/' condition_str '_diff_image_asid_' num2str( asid ) '_sid_' num2str(sid) '_roi_session_' num2str(roi_session) ];
clicky_two_condition_bdata(ref_img, PLANE_OF_INTEREST, TRIAL_TYPE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% Testing new unsupervised comparison of clustered data
roi_session = 0;

PLANE_OF_INTEREST = 10;
TRIAL_TYPE_OF_INTEREST = 2;
ref_img = mean(squeeze(cdata_raw{ 1 }(1,:,:,PLANE_OF_INTEREST,:)),3);

diff_image_path = [ analysis_path '/' condition_str '_diff_image_asid_' num2str( asid ) '_sid_' num2str(sid) '_roi_session_' num2str(roi_session) ];
display_two_condition_difference_image_debug(ref_img, PLANE_OF_INTEREST, TRIAL_TYPE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% test dithering filter
img_load = load('/tmp/diff_img.mat');
img = img_load.diff_img;
img_filt = filter_dithered_image( img );

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
imagesc(img);
axis image;
colormap jet;
caxis([-0.1 0.5]);
title('Before');

subplot(2,1,2);
imagesc(img_filt);
axis image;trial_exclusion_list = {[],[],[]};
colormap jet;
caxis([-0.1 1.0]);
title('After');

%% Collect and display time courses in an ROI, both condition        
rois = get_rois_from_volume_v2( asid, squeeze(cdata_raw{ 1 }(1,:,:,:,:,:)), analysis_path );

tic; [ btraces_per_condition, ctraces_in_roi_per_condition ] = collect_two_behavioral_condition_traces( condition_trials, cdata_raw, bdata_vel, VPS, rois, trial_exclusion_list, btrial_meta ); toc;

avg_trace_filepath = [ analysis_path '/' condition_str '_avg_traces_asid_' num2str( asid ) '_sid_' num2str(sid) ];
display_two_behavioral_condition_traces( condition_trials_str, btraces_per_condition, ctraces_in_roi_per_condition, bdata_vel_time, frame_start_offsets_per_plane, VPS, avg_trace_filepath );

diff_avg_trace_filepath = [ analysis_path '/' condition_str '_avg_diff_traces_asid_' num2str( asid ) '_sid_' num2str(sid) ];
display_two_behavioral_condition_diff_traces( condition_trials_str, btraces_per_condition, ctraces_in_roi_per_condition, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_avg_trace_filepath );
[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data( b_time, b_rawdata );

settings = sensor_settings;
global file_writer_cnt;
file_writer_cnt = 1;

% Create the reverse mapping from 
external_trial_id_to_internal_ordinal_map = get_external_trial_id_to_internal_ordinal_map(btrial_meta);

% Create a map of frame start offsets per plane. (This is key because
% there's about 150 ms interval from first plane to first plane in a
% volume.
planes = size( cdata_raw{ 1 }, 5 );
VPS = cdata_meta.volume_rate;
frame_start_offsets_per_plane = generate_frame_start_offsets_per_plane( planes, b_rawdata, b_time );

ref_imgs = generate_ref_imgs(cdata_raw);
%% Generate time courses for each plane, using the ROIs of a session

cur_bdata_vel = squeeze(bdata_vel{ cur_trial_type }( trial_ord,:,: ));

cur_tbt_filename_prefix = [ analysis_path '/time_courses_in_volume_asid_' num2str(asid) '_sid_' num2str(sid) '_' cur_trial_type_str '_tid_' num2str(cur_trial_id)];
generate_volume_time_courses_with_rois( cur_cdata, cur_bdata_vel, bdata_vel_time, new_rois, frame_start_offsets_per_plane, VPS, cur_tbt_filename_prefix );

%% Add an ROI to a plane
plane_to_append = 5;
new_rois = add_rois_from_volume(cur_cdata, rois_v2, plane_to_append);

%%  Generate ROIs 
asid = 0; % roi analysis session id

rois = get_rois_from_volume_v2( asid, squeeze(cdata_raw{ 1 }(1,:,:,:,:,:)), analysis_path );

a_const = get_analysis_constants;
cur_trial_type = a_const.RIGHT;
trial_ord = 4;

cur_cdata     = squeeze(cdata_raw{ cur_trial_type }(trial_ord,:,:,:,:,:));
cur_trial_id = squeeze(btrial_meta{ cur_trial_type }(trial_ord, 2));

cur_tbt_filename_prefix = [ analysis_path '/roi_avg_volume_asid_' num2str(asid) '_sid_' num2str(sid) '_' cur_trial_type_str '_tid_' num2str(cur_trial_id)];
generate_volume_avg_with_rois( cur_cdata, rois, cur_tbt_filename_prefix );

%% Generate trial_by_trial data
VPS = cdata_meta.volume_rate;
%tic; generate_trial_by_trial_composite_behaviour_and_calcium_panels( asid, sid, cdata_raw, bdata_vel, btrial_meta, bdata_vel_time, VPS, analysis_path, new_rois ); toc;
tic; generate_trial_by_trial_composite_behaviour_and_calcium_panels( asid, sid, cdata_raw, bdata_vel, btrial_meta, bdata_vel_time, VPS, analysis_path, rois ); toc;

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
ac = get_analysis_constants;
one_trial_bdata = squeeze(b_rawdata{ ac.BOTH }(1,:,:));

right_odor_stim = squeeze(one_trial_bdata(:,6));
left_odor_stim  = squeeze(one_trial_bdata(:,7));

figure;
hold on;
plot(b_time, left_odor_stim);
first_stim = find((left_odor_stim > 2.0), 1, 'first');
last_stim = find((left_odor_stim > 2.0), 1, 'last');
plot(b_time(first_stim), 5.0, 'xr');
plot(b_time(last_stim), 5.0, 'xr');

%%
ac = get_analysis_constants;
one_trial_bdata = squeeze(b_rawdata{ ac.BOTH }(1,:,:));
frame_clock     = squeeze(one_trial_bdata(:,5));

figure;
hold on;
plot(b_time, frame_clock);

FSTATE_HIGH = 1;
FSTATE_LOW = 2; 
cur_state = FSTATE_LOW;
FRAME_STATE_CHANGE_THRESHOLD = 2.0;

frame_begins_t = [];
frame_ends_t = [];

for t=1:length(frame_clock)-1

    cur_frame_signal = frame_clock(t);
    next_frame_signal = frame_clock(t+1);
    
    if( cur_state == FSTATE_LOW )
        % disp(['cfs: ' num2str(cur_frame_signal) ' nfs: ' num2str(next_frame_signal)]);
        if((cur_frame_signal - next_frame_signal) < -1.0*FRAME_STATE_CHANGE_THRESHOLD )
            frame_begins_t(end+1) = b_time(t+1);
            cur_state = FSTATE_HIGH;
        end
    elseif( cur_state == FSTATE_HIGH )        
        if((cur_frame_signal - next_frame_signal) > FRAME_STATE_CHANGE_THRESHOLD )
            frame_ends_t(end+1) = b_time(t);
            cur_state = FSTATE_LOW;
        end
    else
        disp(['ERROR: State: ' num2str(cur_state) ' is not recognized.']);
    end
end

plot(frame_begins_t, 5.0, 'xr');
plot(frame_ends_t, 5.0, 'xb');

FPLANES = 16;

VOLUMES = length(frame_begins_t) / FPLANES;

frame_begins_in_vol = reshape( frame_begins_t, [ FPLANES, VOLUMES ] );

frame_start_offsets_per_plane_per_vol = zeros(FPLANES, VOLUMES);

for v = 1:size(frame_begins_in_vol,2)
    
    first_plane_time = frame_begins_in_vol(1,v);
    
    for p = 1:FPLANES
        frame_start_offsets_per_plane_per_vol(p,v) = frame_begins_in_vol(p,v) - first_plane_time;
    end
end
    
frame_start_offsets_per_plane = zeros(1, FPLANES);
for p = 1:FPLANES
    frame_start_offsets_per_plane(p) = squeeze(mean(frame_start_offsets_per_plane_per_vol(p,:),2));
end

if 1
f = figure;
hold on;
clear avg_offset err_in_offset;
for p = 1:size(frame_begins_in_vol,1)
    avg_offset(p) = squeeze(mean(frame_start_offsets_per_plane_per_vol(p,:),2));
    err_in_offset(p) = squeeze(std(frame_start_offsets_per_plane_per_vol(p,:)));       
end

errorbar( [1:size(frame_begins_in_vol,1)], avg_offset, err_in_offset );
xlabel('Planes');
ylabel('Frame start offset (s)');

saveas(f, [analysis_path '/frame_start_offsets_per_plane.fig']);
saveas(f, [analysis_path '/frame_start_offsets_per_plane.png']);
end

cur_plane = 4;
start_times_in_plane = squeeze(frame_begins_in_vol(cur_plane,:));

figure;
hold on;

plot(start_times_in_plane, 4, 'x');
plot([0:size(frame_begins_in_vol,2)-1]./VPS+frame_start_offsets_per_plane(cur_plane), 4, 'o');
xlim([0 6.5]);

num_points = length(find(start_times_in_plane <= 3.0));

%%%%%% ATTIC %%%%%%%

%% Temporary test space 
a_const = get_analysis_constants;
cur_trial_type = a_const.RIGHT;
cur_trial_type_str = a_const.task_str{cur_trial_type};diff_avg_trace_filepath = [ analysis_path '/' condition_str '_avg_diff_traces_asid_' num2str( asid ) '_sid_' num2str(sid) ];
display_two_behavioral_condition_diff_traces( condition_trials_str, btraces_per_condition, ctraces_in_roi_per_condition, bdata_vel_time, VPS, diff_avg_trace_filepath );

trial_ord = 4;

cur_cdata     = squeeze(cdata_raw{ cur_trial_type }(trial_ord,:,:,:,:,:));
display_two_behavioral_condition_diff_traces( condition_trials_str, btraces_per_condition, ctraces_in_roi_per_condition, bdata_vel_time, VPS, avg_trace_filepath );
% 
rois_v2 = get_rois_from_volume_v2(cur_cdata);

cur_trial_id = squeeze(btrial_meta{ cur_trial_type }(trial_ord, 2));

cur_tbt_filename_prefix = [ analysis_path '/roi_avg_volume_asid_' num2str(asid) '_sid_' num2str(sid) '_' cur_trial_type_str '_tid_' num2str(cur_trial_id)];
generate_volume_avg_with_rois( cur_cdata, rois, cur_tbt_filename_prefix );

%% Play with frame clock
ac = get_analysis_constants;
one_trial_bdata = squeeze(b_rawdata{ ac.BOTH }(1,:,:));
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

%% Save imaging data for faster load
datamat_path = [datapath slash 'data_mat'];

if(~exist(datamat_path, 'dir'))
    mkdir( datamat_path );
end

cdata_mat_path = [datamat_path '/cdata_all.mat'];
tic; savefast(cdata_mat_path, 'cdata_raw', 'cdata_meta', 'ctrial_meta'); toc
