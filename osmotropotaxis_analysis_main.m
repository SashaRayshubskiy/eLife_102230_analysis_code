%% Step 1: Load imaging and behavioral data
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

datapath = '/data/drive1/sasha/180210_gcamp6f_ss730_LexAOp_CsCr_83b-lexA_01/';
%datapath = '/data/drive0/tots/160824_op6s_VT040354_01/';

analysis_path = [datapath slash 'analysis'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

sid = [ 1 ];

aconstants = get_analysis_constants;
trial_type_cnt = 2; 

% Load behavioral data
bdata_path = [datapath  slash 'ball' slash ];
tic; [ b_rawdata, b_time, btrial_meta ] = load_behavioral_data(sid, bdata_path, trial_type_cnt ); toc

% Load imaging data
cdata_path = [datapath  slash '2p' slash ];
dx = 1;
dy = 1;
dt = 1;
down_params(1) = dx;
down_params(2) = dy;

if( dt == 1 )
    tic; [ cdata_raw, cdata_meta, ctrial_meta ] = load_imaging_data(sid, cdata_path, trial_type_cnt, dx, dy, dt ); toc
else
    tic; [ cdata_raw, cdata_meta, ctrial_meta ] = load_imaging_data_single_plane(sid, cdata_path, trial_type_cnt, dx, dy, dt ); toc
end

% Check that the behavioral and imaging trials match up
%check_bdata_and_cdata_trial_integrity( btrial_meta, ctrial_meta );

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

%% Step 2: Get behavioral data analysis
display_avg_velocity_exclude_zero_vel_RL_only(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);
display_avg_velocity_left_right_only(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);
display_per_trial_velocity(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);

%% Step 3: Generate turning metadata and select expected turning vs. ignoring trials
turn_metadata = generate_turning_metadata( sid, bdata_vel_time, bdata_vel, analysis_path );

[ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list_v2( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path );

asid= 0;
avg_cond_btrace_trace_filepath = [ analysis_path '/' condition_str '_asid_' num2str( asid ) '_sid_' num2str(sid) ];
with_single_trials = 1;
display_two_condition_trials( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, with_single_trials );

with_single_trials = 0;
display_two_condition_trials_avg( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, with_single_trials );

%% Step 4: Generate df/f time courses for both conditions (expected turning vs ingnoring stim)

if (dt == 1 )
    tic; [ btraces_per_condition, avg_df_f_per_condition_per_plane ] = collect_two_behavioral_condition_and_df_f_per_cond( condition_trials, cdata_raw, bdata_vel, VPS, trial_exclusion_list, btrial_meta ); toc;
else
    FPS = VPS/dt;
    tic; [ btraces_per_condition, avg_df_f_per_condition_per_plane ] = collect_two_behavioral_condition_and_df_f_per_cond_single_plane( condition_trials, cdata_raw, bdata_vel, FPS, trial_exclusion_list, btrial_meta ); toc;
end

%% Step 5: Examine mean image from every plane

a_const = get_analysis_constants;

cur_trial_type = a_const.RIGHT;
cur_trial_type_str = a_const.task_str{cur_trial_type};
trial_id = 4;

cur_cdata     = squeeze(cdata_raw{ cur_trial_type }(trial_id,:,:,:,:,:));
figsave_prefix = [analysis_path '/avg_volume_' cur_trial_type_str '_tid_' num2str(trial_id) ];
display_avg_volume( down_params, cur_cdata, figsave_prefix );

%% Step 6: Get ROIs and show avg tc for each

ac = get_analysis_constants();

roi_session = 100;
PLANE_OF_INTEREST = 2;
ref_img = squeeze(mean(mean(squeeze(cdata_raw{ 1 }(40:45,:,:,PLANE_OF_INTEREST,:)),4),1));

diff_image_path = [ analysis_path '/' condition_str '_sid_' num2str(sid) '_plane_' num2str(PLANE_OF_INTEREST) '_both_sides_roi_session_' num2str(roi_session) ];
clicky_two_condition_bdata_LR(ref_img, PLANE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% Display behavioral data
display_avg_velocity(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);
%display_avg_velocity_exclude_zero_vel(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path ); 

%% Display single trial trajectories
traj = get_single_trial_trajectories(sid, bdata_vel_time, bdata_vel);

% display_single_trial_trajectories( sid, bdata_vel_time, traj, analysis_path );
display_single_trial_trajectories_v2( sid, bdata_vel_time, traj, analysis_path );

%% Display behavioral data
%datapath_tmp = '/data/drive0/sasha//';
slash = '/';
datapath_tmp = '/data/drive1/sasha/180220_RIW5_CsCr_01/';

analysis_path_tmp = [datapath_tmp slash 'analysis'];

if(~exist(analysis_path_tmp, 'dir'))
    mkdir(analysis_path_tmp);
end

sid_tmp = [ 1 ];
trial_type_cnt_tmp = 1; 

bdata_path_tmp = [datapath_tmp  slash 'ball' slash ];
tic; [ b_rawdata_tmp, b_time_tmp, btrial_meta_tmp ] = load_behavioral_data(sid_tmp, bdata_path_tmp, trial_type_cnt_tmp ); toc

[bdata_vel_time_tmp, bdata_vel_tmp] = reformat_raw_behavioral_data( b_time_tmp, b_rawdata_tmp );

if( trial_type_cnt_tmp == 3 )
    display_avg_velocity(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);
    display_avg_velocity_exclude_zero_vel(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp); 
elseif( trial_type_cnt_tmp == 2 )
    %display_avg_velocity_exclude_zero_vel_RL_only(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);    
    %display_avg_velocity_exclude_zero_vel_both_only(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);    
    display_avg_velocity_left_right_only(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);
    %display_avg_velocity_left_right_both(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);
    %display_per_trial_velocity(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);
elseif( trial_type_cnt_tmp == 1 )
    display_avg_velocity_exclude_zero_vel_one_side_only(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);    
    display_avg_velocity_one_side_only(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);
    display_2p_stim_per_trial_velocity(sid_tmp, b_rawdata_tmp, bdata_vel_tmp, bdata_vel_time_tmp, analysis_path_tmp);
end

%%
CLUSTER_DATA = 1;
if (CLUSTER_DATA == 1)
    % Generate turning clusters using standard unsupervised clustering tools.
    MAX_CLUSTERS = 3;

    plot_idx = 2;
    
    clust = generate_turning_clusters_v2( sid_tmp, bdata_vel_time_tmp, bdata_vel_tmp, analysis_path_tmp, MAX_CLUSTERS, plot_idx);
end

%% Generate turning clusters using standard unsupervised clustering tools.
MAX_CLUSTERS = 5;
clust = generate_turning_clusters( sid, bdata_vel_time, bdata_vel, analysis_path, MAX_CLUSTERS );

%% Test the idea that activation of 'LAL' turning region between the 2 sides corresponds to turn magnitude.

generate_turning_magnitude_vs_bilateral_calcium_delta_response_plot( bdata_vel_time, bdata_vel, cdata_raw, frame_start_offsets_per_plane, VPS, analysis_path );

%% Generate turning metadataum

turn_metadata = generate_turning_metadata( sid, bdata_vel_time, bdata_vel, analysis_path );

%% Generate expected vs. ignored
% [ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list( bdata_vel_time, bdata_vel );
[ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list_v2( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path );

%% Generate expected vs. ignored
% [ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list( bdata_vel_time, bdata_vel );
[ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list_v3( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path );

%% Generate stationary vs. walking 
[ condition_trials, condition_trials_str, condition_str ] = generate_stationary_vs_walking( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path );

%% Generate large vs. small counter turn

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
asid= 0;
avg_cond_btrace_trace_filepath = [ analysis_path '/' condition_str '_asid_' num2str( asid ) '_sid_' num2str(sid) ];
with_single_trials = 1;
display_two_condition_trials( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, with_single_trials );

with_single_trials = 0;
display_two_condition_trials_avg( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, with_single_trials );

%% Display behavioral 2 condition trials.

%condition_trials{2}(2) = 14;
asid= 0;
avg_cond_btrace_trace_filepath = [ analysis_path '/' condition_str '_asid_' num2str( asid ) '_sid_' num2str(sid) ];
with_single_trials = 1;
display_two_condition_trials_fwd( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, with_single_trials );

with_single_trials = 0;
display_two_condition_trials_avg_fwd( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, with_single_trials );

%% Display trial running trajectories for each analysis condition

traj = get_single_trial_trajectories(sid, bdata_vel_time, bdata_vel);
    
display_two_condition_single_trial_trajectories( sid, condition_trials, condition_trials_str, bdata_vel_time, traj, analysis_path );

%% Get covariance plots from 2 symmetric ROIs
ac = get_analysis_constants();

PLANE_OF_INTEREST = 9;

roi_session = 0;
rois_path = [analysis_path '/roi_session_' num2str(roi_session) '.mat' ];
rois_1 = get_rois_from_volume( PLANE_OF_INTEREST, squeeze(cdata_raw{ 1 }(1,:,:,:,:,:)), rois_path );

% generate_df_f_vs_turning_magnitude_trial_plot( sid, rois_1, PLANE_OF_INTEREST, cdata_raw, turn_metadata, frame_start_offsets_per_plane, VPS, btrial_meta, analysis_path );
generate_df_f_delta_vs_turning_magnitude_trial_plot( sid, rois_1, PLANE_OF_INTEREST, cdata_raw, turn_metadata, frame_start_offsets_per_plane, VPS, btrial_meta, analysis_path );

%% Create a differece image for each plane
if (dt == 1 )
    tic; [ btraces_per_condition, avg_df_f_per_condition_per_plane ] = collect_two_behavioral_condition_and_df_f_per_cond( condition_trials, cdata_raw, bdata_vel, VPS, trial_exclusion_list, btrial_meta ); toc;
else
    FPS = VPS/dt;
    tic; [ btraces_per_condition, avg_df_f_per_condition_per_plane ] = collect_two_behavioral_condition_and_df_f_per_cond_single_plane( condition_trials, cdata_raw, bdata_vel, FPS, trial_exclusion_list, btrial_meta ); toc;
end


%%
roi_session = 30;

diff_image_path = [ analysis_path '/' condition_str '_diff_image_asid_' num2str( asid ) '_sid_' num2str(sid) '_roi_session_' num2str(roi_session) ];
display_two_condition_difference_image(down_params, ref_imgs, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% Clicky showing avg data for 2, with ialboth left and right trials. Pick 4 rois on each side.

ac = get_analysis_constants();

roi_session = 2;
PLANE_OF_INTEREST = 13;
ref_img = squeeze(mean(mean(squeeze(cdata_raw{ 1 }(1:5,:,:,PLANE_OF_INTEREST,:)),4),1));

diff_image_path = [ analysis_path '/' condition_str '_sid_' num2str(sid) '_plane_' num2str(PLANE_OF_INTEREST) '_roi_session_' num2str(roi_session) ];
clicky_two_condition_bdata_8roi_LR(ref_img, PLANE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% Clicky showing avg data for 2 conditionsquick

ac = get_analysis_constants();

roi_session = 11;
PLANE_OF_INTEREST = 10;
TRIAL_TYPE_OF_INTEREST = ac.RIGHT;
ref_img = squeeze(mean(mean(squeeze(cdata_raw{ 1 }(40:45,:,:,PLANE_OF_INTEREST,:)),4),1));

diff_image_path = [ analysis_path '/' condition_str '_sid_' num2str(sid) '_plane_' num2str(PLANE_OF_INTEREST) '_side_' num2str(TRIAL_TYPE_OF_INTEREST) '_roi_session_' num2str(roi_session) ];
clicky_two_condition_bdata(ref_img, PLANE_OF_INTEREST, TRIAL_TYPE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% Clicky showing both left and right trials
ac = get_analysis_constants();

roi_session = 0;
PLANE_OF_INTEREST = 12;
ref_img = squeeze(mean(mean(squeeze(cdata_raw{ 1 }(10:20,:,:,PLANE_OF_INTEREST,:)),4),1));

diff_image_path = [ analysis_path '/' condition_str '_sid_' num2str(sid) '_plane_' num2str(PLANE_OF_INTEREST) '_both_sides_roi_session_' num2str(roi_session) ];
clicky_two_condition_bdata_LR(ref_img, PLANE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% Clicky showing avg data for 2 conditionsquick, version 2
ac = get_analysis_constants();

roi_session = 0;
PLANE_OF_INTEREST = 9;
ref_img = squeeze(mean(mean(squeeze(cdata_raw{ 3 }(1:5,:,:,PLANE_OF_INTEREST,:)),4),1));

diff_image_path = [ analysis_path '/' condition_str '_tc_left_right_asid_' num2str( asid ) '_sid_' num2str(sid) '_roi_session_' num2str(roi_session) ];
clicky_two_condition_bdata_v2(ref_img, PLANE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% SINGLE PLANE version
ac = get_analysis_constants();

roi_session = 1;
TRIAL_TYPE_OF_INTEREST = ac.RIGHT;
ref_img = squeeze(mean(mean(cdata_raw{ 1 }(1:5,:,:,:),4),1));

diff_image_path = [ analysis_path '/' condition_str '_diff_image_sid_' num2str(sid) '_roi_session_' num2str(roi_session) ];
FPS  = VPS / dt;
clicky_two_condition_bdata_single_plane(ref_img, TRIAL_TYPE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, FPS, diff_image_path );

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
                        cur_trial_type_str = a_const.task_str{ cur_trial_type };
                        trial_ord = 4;

                        cur_cdata     = squeeze(cdata_raw{ cur_trial_type }(trial_ord,:,:,:,:,:));
                        cur_trial_id = squeeze(btrial_meta{ cur_trial_type }(trial_ord, 2));
display_avg_velocity_exclude_zero_vel_RL_only(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);
display_avg_velocity_left_right_only(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);
display_per_trial_velocity(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);

%% Generate turning clusters using standard unsupervised clustering tools.
MAX_CLUSTERS = 5;
clust = generate_turning_clusters( sid, bdata_vel_time, bdata_vel, analysis_path, MAX_CLUSTERS );


%% Test the idea that activation of 'LAL' turning region between the 2 sides corresponds to turn magnitude.

generate_turning_magnitude_vs_bilateral_calcium_delta_response_plot( bdata_vel_time, bdata_vel, cdata_raw, frame_start_offsets_per_plane, VPS, analysis_path );

%% Generate turning metadata

turn_metadata = generate_turning_metadata( sid, bdata_vel_time, bdata_vel, analysis_path );

%% Generate expected vs. ignored
% [ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list( bdata_vel_time, bdata_vel );
[ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list_v2( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path );

%% Generate large vs. small counter turn

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
asid = 0;
avg_cond_btrace_trace_filepath = [ analysis_path '/' condition_str '_11_asid_' num2str( asid ) '_sid_' num2str(sid) ];
with_single_trials = 1;
display_two_condition_trials( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, with_single_trials );

with_single_trials = 0;
display_two_condition_trials_avg( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, with_single_trials );

%% Display trial running trajectories for each analysis condition

traj = get_single_trial_trajectories(sid, bdata_vel_time, bdata_vel);

display_two_condition_single_trial_trajectories( sid, condition_trials, condition_trials_str, bdata_vel_time, traj, analysis_path );

%% Get covariance plots from 2 symmetric ROIs
ac = get_analysis_constants();

PLANE_OF_INTEREST = 9;

roi_session = 0;
rois_path = [analysis_path '/roi_session_' num2str(roi_session) '.mat' ];
rois_1 = get_rois_from_volume( PLANE_OF_INTEREST, squeeze(cdata_raw{ 1 }(1,:,:,:,:,:)), rois_path );

% generate_df_f_vs_turning_magnitude_trial_plot( sid, rois_1, PLANE_OF_INTEREST, cdata_raw, turn_metadata, frame_start_offsets_per_plane, VPS, btrial_meta, analysis_path );
generate_df_f_delta_vs_turning_magnitude_trial_plot( sid, rois_1, PLANE_OF_INTEREST, cdata_raw, turn_metadata, frame_start_offsets_per_plane, VPS, btrial_meta, analysis_path );

%% Create a differece image for each plane
if (dt == 1 )
    tic; [ btraces_per_condition, avg_df_f_per_condition_per_plane ] = collect_two_behavioral_condition_and_df_f_per_cond( condition_trials, cdata_raw, bdata_vel, VPS, trial_exclusion_list, btrial_meta ); toc;
else
    FPS = VPS/dt;
    tic; [ btraces_per_condition, avg_df_f_per_condition_per_plane ] = collect_two_behavioral_condition_and_df_f_per_cond_single_plane( condition_trials, cdata_raw, bdata_vel, FPS, trial_exclusion_list, btrial_meta ); toc;
end

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

%% Clicky showing avg data for 2, with ialboth left and right trials. Pick 4 rois on each side.

ac = get_analysis_constants();

roi_session = 2;
PLANE_OF_INTEREST = 9;
ref_img = squeeze(mean(mean(squeeze(cdata_raw{ 1 }(1:5,:,:,PLANE_OF_INTEREST,:)),4),1));

diff_image_path = [ analysis_path '/' condition_str '_sid_' num2str(sid) '_plane_' num2str(PLANE_OF_INTEREST) '_roi_session_' num2str(roi_session) ];
clicky_two_condition_bdata_8roi_LR(ref_img, PLANE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% Clicky showing avg data for 2 conditionsquick

ac = get_analysis_constants();

roi_session = 0;
PLANE_OF_INTEREST = 9;
TRIAL_TYPE_OF_INTEREST = ac.RIGHT;
ref_img = squeeze(mean(mean(squeeze(cdata_raw{ 1 }(100:145,:,:,PLANE_OF_INTEREST,:)),4),1));

diff_image_path = [ analysis_path '/' condition_str '_sid_' num2str(sid) '_plane_' num2str(PLANE_OF_INTEREST) '_side_' num2str(TRIAL_TYPE_OF_INTEREST) '_roi_session_' num2str(roi_session) ];
clicky_two_condition_bdata(ref_img, PLANE_OF_INTEREST, TRIAL_TYPE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% Plot calcium in 2 ROIs along with behavioral data
ac = get_analysis_constants();

roi_session = 0;
PLANE_OF_INTEREST = 12;
TRIAL_TYPE_OF_INTEREST = ac.LEFT;
ref_img = squeeze(mean(mean(squeeze(cdata_raw{ 1 }( 100:145,:,:,PLANE_OF_INTEREST,:)),4),1) );

diff_image_path = [ analysis_path '/' condition_str '_sid_' num2str(sid) '_plane_' num2str(PLANE_OF_INTEREST) '_side_' num2str(TRIAL_TYPE_OF_INTEREST) '_roi_session_' num2str(roi_session) ];
plot_gcamp_in_roi_with_variable_turning_behav( condition_trials, bdata_vel, ref_img, PLANE_OF_INTEREST, TRIAL_TYPE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );


%% Clicky showing avg data for 2 conditionsquick, version 2
ac = get_analysis_constants();

roi_session = 0;
PLANE_OF_INTEREST = 8;
ref_img = squeeze(mean(mean(squeeze(cdata_raw{ 3 }(1:5,:,:,PLANE_OF_INTEREST,:)),4),1));

diff_image_path = [ analysis_path '/' condition_str '_tc_left_right_asid_' num2str( asid ) '_sid_' num2str(sid) '_roi_session_' num2str(roi_session) ];
clicky_two_condition_bdata_v2(ref_img, PLANE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_image_path );

%% SINGLE PLANE version
ac = get_analysis_constants();

roi_session = 1;
TRIAL_TYPE_OF_INTEREST = ac.RIGHT;
ref_img = squeeze(mean(mean(cdata_raw{ 1 }(1:5,:,:,:),4),1));

diff_image_path = [ analysis_path '/' condition_str '_diff_image_sid_' num2str(sid) '_roi_session_' num2str(roi_session) ];
FPS  = VPS / dt;
clicky_two_condition_bdata_single_plane(ref_img, TRIAL_TYPE_OF_INTEREST, condition_trials_str, btraces_per_condition, avg_df_f_per_condition_per_plane, bdata_vel_time, FPS, diff_image_path );

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
asid = 0;
rois = get_rois_from_volume_v2( asid, squeeze(cdata_raw{ 1 }(1,:,:,:,:,:)), analysis_path );

%%
tic; [ btraces_per_condition, ctraces_in_roi_per_condition ] = collect_two_behavioral_condition_traces( condition_trials, cdata_raw, bdata_vel, VPS, rois, trial_exclusion_list, btrial_meta ); toc;

avg_trace_filepath = [ analysis_path '/' condition_str '_avg_traces_asid_' num2str( asid ) '_sid_' num2str(sid) ];
display_two_behavioral_condition_traces( condition_trials_str, btraces_per_condition, ctraces_in_roi_per_condition, bdata_vel_time, frame_start_offsets_per_plane, VPS, avg_trace_filepath );
    
diff_avg_trace_filepath = [ analysis_path '/' condition_str '_avg_diff_traces_asid_' num2str( asid ) '_sid_' num2str(sid) ];
display_two_behavioral_condition_diff_traces( condition_trials_str, btraces_per_condition, ctraces_in_roi_per_condition, bdata_vel_time, frame_start_offsets_per_plane, VPS, diff_avg_trace_filepath );

%% Generate time courses for each plane, using the ROIs of a session

cur_bdata_vel = squeeze(bdata_vel{ cur_trial_type }( trial_ord,:,: ));

cur_tbt_filename_prefix = [ analysis_path '/time_courses_in_volume_asid_' num2str(asid) '_sid_' num2str(sid) '_' cur_trial_type_str '_tid_' num2str(cur_trial_id)];
generate_volume_time_courses_with_rois( cur_cdata, cur_bdata_vel, bdata_vel_time, new_rois, frame_start_offsets_per_plane, VPS, cur_tbt_filename_prefix );

%% Add an ROI to a plane
plane_to_append = 5;
                        new_rois = add_rois_from_volume(cur_cdata, rois_v2, plane_to_append);

%% Generate ROIs
asid = 0; % roi analysis session id

rois = get_rois_from_volume_v2( asid, squeeze(cdata_raw{ 1 }(1,:,:,:,:,:)), analysis_path );

a_const = get_analysis_constants;
cur_trial_type = a_const.LEFT;
cur_trial_type_str = a_const.task_str{ cur_trial_type };
trial_ord = 76;

cur_cdata     = squeeze(cdata_raw{ cur_trial_type }(trial_ord,:,:,:,:,:));
cur_trial_id = squeeze(btrial_meta{ cur_trial_type }(trial_ord, 2));

cur_tbt_filename_prefix = [ analysis_path '/roi_avg_volume_asid_' num2str(asid) '_sid_' num2str(sid) '_' cur_trial_type_str '_tid_' num2str(cur_trial_id)];
generate_volume_avg_with_rois( cur_cdata, rois, cur_tbt_filename_prefix );

%% Generate trial_by_trial data for single plane imaging runs

asid = 0; % roi analysis session id

rois = get_rois_from_single_plane( asid, squeeze(cdata_raw{ 1 }(1,:,:,:)), analysis_path );

FPS = cdata_meta.frame_rate / dt;

tic; generate_tt_trial_single_plane( asid, sid, cdata_raw, bdata_vel, btrial_meta, bdata_vel_time, FPS, analysis_path, rois ); toc;


%% Generate trial_by_trial data
VPS = cdata_meta.volume_rate;
%tic; generate_trial_by_triwww.googl.ecomwwwal_composite_behaviour_and_calcium_panels( asid, sid, cdata_raw, bdata_vel, btrial_meta, bdata_vel_time, VPS, analysis_path, new_rois ); toc;
tic; generate_trial_by_trial_composite_behaviour_and_calcium_panels( asid, sid, cdata_raw, bdata_vel, btrial_meta, bdata_vel_time, VPS, analysis_path, rois ); toc;

%% Generate all raw cdata on the same axis as behavioral data.

PLANE = 15;
%plot_all_behaviour_with_gcamp_in_roi( btrial_meta, cdata_raw, bdata_vel, bdata_vel_time, frame_start_offsets_per_plane, PLANE, VPS );

plot_all_behaviour_with_gcamp_in_roi_in_2D( btrial_meta, cdata_raw, bdata_vel, bdata_vel_time, frame_start_offsets_per_plane, PLANE, VPS );


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

%% Generate stack

%stack_datapath = '/data/drive0/sasha/160522_op6s_R13G10_03/2p/';
stack_datapath = '/data/drive2/sasha/180624_R12D09_GFP_stack_01/2p/';

stack_files = {'stack_00001.tif', 'stack_00002.tif', 'stack_00003.tif', ... 
               'stack_00004.tif', 'stack_00005.tif', 'stack_00006.tif', ...
               'stack_00007.tif', 'stack_00008.tif', 'stack_00009.tif', ...
                };

for s = 1:length(stack_files)
    cur_path = [stack_datapath stack_files{ s }];
    cur_stack = generate_stack_movie( cur_path );    
    if( s == 1 )
        stack = cur_stack;
    else
        stack = stack + cur_stack;
    end
    
    disp(['Loaded stack: ' cur_path]);
end

stack_f = squeeze(mean(stack,4)) ./ length(stack_files);


%% Visualize stack

stack_path = [ stack_datapath 'stack' ];
display_stack_one_chan( stack_f, stack_path );



%%
frame = 300;

figure;

for f = 280:310
imagesc(squeeze(stack(:,:,1,f)));
colormap gray;
caxis([0 15000]);
axis image;
title(['frame: ' num2str(f)]);
waitforbuttonpress;
end

%%
%%% Display stack
display_stack_rgb_1( stack, stack_datapath );
%display_stack_rgb( stack, stack_datapath );
%channel = 1;
%display_stack( channel, stack, stack_path );

%channel = 2;
%display_stack( channel, stack, stack_path );

%%
cur_path = '/data/drive0/sasha/160719_ablation_gfp_vt043158_04/2p/';
stack_path = [ cur_path 'post_left_LAL_ablation_00001.tif' ];

stack = generate_stack_movie( stack_path );

stack_save_path = [cur_path 'stack_'];
display_stack_one_chan(stack, stack_save_path );

%%
%cur_path = '/data/drive0/sasha/160721_ablation_gfp_vt025718_02/2p/';
cur_path = '/data/drive1/sasha/180219_ss730_2xGFP_01/2p/';

stack_pre = 'left_A2_stack_zoom_to_LAL_EPA_00002';

N = 1;
cur_stack = [];
for i=1:N

    stack_path = [ cur_path stack_pre '.tif' ];

    stack_tmp = generate_stack_movie( stack_path );
    
    if( i==1 )
        cur_stack = stack_tmp;
    else
        cur_stack = cur_stack + stack_tmp;
    end

end

stack = cur_stack ./ N;

%% Display stack 
stack_save_path = [cur_path stack_pre 'stack_'];
display_stack_one_chan(stack, stack_save_path );

%% Read in a snapshot and display it.

load_path = '/data/drive0/sasha/160713_test_2p_ablation_01/';
load_filename = [ 'sid_0_cavitation_test_left__00001.tif' ];

load_file = [load_path load_filename];

raw_data = open_tif_fast_single_plane( load_file );

f = figure;
imagesc(squeeze(mean(squeeze(raw_data),3)));
colormap gray;

saveas(f, [load_path load_filename '.fig']);
saveas(f, [load_path load_filename '.png']);
