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

datapath = '/data/drive1/sasha/180207_6f_R60D05_LexAOp_CsCr_83b_lexA_03/';
%datapath = '/data/drive0/tots/160824_op6s_VT040354_01/';

analysis_path = [datapath slash 'analysis'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

sid = [ 0 ];

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finished data load.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Display bump dynamics for each trial    end_t   = start_t + nframes;

display_bump_dynamics_per_trial( cdata_raw, bdata_vel_time, bdata_vel, VPS, analysis_path );


%%
display_avg_tc(cdata_raw, VPS, analysis_path);

%% Show average tc for trials that have a counter-turn only
display_avg_tc_counter_turn_only(cdata_raw, bdata_vel_time, bdata_vel, VPS, analysis_path);

%% Show average over planes and pixels for each trial

f = figure;
for trial_type = 1:2
    for tr = 1:size(cdata_raw{trial_type},1)
        
        cur_trial = squeeze(cdata_raw{trial_type}(tr,:,:,:,:));
        
        cur_trial_all_plane_tc = squeeze(mean(squeeze(mean(squeeze(mean(cur_trial(:,:,4:end,:),3))))));
        
        cur_trial_all_plane_tc_corr = cur_trial_all_plane_tc - repmat(mean(cur_trial_all_plane_tc(1:10)), [1 length(cur_trial_all_plane_tc)]);
        
        hold on;
        plot(cur_trial_all_plane_tc_corr);
    end
end

saveas(f, [analysis_path '/average_img_per_trial.fig' ]);
saveas(f, [analysis_path '/average_img_per_trial.png' ]);

%% Display behavioral data
display_avg_velocity(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path);
%display_avg_velocity_exclude_zero_vel(sid, b_rawdata, bdata_vel, bdata_vel_time, analysis_path ); 

%% Display single trial trajectories
traj = get_single_trial_trajectories(sid, bdata_vel_time, bdata_vel);

% display_single_trial_trajectories( sid, bdata_vel_time, traj, analysis_path );
display_single_trial_trajectories_v2( sid, bdata_vel_time, traj, analysis_path );

%% Display behavioral data
%datapath_tmp = '/data/drive0/sasha//';
datapath_tmp = '/data/drive1/sasha/180129_op6s_ss730_LexAOp_CsCr_83b-lexA_03/';

analysis_path_tmp = [datapath_tmp slash 'analysis'];

if(~exist(analysis_path_tmp, 'dir'))
    mkdir(analysis_path_tmp);
end

sid_tmp = [ 0 ];
trial_type_cnt_tmp = 2;

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

%% Merge sids
datapath = '/data/drive1/sasha/180205_6f_R60D05_LexAOp_CsCr_83b_lexA_01/';
sid_1 = 0;
sid_2 = 1;

merge_sids( datapath, sid_1, sid_2);
