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

datapath = '/data/drive2/sasha/180625_R12D09_P2X2_op6s_60D05_01/';

analysis_path = [datapath slash 'analysis'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

sid = [ 1 ];

aconstants = get_analysis_constants;
trial_type_cnt = 1; 

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

%%
display_PB_roi_dynamics( cdata_raw, bdata_vel_time, bdata_vel, VPS, analysis_path );

