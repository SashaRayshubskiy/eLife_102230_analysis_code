%% Load imaging and behavioral data
global slash;

clear all;

if isunix() == 1
    slash = '/';
else
    slash = '\';
end

% sid 0
% datapath = '/data/drive3/sasha/190504_VT7338_lexA_LexAOp_GCaMP7f_01';

% datapath = '/data/drive3/sasha/190506_VT7338_lexA_LexAOp_GCaMP7f_02';
datapath = '/data/drive3/sasha/190521_VT007338_Gal4_UAS_GCaMP7f_01';

% datapath = '/data/drive3/sasha/190507_VT7338_lexA_LexAOp_GCaMP7f_03';

analysis_path = [datapath slash 'analysis'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

sid = [ 3 ];

aconstants = get_analysis_constants;
trial_type_cnt = 1; 

% Load behavioral data
bdata_path = [datapath  slash 'ball' slash ];
tic; [ b_rawdata, b_time, btrial_meta ] = load_behavioral_data(sid, bdata_path, trial_type_cnt ); toc

% Get behavioral data that is usable for analysis
%[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data( b_time, b_rawdata );
[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data_berg1( b_time, b_rawdata );

% Get ephys data
% [ephys_time, ephys_data] = reformat_ephys_data_berg1( b_time, b_rawdata );

% Get pico monitor data
% [ pico_stim_data ] = reformat_pico_stim_data_berg1( b_rawdata );

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

imaging_data_1       = squeeze(cdata_raw{1}( :, :, :, 1, START_PLANE:end, : ));
% Get max for each plane
% imaging_data = squeeze(max(imaging_data_1, [], 4));

%% Examine ROIs
[ rois ] = get_PFLm_rois( imaging_data_1, analysis_path );


%% Display EPG bump, yaw, and ephys on the same plot for quick view of experiment
[ PFL_LAL_dF_F_per_trial, PFL_FB_dF_F_per_trial ] = display_PFLm_dynamics( imaging_data_1, bdata_vel_time, bdata_vel, VPS, analysis_path, sid );


%% Needs the above
% Analyze the difference between left and right PFLm tuft and yaw
analyze_left_right_PFLm_diff_and_yaw( imaging_data_1, PFL_FB_dF_F_per_trial, PFL_LAL_dF_F_per_trial, bdata_vel_time, bdata_vel, analysis_path, VPS, sid );

%%
% Analyze the difference between left and right PFLm tuft and fwd
analyze_left_right_PFLm_diff_and_fwd( imaging_data_1, PFL_FB_dF_F_per_trial, PFL_LAL_dF_F_per_trial, bdata_vel_time, bdata_vel, analysis_path, VPS, sid );


%% Trigger on events of either left or right PFL.LAL bias, examing the average turning signal.

display_PFL_left_right_triggered_yaw( PFL_LAL_dF_F_per_trial, bdata_vel, imaging_data_1, VPS, analysis_path );

%% Plot yaw vs. left-right difference in PFLm.LAL GCaMP signal

ac = get_analysis_constants;

num_trials = size( imaging_data_1, 1 );
nframes = size( imaging_data_1, 5 );

t = [ 0 : nframes-1 ] ./ VPS;

FFT_FILTER_CUTOFF = 1;

DT = 8; % Imaging frame rate is approximately 80 ms per frame

f = figure;

yaw_values = [];
left_right_PFL_diff = [];

for tr = 1:num_trials

    subplot( 1, 1, 1 );

    cur_PFL_LAL = squeeze(PFL_LAL_dF_F_per_trial(tr, :,:));        
    
    left_PFLm_LAL_axon_tract  = squeeze( cur_PFL_LAL(1,:));
    right_PFLm_LAL_axon_tract  = squeeze( cur_PFL_LAL(3,:));
    
    left_PFLm_LAL_axon_tract_filt   = fft_filter( left_PFLm_LAL_axon_tract', FFT_FILTER_CUTOFF, VPS );    
    right_PFLm_LAL_axon_tract_filt  = fft_filter( right_PFLm_LAL_axon_tract', FFT_FILTER_CUTOFF, VPS );    

    left_right_PB_LAL_delta_tract = left_PFLm_LAL_axon_tract_filt - right_PFLm_LAL_axon_tract_filt;
    left_right_PB_LAL_delta_tract = left_right_PB_LAL_delta_tract - mean(left_right_PB_LAL_delta_tract);
    
    cur_yaw = squeeze(bdata_vel{ 1 }( tr, ac.VEL_YAW, : ));
    cur_yaw = cur_yaw - mean( cur_yaw );
    
    cur_yaw_len = floor(length(cur_yaw)/DT) * DT;
    
    cur_yaw_d = squeeze(mean(reshape( cur_yaw(1:cur_yaw_len), [DT, cur_yaw_len/DT])));
    
    min_len = min(length(cur_yaw_d), length(left_right_PB_LAL_delta_tract));
    
    yaw_values = horzcat( yaw_values, cur_yaw_d(1:min_len));
    left_right_PFL_diff = horzcat( left_right_PFL_diff, left_right_PB_LAL_delta_tract(1:min_len)');

    assert(length(yaw_values) == length( left_right_PFL_diff ) );
end

hold on;
plot( left_right_PFL_diff, yaw_values, 'o', 'MarkerSize', 3, 'color', rgb('DarkGray') );

ylabel('Yaw');
xlabel('Left-Right PFL LAL dF/F');

saveas( f, [analysis_path '/left_right_PFL_LAL_diff_vs_yaw.fig'] );
saveas( f, [analysis_path '/left_right_PFL_LAL_diff_vs_yaw.png'] );















