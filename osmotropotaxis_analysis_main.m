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

sid = 4;

% Load behavioral data
bdata_path = [datapath  slash 'ball' slash ];
tic; [ b_rawdata, b_time, btrial_meta ] = load_behavioral_data(sid, bdata_path); toc

% Load imaging data
cdata_path = [datapath  slash '2p' slash ];
tic; [ cdata_raw, cdata_meta, ctrial_meta ] = load_imaging_data(sid, cdata_path); toc

% Check that the behavioral and imaging trials match up
check_bdata_and_cdata_trial_integrity( btrial_meta, ctrial_meta );

% Get behavioral data that is usable for analysis
[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data( b_time, b_rawdata );

%% Display behavioral data
