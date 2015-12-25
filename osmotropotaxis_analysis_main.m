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
[ b_rawdata, b_time, btrial_metadata ] = load_behavioral_data(sid, bdata_path);

% Load imaging data
cdata_path = [datapath  slash '2p' slash ];
[ cdata_raw, ctrial_metadata ] = load_imaging_data(sid, cdata_path);

% Check that the behavioral and imaging trials match up
check_bdata_and_cdata_trial_integrity( btrial_metadata, ctrial_metadata );

% 
[bdata_vel, bdata_meta] = reformat_raw_behavioral_data( b_rawdata );

%% Display behavioral data
