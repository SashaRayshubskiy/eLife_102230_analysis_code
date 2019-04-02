function [ EB_data ] = load_EB_data( datapath, cur_sid )

trial_type_cnt = 1;

% Load behavioral data
bdata_path = [datapath '/ball/' ];
tic; [ b_rawdata, b_time, btrial_meta ] = load_behavioral_data( cur_sid, bdata_path, trial_type_cnt ); toc

% Get behavioral data that is usable for analysis
%[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data( b_time, b_rawdata );
[ EB_data.bdata_vel_time, EB_data.bdata_vel ] = reformat_raw_behavioral_data_berg1( b_time, b_rawdata );

% Get ephys data
[ EB_data.ephys_time, EB_data.ephys_data ] = reformat_ephys_data_berg1( b_time, b_rawdata );

% Get pico monitor data
[ EB_data.pico_stim_data ] = reformat_pico_stim_data_berg1( b_rawdata );

% Load imaging data
cdata_path = [datapath  '/2p/' ];
dx = 1;
dy = 1;
dt = 1;

tic; [ EB_data.cdata_raw, cdata_meta, ctrial_meta ] = load_imaging_data_hack( cur_sid, cdata_path, trial_type_cnt, dx, dy, dt ); toc

EB_data.VPS = cdata_meta.volume_rate;

end

