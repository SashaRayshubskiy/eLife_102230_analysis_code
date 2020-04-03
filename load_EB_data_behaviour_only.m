function [ EB_data ] = load_EB_data_behaviour_only( datapath, cur_sid )

trial_type_cnt = 1;

% Load behavioral data
bdata_path = [datapath '/ball/' ];
tic; [ b_rawdata, b_time, btrial_meta ] = load_behavioral_data( cur_sid, bdata_path, trial_type_cnt ); toc

% Get behavioral data that is usable for analysis
%[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data( b_time, b_rawdata );
[ EB_data.bdata_vel_time, EB_data.bdata_vel ] = reformat_raw_behavioral_data_berg1( b_time, b_rawdata );

end

