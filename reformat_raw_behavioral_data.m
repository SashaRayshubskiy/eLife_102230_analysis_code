function [ bdata_vel, bdata_meta ] = reformat_raw_behavioral_data( bdata_raw )
% * bdata_vel = {trial_type_id} [analysis_trial_id, vel_fwd, vel_side, vel_yaw]
%
%  Note: this data is downscaled to be at 100 Hz (the sampling rate of the
%  ADNS 9800 optical sensor.
%  Also note that the analysis trial id is a combination of sid and tid.
%  It's a unique identifier of the trial for this particular session.
%
% * bdata_meta = { trial_type_id } [ analysis_trial_id, frame_clock, stim_left, stim_right ]
% 

trial_type_cnt = size(bdata_raw, 1);

for i=1:trial_type_cnt
    trial_per_type_cnt = size(bdata_raw{i},1);
    bdata_vel{i} = zeros(trial_per_type_cnt, 3, ?? size of the downsampled behavioral array, from the behavioral settings file ??, 'double' );
    bdata_meta{i} = zeros(trial_per_type_cnt, 3, ?? size of behavioral daq sampling, in the behavioral settings file ??, 'double' );
end

for i=1:trial_type_cnt

    trial_per_type_cnt = size(bdata_raw{i},1);
    for j = 1:trial_per_type_cnt
    
    bdata_vel{i}
    
    end    
end
end

