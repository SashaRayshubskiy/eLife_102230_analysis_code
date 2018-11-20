function [ bdata_time_out, bdata_vel ] = reformat_raw_behavioral_data_berg1( bdata_time_in, bdata_raw )
% * bdata_vel = {trial_type_id} [analysis_trial_id, vel_fwd, vel_side, vel_yaw]
%
%  Note: this data is downscaled to be at 100 Hz (the sampling rate of the
%  ADNS 9800 optical sensor.
%  Also note that the analysis trial id is a combination of sid and tid.
%  It's a unique identifier of the trial for this particular session.
% 

trial_type_cnt = size(bdata_raw, 2);

settings = sensor_settings();

for i=1:trial_type_cnt
    trial_per_type_cnt = size(bdata_raw{i},1);
    
    sensor_sampling = settings.sensorPollFreq;
    bdaq_sampling = size(bdata_raw{i},2);
    
    [ bdata_time_out, dummy1, dummy2, dummy3 ] = get_velocity_berg1(bdata_time_in, squeeze(bdata_raw{i}(1,:,:)));
    
    bdata_vel{ i } = zeros(trial_per_type_cnt, 3, length(bdata_time_out), 'double' );
    
    for j = 1:trial_per_type_cnt    
        [ bdata_time_out, vel_fwd, vel_side, vel_yaw ] = get_velocity_berg1(bdata_time_in, squeeze(bdata_raw{i}(j, :,:)));
        bdata_vel{i}(j, 1, :) = vel_fwd;    
        bdata_vel{i}(j, 2, :) = vel_side;    
        bdata_vel{i}(j, 3, :) = vel_yaw;
    end        
end

end

