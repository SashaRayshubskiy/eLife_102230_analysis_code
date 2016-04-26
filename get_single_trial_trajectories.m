function [ traj ] = get_single_trial_trajectories( sid, bdata_vel_time, bdata_vel )

settings = sensor_settings;
ac = get_analysis_constants;

dt = 1.0/settings.sensorPollFreq;

traj = cell(1,length(bdata_vel));

for trial_type = 1:length(bdata_vel)
    for trial = 1:size(bdata_vel{trial_type},1)

        vel_fwd = bdata_vel{trial_type}( trial, ac.VEL_FWD, : );
        vel_yaw = bdata_vel{trial_type}( trial, ac.VEL_YAW, : );
        vel_side = bdata_vel{trial_type}( trial, ac.VEL_SIDE, : );
        
        [disp_x, disp_y, theta] = calculate_fly_position_with_yaw( vel_fwd, vel_side, vel_yaw, dt, 0, 0, 0);
        %[disp_x, disp_y] = calculate_fly_position_no_yaw(vel_forward, vel_side, dt, 0, 0);
        
        traj{ trial_type }( trial, 1, : ) = disp_x;
        traj{ trial_type }( trial, 2, : ) = disp_y;        
    end
end

end

