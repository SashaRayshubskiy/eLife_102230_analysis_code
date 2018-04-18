function [ data ] = load_causality_experiment_behavioral_data( cur_sid, datapath, pre_generated_filepath )

USE_CALIBRATION = 1;
USE_PRELOADED_DATA = 1;

if(( USE_PRELOADED_DATA == 1) && (exist( pre_generated_filepath, 'file') ) )
    data = load(pre_generated_filepath);
else
    
    disp(['About to load: ' datapath ]);
    
    [ bdata_raw, bdata_time, trial_metadata ] = load_behavioral_data( cur_sid, datapath, 1 );
    
    TIME_OFFSET = bdata_time(1);
    
    % Transform raw data into velocity
    t_volts = [];
    volts_all = [];
    
    t_vel_all = [];
    yaw_vel_all = [];
    forward_vel_all = [];
            
    for trial = 1:size(bdata_raw{1},1)
        
        tid = trial_metadata{1}(trial,2);
        
        t = bdata_time;
        D = squeeze(bdata_raw{1}(trial,:,:));
        
        %[currentA, voltageA, currentB, voltageB] = get_dual_scaled_voltage_and_current( D );
        
        %volts_all(trial, :) = voltageB;
        %t_volts(trial,:) = t-TIME_OFFSET;
        
        [ t_vel, vel_forward, vel_side, vel_yaw ] = get_velocity_ephysrig(t, D, datapath, USE_CALIBRATION );
        
        t_vel_all( trial, : ) = t_vel-TIME_OFFSET;
        
        yaw_vel_all( trial, : ) = vel_yaw;
        forward_vel_all( trial, : ) = vel_forward;
    end
    
    %save(pre_generated_filepath, 't_volts', 'volts_all', 't_vel_all', 'yaw_vel_all', 'forward_vel_all');
    save(pre_generated_filepath, 't_vel_all', 'yaw_vel_all', 'forward_vel_all');
    data = load(pre_generated_filepath);
end
end

