function [ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list( bdata_vel_time, bdata_vel )
% Version 1. Simple threshold during stim

ac = get_analysis_constants;
settings = sensor_settings;

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

EXPECTED = 1;
IGNORED = 2;
condition_trials_str = {'expected_turning', 'ignoring_stim' };
condition_str = 'expected_vs_ignoring_turning';

LEFT_TURN  = 1;
RIGHT_TURN = 2;
NO_TURN    = 3;

TURN_THRESHOLD = 0.001; % ???
FWD_VELOCITY_THRESHOLD = 0.001;

condition_trials = cell(3,2);

for trial_type = 1:size( bdata_vel, 2 )
    for trial_ord = 1:size( bdata_vel{trial_type}, 1 )
        
        cur_yaw_tc = bdata_vel{ trial_type }( trial_ord, ac.VEL_YAW, : );
        cur_fwd_tc = bdata_vel{ trial_type }( trial_ord, ac.VEL_FWD, : );
        
        fwd_vel = cur_yaw_tc( find( bdata_vel_time < (prestim+stim)) );
        avg_fwd_vel = mean( fwd_vel );
                
        if( avg_fwd_vel < FWD_VELOCITY_THRESHOLD )
            continue;
        end
        
        yaw_during_stim = cur_yaw_tc( find( (bdata_vel_time > prestim) & (bdata_vel_time <= (prestim+stim))) );
        avg_yaw_during_stim = mean(yaw_during_stim);
                   
        turn_status = -1;
        if( avg_yaw_during_stim > TURN_THRESHOLD )
            turn_status = RIGHT_TURN;
        elseif( avg_yaw_during_stim < -1.0*TURN_THRESHOLD )  
            turn_status = LEFT_TURN;        
        else
            turn_status = NO_TURN;
        end
        
        if( trial_type == ac.BOTH )
            if( turn_status == NO_TURN )
                condition_trials{ trial_type, EXPECTED }( end+1 ) = trial_ord;
            else
                condition_trials{ trial_type, IGNORED }( end+1 ) = trial_ord;
            end
        elseif( trial_type == ac.RIGHT )            
            if( turn_status == RIGHT_TURN )
                condition_trials{ trial_type, EXPECTED }( end+1 ) = trial_ord;
            else
                condition_trials{ trial_type, IGNORED }( end+1 ) = trial_ord;
            end        
        elseif( trial_type == ac.LEFT )
            if( turn_status == LEFT_TURN )
                condition_trials{ trial_type, EXPECTED }( end+1 ) = trial_ord;
            else
                condition_trials{ trial_type, IGNORED }( end+1 ) = trial_ord;
            end        
        else
            disp(['ERROR: trial_type: ' num2str(trial_type)]);
        end
    end
end

end

