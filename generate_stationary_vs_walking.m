function [ condition_trials, condition_trials_str, condition_str ] = generate_stationary_vs_walking( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path )
% Version 2. Use turning metadata to categorize turning: 
% [ turn_t, turn_mag, counter_turn_t, counter_turn_mag ]

ac = get_analysis_constants;
settings = sensor_settings;

prestim = settings.pre_stim;
stim    = settings.stim;

STATIONARY = 1;
WALKING = 2;
condition_trials_str = {'stationary', 'walking' };
condition_str = 'stationary_vs_walking';

FWD_VELOCITY_THRESHOLD = 0.04;

trial_cnt = size( bdata_vel, 2 );
condition_trials = cell(trial_cnt,2);

fwd_vels = [];

for trial_type = 1:trial_cnt
    
    for trial_ord = 1:size( bdata_vel{trial_type}, 1 )
        
        cur_fwd_tc = bdata_vel{ trial_type }( trial_ord, ac.VEL_FWD, : );
                        
        %fwd_vel = cur_fwd_tc( find( bdata_vel_time < (prestim+stim)) );
        fwd_vel = cur_fwd_tc( find( bdata_vel_time < (prestim)) );
        
        walking_mag = mean( fwd_vel );

        fwd_vels(end+1) = walking_mag;

        %         f = figure()
        %         plot(squeeze(fwd_vel));
        %         title(['walking mag: ' num2str(walking_mag)]);
        %         waitforbuttonpress()
        %         close(f)

        if( walking_mag < FWD_VELOCITY_THRESHOLD )
            condition_trials{ trial_type, STATIONARY }( end+1 ) = trial_ord;
        else
            condition_trials{ trial_type, WALKING }( end+1 ) = trial_ord;
        end
            
    end
end

% 
% figure;
% hist(fwd_vels);

end

