function [ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list_v3( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path )
% Version 2. Use turning metadata to categorize turning: 
% [ turn_t, turn_mag, counter_turn_t, counter_turn_mag ]

%%% THIS IS A TEMP 

ac = get_analysis_constants;
settings = sensor_settings;

prestim = settings.pre_stim;
stim    = settings.stim;

EXPECTED = 1;
IGNORED = 2;
condition_trials_str = {'expected_turning', 'ignoring_stim' };
condition_str = 'expected_vs_ignoring_turning';

LEFT_TURN  = 1;
RIGHT_TURN = 2;
NO_TURN    = 3;

TURN_THRESHOLD = 0.02; % ???
FWD_VELOCITY_THRESHOLD = 0.001;
%TURN_THRESHOLD = 0.0000001; % ???
%FWD_VELOCITY_THRESHOLD = 0.00000001;

trial_cnt = size( bdata_vel, 2 );
condition_trials = cell(trial_cnt,2);

for trial_type = 1:trial_cnt
    
    [Nbins, edges] = histcounts(turn_metadata{ trial_type }(:,2));    
        
    % Set the cutoff to the left of the first bin for
    turn_cutoff = TURN_THRESHOLD;
    if( trial_type == ac.LEFT )
        turn_cutoff = -1.0 * TURN_THRESHOLD;
        %turn_cutoff = edges(find(edges == 0) - 1);
    elseif( ( trial_type == ac.RIGHT ) || ( trial_type == ac.BOTH ))
        turn_cutoff = TURN_THRESHOLD;
        %turn_cutoff = edges(find(edges == 0) + 1);        
    end
    
    avg_turn_mag = mean(turn_metadata{ trial_type }( :, 2 ));
    
    for trial_ord = 1:size( bdata_vel{trial_type}, 1 )
        
        cur_avg_yaw = mean(bdata_vel{ trial_type }( trial_ord, ac.VEL_YAW, : ));
        
        %turn_magnitude = turn_metadata{ trial_type }( trial_ord, 2 ) - cur_avg_yaw;
        turn_magnitude = turn_metadata{ trial_type }( trial_ord, 2 );
        
        disp(['turn_mag: ' num2str(turn_magnitude) ' avg_turn_mag: ' num2str(cur_avg_yaw)]);
        
        turn_t = turn_metadata{ trial_type }( trial_ord, 1 );
        
        cur_fwd_tc = bdata_vel{ trial_type }( trial_ord, ac.VEL_FWD, : );
                        
        fwd_vel = cur_fwd_tc( find( bdata_vel_time < (prestim+stim)) );
        avg_fwd_vel = mean( fwd_vel );
        
        if 0
        if( abs(avg_fwd_vel) < FWD_VELOCITY_THRESHOLD )
            continue;
        end
        
%         if( turn_t < 3.1 )
%             continue;
%         end
        end
        
        % yaw_during_stim = cur_yaw_tc( find( (bdata_vel_time > prestim) & (bdata_vel_time <= (prestim+stim))) );
        % avg_yaw_during_stim = mean(yaw_during_stim);
                   
        if( trial_type == ac.LEFT )
            if( turn_magnitude < turn_cutoff )
                condition_trials{ trial_type, EXPECTED }( end+1 ) = trial_ord;
            else
                condition_trials{ trial_type, IGNORED }( end+1 ) = trial_ord;
            end
        elseif( ( trial_type == ac.RIGHT ) || ( trial_type == ac.BOTH ))
            if( turn_magnitude > turn_cutoff )
                condition_trials{ trial_type, EXPECTED }( end+1 ) = trial_ord;
            else
                condition_trials{ trial_type, IGNORED }( end+1 ) = trial_ord;
            end        
        end
    end
end

end

