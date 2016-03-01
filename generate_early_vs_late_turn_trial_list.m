function [ condition_trials, condition_trials_str, condition_str ] = generate_early_vs_late_turn_trial_list( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path )

ac = get_analysis_constants();
condition_trials_str = { 'early_turns', 'late_turns' };
condition_str = 'early_vs_late_turning';

settings = sensor_settings;
prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

total_time = prestim + stim + poststim;

trial_type_cnt = size( bdata_vel, 2 );
condition_trials = cell(trial_type_cnt,2);

EARLY_TURN = 1;
LATE_TURN = 2;

TURN_THRESHOLD = 0.03;
FWD_VELOCITY_THRESHOLD = 0.001;

for trial_type = 1:trial_type_cnt
    
    % second column is the turn magnitude
    [Nbins, edges] = histcounts( turn_metadata{ trial_type }(:,1) );    
        
    early_turn_range = [ 3.1, 3.4 ];
    late_turn_range =  [ 3.4, 4.0 ];    
    
    for trial_ord = 1:size( bdata_vel{trial_type}, 1 )
        
        turn_t = turn_metadata{ trial_type }( trial_ord, 1 );
        turn_magnitude = turn_metadata{ trial_type }( trial_ord, 2 );
        
        cur_fwd_tc = bdata_vel{ trial_type }( trial_ord, ac.VEL_FWD, : );
                        
        fwd_vel = cur_fwd_tc( find( bdata_vel_time < (prestim+stim)) );
        avg_fwd_vel = mean( fwd_vel );
                
        if( avg_fwd_vel < FWD_VELOCITY_THRESHOLD )
            continue;
        end
        
        if( abs(turn_magnitude) < TURN_THRESHOLD )
            continue;
        end
        
        if( (turn_t > early_turn_range(1)) && (turn_t < early_turn_range(2)) )
            condition_trials{ trial_type, EARLY_TURN }( end+1 ) = trial_ord;
        elseif ( (turn_t > late_turn_range(1)) && (turn_t < late_turn_range(2)) )
            condition_trials{ trial_type, LATE_TURN }( end+1 ) = trial_ord;        
        end
    end
end

end