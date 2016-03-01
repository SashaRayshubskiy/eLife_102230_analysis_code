function [ condition_trials, condition_trials_str, condition_str ] = generate_quick_vs_delayed_counter_turn_trial_list( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path )

ac = get_analysis_constants();
condition_trials_str = { 'quick_counter_turns', 'delayed_counter_turns' };
condition_str = 'quick_vs_delayed_counter_turning';

settings = sensor_settings;
prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

total_time = prestim + stim + poststim;

trial_type_cnt = size( bdata_vel, 2 );
condition_trials = cell(trial_type_cnt,2);

QUICK_CTURN = 1;
DELAYED_CTURN = 2;

TURN_THRESHOLD = 0.03;
FWD_VELOCITY_THRESHOLD = 0.001;

for trial_type = 1:trial_type_cnt
    
    % second column is the turn magnitude
    [Nbins, edges] = histcounts( turn_metadata{ trial_type }(:,5) );    
        
    quick_cturn_range = [0.2 0.6];
    delayed_cturn_range = [0.8 1.4];
    
    for trial_ord = 1:size( bdata_vel{trial_type}, 1 )
        
        turn_magnitude = turn_metadata{ trial_type }( trial_ord, 2 );
        
        turn_delay = turn_metadata{ trial_type }( trial_ord, 5 );
        
        cur_fwd_tc = bdata_vel{ trial_type }( trial_ord, ac.VEL_FWD, : );
                        
        fwd_vel = cur_fwd_tc( find( bdata_vel_time < (prestim+stim)) );
        avg_fwd_vel = mean( fwd_vel );
                
        if( avg_fwd_vel < FWD_VELOCITY_THRESHOLD )
            continue;
        end
        
        if( abs(turn_magnitude) < TURN_THRESHOLD )
            continue;
        end
        
        if( (turn_delay > quick_cturn_range(1)) && (turn_delay < quick_cturn_range(2)) )
            condition_trials{ trial_type, QUICK_CTURN }( end+1 ) = trial_ord;
        elseif ( (turn_delay > delayed_cturn_range(1)) && (turn_delay < delayed_cturn_range(2)) )
            condition_trials{ trial_type, DELAYED_CTURN }( end+1 ) = trial_ord;        
        end
    end
end

end