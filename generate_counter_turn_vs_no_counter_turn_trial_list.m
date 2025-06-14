function [ condition_trials, condition_trials_str, condition_str ] = generate_counter_turn_vs_no_counter_turn_trial_list( sid, bdata_vel_time, bdata_vel, turn_metadata, analysis_path )

ac = get_analysis_constants();
condition_trials_str = { 'counter_turns', 'no_counter_turns' };
condition_str = 'counter_turns_vs_no_counter_turning';

settings = sensor_settings;
prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

total_time = prestim + stim + poststim;

trial_type_cnt = size( bdata_vel, 2 );
condition_trials = cell(trial_type_cnt,2);

COUNTER_TURN = 1;
NO_COUNTER_TURN = 2;

FWD_VELOCITY_THRESHOLD = 0.001;

for trial_type = 1:trial_type_cnt
    
    % fourth column is the counter turn magnitude
    [Nbins, edges] = histcounts( turn_metadata{ trial_type }(:,4) );    
        
    % Take out bins near zero (and to the opposite side of the turn)
    % divide the rest of the bin space into 3 groups: small, middle, large
    % turns.
    if (trial_type == ac.LEFT )
        left_total_range = [1:(find(edges == 0) - 1)];
        
        range_len = length(left_total_range);
        range_bin = floor(range_len / 2.0);
                
        %large_turn_range = [ edges(1), edges(range_bin) ];
        %small_turn_range = [ edges(range_bin) edges(2*range_bin)];
        ct_large_turn_range = [ 0.02, 0.2 ];
        ct_small_turn_range = [ 0 0.02];

        t_turn_range = [ -0.4, -0.03 ];
    
    elseif( trial_type == ac.RIGHT )
        first_nonzero_turn_idx = (find(edges == 0) + 1);
        right_total_range = [first_nonzero_turn_idx:length(edges)];

        range_len = length(right_total_range);
        range_bin = floor(range_len / 2.0);
        
        %large_turn_range = [edges(first_nonzero_turn_idx + range_bin) edges(first_nonzero_turn_idx + 2*range_bin-1)];
        %small_turn_range = [edges(first_nonzero_turn_idx), edges(first_nonzero_turn_idx + range_bin)];        
        ct_large_turn_range = [-0.2 -0.02];
        ct_small_turn_range = [-0.02 0];
        
        t_turn_range = [ 0.02, 0.2 ];
      end    
    
    for trial_ord = 1:size( bdata_vel{trial_type}, 1 )
        
        turn_magnitude = turn_metadata{ trial_type }( trial_ord, 2 );
        counter_turn_magnitude = turn_metadata{ trial_type }( trial_ord, 4 );
        
        cur_fwd_tc = bdata_vel{ trial_type }( trial_ord, ac.VEL_FWD, : );
                        
        fwd_vel = cur_fwd_tc( find( bdata_vel_time < (prestim+stim)) );
        avg_fwd_vel = mean( fwd_vel );
                
%         if( avg_fwd_vel < FWD_VELOCITY_THRESHOLD )
%             continue;
%         end
        
        if( (turn_magnitude > t_turn_range(1)) && (turn_magnitude < t_turn_range(2)) )
            if( (counter_turn_magnitude > ct_large_turn_range(1)) && (counter_turn_magnitude < ct_large_turn_range(2)) )
                condition_trials{ trial_type, COUNTER_TURN }( end+1 ) = trial_ord;
            elseif ( (counter_turn_magnitude > ct_small_turn_range(1)) && (counter_turn_magnitude < ct_small_turn_range(2)) )
                condition_trials{ trial_type, NO_COUNTER_TURN }( end+1 ) = trial_ord;
            end
        end
    end
end

end