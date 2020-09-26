function examine_bump_turn_dynamics_NN( full_fwd, full_yaw, full_bump, full_ephys, full_yaw_t, full_ephys_t, yaw_win_all, timebase_yaw_win )



% [ ] Compute bump jump magnitude, mean +/- SD across all trials that made it into the final dataset @Sasha R 
BUMP_JUMP_WINDOW = [0.5 1.0]; % from stim 

bump_jump_magnitudes = [];

for cond = [ 1 : size( full_bump, 1 ) ]

    for d = [ 1 : size( full_bump, 2 ) ]
        
        cur_bump_set   = full_bump{ cond, d };
        cur_bump_t_set = full_yaw_t{ cond, d };
        
        for s = [1:size(cur_bump_set,1)]
        
            cur_bump   = cur_bump_set( s, : );
            cur_bump_t = cur_bump_t_set( s, : );
            
            bump_jump_win = find( cur_bump_t >= BUMP_JUMP_WINDOW(1) & (cur_bump_t <= BUMP_JUMP_WINDOW(2) ));
        
            cur_mean_bump_jump_mag = abs(mean( cur_bump( bump_jump_win ) ));
            
            bump_jump_magnitudes( end+1 ) = cur_mean_bump_jump_mag;
        end        
    end    
end

savepath = '/data/drive2/sasha/NN_analysis/';
f = figure;

bump_jump_mag_avg = mean( bump_jump_magnitudes );
bump_jump_mag_std = std( bump_jump_magnitudes, 1 );

errorbar( 1, bump_jump_mag_avg, [ bump_jump_mag_avg - bump_jump_mag_std ], [ bump_jump_mag_avg + bump_jump_mag_std ] );
ylabel('BUmp jump magnitude');

saveas(f, [savepath 'bump_jump_magnitude.fig']);
saveas(f, [savepath 'bump_jump_magnitude.png']);


% [ ] scatterplot the signed bump jump magnitude (+ clockwise/- counterclockwise) against the subsequent signed 
% turn magnitude (again, + clockwise/- counterclockwise); we do not expect a very strong correlation, but given 
% Gaby’s comments, it seems important to examine whether there’s any relationship whatsoever (to tell him if the 
% answer to his Question 1 is A+B or just B) @Sasha R 
bump_jumps = [];
turn_post_bump_jump = [];

TURN_WINDOW = [-0.75 0];

for cond = [ 1 : size( full_bump, 1 ) ]

    for d = [ 1 : size( full_bump, 2 ) ]
        
        cur_bump_set   = full_bump{ cond, d };
        cur_bump_t_set = full_yaw_t{ cond, d };
        
        cur_yaw_in_win_set    = yaw_win_all{ cond, d };
        cur_yaw_in_win_t_set    = timebase_yaw_win{ cond, d };        
        
        for s = [ 1 : size(cur_bump_set,1) ]
        
            cur_bump   = cur_bump_set( s, : );
            cur_bump_t = cur_bump_t_set( s, : );
            
            cur_yaw_in_win = cur_yaw_in_win_set( s, : );
            cur_yaw_in_win_t = cur_yaw_in_win_t_set( s, : );
            
            % Compute avg bump (signed)
            bump_jump_win = find( cur_bump_t >= BUMP_JUMP_WINDOW(1) & (cur_bump_t <= BUMP_JUMP_WINDOW(2) ));        
            cur_mean_bump_jump = mean( cur_bump( bump_jump_win ) );            
            bump_jumps( end+1 ) = cur_mean_bump_jump;
                        
            % Compute avg turn (signed)
            cur_turn_win_idx = find( (cur_yaw_in_win_t >= TURN_WINDOW(1)) & (cur_yaw_in_win_t <= TURN_WINDOW(2) ));
            avg_turn = mean(cur_yaw_in_win( cur_turn_win_idx ));
            turn_post_bump_jump(end+1) = avg_turn;            
        end        
    end    
end

f = figure;

scatter(bump_jumps, turn_post_bump_jump);

saveas(f, [savepath 'bump_jump_vs_turn_velocity.fig']);
saveas(f, [savepath 'bump_jump_vs_turn_velocity.png']);


% [ ] Quantify probability of behavioral stopping after the ATP injection, and also the across-trial mean (+/- SD) 
% duration of stopping after the ATP injection (Question 1B) @Sasha R 

f = figure;

STOPPING_WINDOW_POST_STIM = [0.5 1.5];
STOPPING_THRESHOLD = 1.0;

stopping_event_cnt = 0;
total_event_cnt = 0;

stop_durations = [];

for cond = [ 1 : size( full_bump, 1 ) ]

    for d = [ 1 : size( full_bump, 2 ) ]
        
        cur_bump_set   = full_bump{ cond, d };
        cur_bump_t_set = full_yaw_t{ cond, d };
            
        cur_fwd_set = full_fwd{ cond, d };
        
        for s = [ 1 : size(cur_bump_set,1) ]
            cur_fwd = cur_fwd_set(s,:);
            cur_bump_t = cur_bump_t_set( s, : );
            
            stopping_win = find( cur_bump_t >= STOPPING_WINDOW_POST_STIM(1) & (cur_bump_t <= STOPPING_WINDOW_POST_STIM(2) ));        
            cur_fwd_vel_in_stopping_win = mean( cur_fwd( stopping_win ) );
            
            if( cur_fwd_vel_in_stopping_win <= STOPPING_THRESHOLD )
                stopping_event_cnt = stopping_event_cnt + 1;
                
                stop_durations(end+1) = compute_stop_duration( cur_fwd(:end), cur_bump_t(2)-cur_bump_t(1) );
                
            end
            total_event_cnt = total_event_cnt + 1;
        end
    end
end

duration_of_stop_avg = mean( stop_durations );
duration_of_stop_sd = std( stop_durations );

text(0,0, ['Probability of stopping: ' num2str(stopping_event_cnt/total_event_cnt)]);
text(0,30, ['Avg duration of stopping: ' num2str(duration_of_stop_avg) ' (SD): ' num2str(duration_of_stop_sd)]);


saveas(f, [savepath 'stopping_analysis.fig']);
saveas(f, [savepath 'stopping_analysis.png']);
end

