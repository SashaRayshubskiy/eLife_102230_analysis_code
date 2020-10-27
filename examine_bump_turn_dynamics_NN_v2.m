function examine_bump_turn_dynamics_NN_v2( full_fwd, full_yaw, full_bump, full_bump_t, full_ephys, full_yaw_t, full_ephys_t, yaw_win_all, timebase_yaw_win, non_return_bump_jump, non_return_yaw_vel_near_max_bump_movement )

savepath = '/data/drive2/sasha/NN_analysis/';

BUMP_JUMP_WINDOW = [0.75 1.25]; % from stim 
BUMP_TO_DEG_CONVERSION_FACTOR = 45;

% [ ] scatterplot the signed bump jump magnitude (+ clockwise/- counterclockwise) against the subsequent signed 
% turn magnitude (again, + clockwise/- counterclockwise); we do not expect a very strong correlation, but given 
% Gaby’s comments, it seems important to examine whether there’s any relationship whatsoever (to tell him if the 
% answer to his Question 1 is A+B or just B) @Sasha R 
bump_jumps = [];
turn_post_bump_jump = [];

TURN_WINDOW = [ -0.65 0.0 ];
cnt = 0;
cnt2 = 0;

neg_bump_jump = [];
neg_vel       = [];
pos_bump_jump = [];
pos_vel       = [];

for cond = [ 1 : size( full_bump, 1 ) ]

    for d = [ 1 : size( full_bump, 2 ) ]
        
        cur_bump_set   = full_bump{ cond, d };
        cur_bump_t_set = full_bump_t{ cond, d };
        
        cur_yaw_in_win_set    = yaw_win_all{ cond, d };
        cur_yaw_in_win_t  = timebase_yaw_win{ cond, d };        
        
        for s = [ 1 : size(cur_bump_set,1) ]
        
            cur_bump   = cur_bump_set( s, : );
            cur_bump_t = cur_bump_t_set( s, : );
            
            cur_yaw_in_win = cur_yaw_in_win_set( s, : );
                        
            % Compute avg bump (signed)
            bump_jump_win = find( cur_bump_t >= BUMP_JUMP_WINDOW(1) & (cur_bump_t <= BUMP_JUMP_WINDOW(2) ));        
            cur_mean_bump_jump = mean( cur_bump( bump_jump_win ) );            
            bump_jumps( end+1 ) = BUMP_TO_DEG_CONVERSION_FACTOR*cur_mean_bump_jump;
            
            % Compute avg turn (signed)st
            cur_turn_win_idx = find( (cur_yaw_in_win_t >= TURN_WINDOW(1)) & (cur_yaw_in_win_t <= TURN_WINDOW(2) ));
            avg_turn = mean(cur_yaw_in_win( cur_turn_win_idx ));
            turn_post_bump_jump(end+1) = avg_turn;            
                        
            if( bump_jumps(end) < 0 )
                neg_bump_jump(end+1) = bump_jumps(end);
                neg_vel(end+1) = avg_turn;
            else
                pos_bump_jump(end+1) = bump_jumps(end);
                pos_vel(end+1) = avg_turn;
            end
            
            if( (bump_jumps(end) < 0) && (turn_post_bump_jump(end) > 0))
                cnt = cnt + 1;
            end
            if( (bump_jumps(end) > 0) && (turn_post_bump_jump(end) < 0))
                cnt2 = cnt2 + 1;
            end
        end        
    end    
end

disp( [ num2str(cnt) ' ' num2str(cnt2) ] );

[ h, p ] = ttest2( neg_vel, pos_vel ); 

disp(['yaw p-val neg vs. pos bump jump - returning trials: ' num2str(p)]);

non_return_bump_jump = non_return_bump_jump*BUMP_TO_DEG_CONVERSION_FACTOR;

neg_vel_2 = [];
pos_vel_2 = [];
neg_bump_jump_2 = [];
pos_bump_jump_2 = [];
for i = 1:length( non_return_bump_jump )
    if( non_return_bump_jump(i) < 0 )
        neg_vel_2(end+1) = non_return_yaw_vel_near_max_bump_movement( i );
        neg_bump_jump_2(end+1) = non_return_bump_jump(i);
    else
        pos_vel_2(end+1) = non_return_yaw_vel_near_max_bump_movement( i );
        pos_bump_jump_2(end+1) = non_return_bump_jump(i);
    end
end

[ h, p_non_return ] = ttest2( neg_vel_2, pos_vel_2 ); 

disp( ['yaw p-val neg vs. pos bump jump - non returning trials: ' num2str(p_non_return)] );

f = figure;

hold on;
scatter( bump_jumps, turn_post_bump_jump, 'DisplayName', ['yaw p-val neg vs. pos bump jump - returning trials: ' num2str(p)] );
scatter( non_return_bump_jump, non_return_yaw_vel_near_max_bump_movement, 'x', 'DisplayName', ['yaw p-val neg vs. pos bump jump - non returning trials: ' num2str(p_non_return)] );

% Returned trials
neg_bump_jump_avg = mean( neg_bump_jump );
neg_bump_jump_sem = std( neg_bump_jump, 1 ) ./ sqrt( length( neg_bump_jump ) );

pos_bump_jump_avg = mean( pos_bump_jump );
pos_bump_jump_sem = std( pos_bump_jump, 1 ) ./ sqrt( length( pos_bump_jump ) );

neg_vel_avg = mean( neg_vel );
neg_vel_sem = std( neg_vel, 1 ) ./ sqrt( length( neg_vel ) );

pos_vel_avg = mean( pos_vel );
pos_vel_sem = std( pos_vel, 1 ) ./ sqrt( length( pos_vel ) );

y_lower = [ neg_vel_sem, pos_vel_sem ];
y_upper = [ neg_vel_sem, pos_vel_sem ];

x_left = [ neg_bump_jump_sem, pos_bump_jump_sem ];
x_right = [ neg_bump_jump_sem, pos_bump_jump_sem ];

x = [ neg_bump_jump_avg, pos_bump_jump_avg ];
y = [ neg_vel_avg, pos_vel_avg ];

hold on;
errorbar( x, y, y_lower, y_upper, x_left, x_right, 'o', 'DisplayName', 'returning trials' );

% Non Returned trials
neg_bump_jump_avg = mean( neg_bump_jump_2 );
neg_bump_jump_sem = std( neg_bump_jump_2, 1 ) ./ sqrt( length( neg_bump_jump_2 ) );

pos_bump_jump_avg = mean( pos_bump_jump_2 );
pos_bump_jump_sem = std( pos_bump_jump_2, 1 ) ./ sqrt( length( pos_bump_jump_2 ) );

neg_vel_avg = mean( neg_vel_2 );
neg_vel_sem = std( neg_vel_2, 1 ) ./ sqrt( length( neg_vel_2 ) );

pos_vel_avg = mean( pos_vel_2 );
pos_vel_sem = std( pos_vel_2, 1 ) ./ sqrt( length( pos_vel_2 ) );

y_lower = [ neg_vel_sem, pos_vel_sem ];
y_upper = [ neg_vel_sem, pos_vel_sem ];

x_left = [ neg_bump_jump_sem, pos_bump_jump_sem ];
x_right = [ neg_bump_jump_sem, pos_bump_jump_sem ];

x = [ neg_bump_jump_avg, pos_bump_jump_avg ];
y = [ neg_vel_avg, pos_vel_avg ];

errorbar( x, y, y_lower, y_upper, x_left, x_right, 'o', 'DisplayName', 'non returning trials' );

legend();

if 0
plot( x, y, 'x' );

plot( [ x(1) - neg_bump_jump_sem, x(1)], ...
      [ y(1), y(1) ] );
plot( [ x(1) + neg_bump_jump_sem, x(1)], ...
      [ y(1), y(1) ] );

plot( [ x(1), x(1)], ...
      [ y(1) - neg_vel_sem, y(1) ] );
plot( [ x(1), x(1)], ...
      [ y(1) + neg_vel_sem, y(1) ] );

plot( [ x(2) - pos_bump_jump_sem, x(2)], ...
      [ y(2), y(2) ] );
plot( [ x(2) + pos_bump_jump_sem, x(2)], ...
      [ y(2), y(2) ] );

plot( [ x(2), x(2)], ...
      [ y(2) - pos_vel_sem, y(2) ] );
plot( [ x(2), x(2)], ...
      [ y(2) + pos_vel_sem, y(2) ] );

end 
xlabel('Bump jump (deg)');
ylabel('Yaw post bump jump (deg/s)');
xlim([-200 200]);
ylim([-300 300]);
set( gca(), 'TickDir', 'out' );

grid on;

saveas(f, [savepath 'bump_jump_vs_turn_velocity_v3.fig']);
saveas(f, [savepath 'bump_jump_vs_turn_velocity_v3.png']);
saveas(f, [savepath 'bump_jump_vs_turn_velocity_v3.svg']);

end
