function examine_bump_turn_dynamics_NN( full_fwd, full_yaw, full_bump, full_bump_t, full_ephys, full_yaw_t, full_ephys_t, yaw_win_all, timebase_yaw_win )

savepath = '/data/drive2/sasha/NN_analysis/';

% [ ] Compute bump jump magnitude, mean +/- SD across all trials that made it into the final dataset @Sasha R 
BUMP_JUMP_WINDOW = [0.75 1.25]; % from stim 
BUMP_TO_DEG_CONVERSION_FACTOR = 45;

bump_jump_magnitudes = [];
f = figure;
for cond = [ 1 : size( full_bump, 1 ) ]

    for d = [ 1 : size( full_bump, 2 ) ]
        
        cur_bump_set   = full_bump{ cond, d };
        cur_bump_t_set = full_bump_t{ cond, d };
        
        for s = [1:size(cur_bump_set,1)]
        
            cur_bump   = BUMP_TO_DEG_CONVERSION_FACTOR*cur_bump_set( s, : );
            cur_bump_t = cur_bump_t_set( s, : );
            
            bump_jump_win = find( cur_bump_t >= BUMP_JUMP_WINDOW(1) & (cur_bump_t <= BUMP_JUMP_WINDOW(2) ));
        
            hold on; 
            plot( cur_bump_t, cur_bump );
            
            % disp(size(bump_jump_win));
            cur_mean_bump_jump_mag = abs(mean( cur_bump( bump_jump_win ) ));
            
            bump_jump_magnitudes( end+1 ) = cur_mean_bump_jump_mag;
        end        
    end    
end

saveas(f, [savepath 'bump_jumps.fig']);
saveas(f, [savepath 'bump_jumps.png']);

f = figure;

bump_jump_mag_avg = mean( bump_jump_magnitudes );
bump_jump_mag_std = std( bump_jump_magnitudes, 1 );

text(0,0.5, ['Bump jump magnitude: ' num2str(bump_jump_mag_avg) ' (SD): ' num2str(bump_jump_mag_std) ' N trials: ' num2str(length(bump_jump_magnitudes))]);

saveas(f, [savepath 'bump_jump_magnitude.fig']);
saveas(f, [savepath 'bump_jump_magnitude.png']);


% [ ] scatterplot the signed bump jump magnitude (+ clockwise/- counterclockwise) against the subsequent signed 
% turn magnitude (again, + clockwise/- counterclockwise); we do not expect a very strong correlation, but given 
% Gaby’s comments, it seems important to examine whether there’s any relationship whatsoever (to tell him if the 
% answer to his Question 1 is A+B or just B) @Sasha R 
bump_jumps = [];
turn_post_bump_jump = [];

TURN_WINDOW = [ -0.65 0.0 ];
cnt = 0;
cnt2 = 0;

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
            
            % Compute avg turn (signed)
            cur_turn_win_idx = find( (cur_yaw_in_win_t >= TURN_WINDOW(1)) & (cur_yaw_in_win_t <= TURN_WINDOW(2) ));
            avg_turn = mean(cur_yaw_in_win( cur_turn_win_idx ));
            turn_post_bump_jump(end+1) = avg_turn;            

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

f = figure;

scatter(bump_jumps, turn_post_bump_jump);
xlabel('Bump jump (deg)');
ylabel('Yaw post bump jump (deg/s)');
xlim([-200 200]);
ylim([-300 300]);
grid on;

saveas(f, [savepath 'bump_jump_vs_turn_velocity.fig']);
saveas(f, [savepath 'bump_jump_vs_turn_velocity.png']);

% Here's an alternative suggestion for an analysis that should show us the effect of the stimulus on forward 
% walking (if indeed there is any effect at all): separate trials into two sets, (1) the fly is walking prior 
% to stimulus onset, and (2) the fly is not walking prior to stimulus onset within each set, align trials to 
% the stimulus and compute the average forward velocity versus time if it looks like there is potentially an 
% effect of the stimulus, then we can follow up with statistical testing; here I would suggest taking the mean 
% forward velocity over a pre-stim window and a post-stim window for each trial, and then comparing them via a 
% paired t-test or some other paired test (n=# of trials)

STOPPING_WINDOW_POST_STIM = [ 0.5 1.0 ];
STOPPING_WINDOW_PRE_STIM = [ -0.5 0.0 ];
% STOPPING_THRESHOLD = 3.0; % total speed mm/s
STOPPING_THRESHOLD = 0.5; % mm/s
YAW_STOPPING_THRESHOLD = 20; % deg/s
RAW_FWD_ZERO_OFFSET = 0.05;

pre_stim_stopped_fwd = [];
pre_stim_moving_fwd = [];

post_stim_stopped_fwd = [];
post_stim_moving_fwd = [];

idx = 1;
for cond = [ 1 : size( full_bump, 1 ) ]

    for d = [ 1 : size( full_bump, 2 ) ]
                    
        cur_fwd_set = full_fwd{ cond, d };
        cur_yaw_set = full_yaw{ cond, d };
        cur_fwd_t_set = full_yaw_t{ cond, d };
                
        for s = [ 1 : size(cur_fwd_set,1) ]
            cur_fwd   = abs(convert_fwd_to_mms(cur_fwd_set(s,:) - RAW_FWD_ZERO_OFFSET));
            cut_fwd_t = cur_fwd_t_set( s, : );
            cur_yaw   = cur_yaw_set(s,:);
                                        
            yaw_corrected_zero_offset = correct_yaw_zero_offsets( cur_yaw, idx );
            idx = idx + 1;
            
            pre_stim_win = find( cut_fwd_t >= STOPPING_WINDOW_PRE_STIM(1) & (cut_fwd_t <= STOPPING_WINDOW_PRE_STIM(2) )); 
            
            avg_fwd_pre_stim = mean( cur_fwd( pre_stim_win ) );
            
            if( avg_fwd_pre_stim <= STOPPING_THRESHOLD )
                pre_stim_stopped_fwd(end+1,:) = cur_fwd;                
            else
                pre_stim_moving_fwd(end+1,:) = cur_fwd;       
                
                post_stim_win = find( cut_fwd_t >= STOPPING_WINDOW_POST_STIM(1) & (cut_fwd_t <= STOPPING_WINDOW_POST_STIM(2) )); 
                    
                avg_fwd_pre_stim = mean( cur_fwd( post_stim_win ) );
                
                if( avg_fwd_pre_stim <= STOPPING_THRESHOLD )
                    post_stim_stopped_fwd(end+1, :) = cur_fwd;
                else
                    post_stim_moving_fwd(end+1, :) = cur_fwd;
                end
            end
        end
    end
end

f = figure;

subplot(2,1,1);
hold on;
plot( cut_fwd_t, mean( pre_stim_stopped_fwd ) );
plot( cut_fwd_t, mean( pre_stim_moving_fwd ) );
legend([ {['pre stim not moving ' num2str(size(pre_stim_stopped_fwd,1)) ' trials']}, ...
         {['pre stim moving ' num2str(size(pre_stim_moving_fwd,1)) ' trials']}]);
ylabel('Fwd vel (mm/s)');
xlabel('Time (s)');
ylim([0.0 15.0]);
title('All trials');

subplot(2,1,2);
hold on;
plot( cut_fwd_t, mean( post_stim_stopped_fwd ) );
plot( cut_fwd_t, mean( post_stim_moving_fwd ) );
legend([ {['Post stim not moving ' num2str(size(post_stim_stopped_fwd,1)) ' trials']}, ...
         {['Post stim moving ' num2str(size(post_stim_moving_fwd,1)) ' trials']}]);
ylabel('Fwd vel (mm/s)');
xlabel('Time (s)');
ylim([0.0 15.0]);
title('Pre-stim moving trials only');

saveas(f, [savepath 'moving_vs_not_moving_pre_stim_fwd_analysis.fig']);
saveas(f, [savepath 'moving_vs_not_moving_pre_stim_fwd_analysis.png']);

% Same analysis as above, only on yaw

STOPPING_WINDOW_POST_STIM = [ 0.5 1.0 ];
STOPPING_WINDOW_PRE_STIM = [ -0.5 0.0 ];
% STOPPING_THRESHOLD = 3.0; % total speed mm/s
YAW_STOPPING_THRESHOLD = 15; % deg/s

pre_stim_stopped_yaw = [];
pre_stim_moving_yaw = [];

post_stim_stopped_yaw = [];
post_stim_moving_yaw = [];

idx = 1;
for cond = [ 1 : size( full_bump, 1 ) ]

    for d = [ 1 : size( full_bump, 2 ) ]
                    
        %cur_fwd_set = full_fwd{ cond, d };
        cur_yaw_set = full_yaw{ cond, d };
        cur_fwd_t_set = full_yaw_t{ cond, d };
                
        for s = [ 1 : size(cur_yaw_set,1) ]
            %cur_fwd   = abs(convert_fwd_to_mms(cur_fwd_set(s,:) - RAW_FWD_ZERO_OFFSET));
            cut_fwd_t = cur_fwd_t_set( s, : );
            cur_yaw_tmp   = cur_yaw_set(s,:);
                                        
            cur_yaw = abs(correct_yaw_zero_offsets( cur_yaw_tmp, idx ));
            idx = idx + 1;
            
            pre_stim_win = find( cut_fwd_t >= STOPPING_WINDOW_PRE_STIM(1) & (cut_fwd_t <= STOPPING_WINDOW_PRE_STIM(2) )); 
            
            avg_yaw_pre_stim = mean( cur_yaw( pre_stim_win ) );
            
            if( avg_yaw_pre_stim <= YAW_STOPPING_THRESHOLD )
                pre_stim_stopped_yaw(end+1,:) = cur_yaw;                
            else
                pre_stim_moving_yaw(end+1,:) = cur_yaw;       
                
                post_stim_win = find( cut_fwd_t >= STOPPING_WINDOW_POST_STIM(1) & (cut_fwd_t <= STOPPING_WINDOW_POST_STIM(2) )); 
                    
                avg_yaw_post_stim = mean( cur_yaw( post_stim_win ) );
                
                if( avg_yaw_post_stim <= YAW_STOPPING_THRESHOLD )
                    post_stim_stopped_yaw(end+1, :) = cur_yaw;
                else
                    post_stim_moving_yaw(end+1, :) = cur_yaw;
                end
            end
        end
    end
end

f = figure;

YAW_SPEED_LIM = 250.0;

subplot(2,1,1);
hold on;
plot( cut_fwd_t, mean( pre_stim_stopped_yaw ) );
plot( cut_fwd_t, mean( pre_stim_moving_yaw ) );
legend([ {['pre stim not moving ' num2str(size(pre_stim_stopped_yaw,1)) ' trials']}, ...
         {['pre stim moving ' num2str(size(pre_stim_moving_yaw,1)) ' trials']}]);
ylabel('Yaw speed (deg/s)');
xlabel('Time (s)');
ylim([0.0 YAW_SPEED_LIM]);
title('All trials');

subplot(2,1,2);
hold on;
plot( cut_fwd_t, mean( post_stim_stopped_yaw ) );
plot( cut_fwd_t, mean( post_stim_moving_yaw ) );
legend([ {['Post stim not moving ' num2str( size( post_stim_stopped_yaw, 1 )) ' trials']}, ...
         {['Post stim moving ' num2str( size( post_stim_moving_yaw, 1 )) ' trials']}]);
ylabel('Yaw speed (deg/s)');
xlabel('Time (s)');
ylim([0.0 YAW_SPEED_LIM]);
title('Pre-stim moving trials only');

saveas(f, [savepath 'moving_vs_not_moving_pre_stim_yaw_analysis.fig']);
saveas(f, [savepath 'moving_vs_not_moving_pre_stim_yaw_analysis.png']);


STOPPING_WINDOW_POST_STIM = [ 0.5 1.0 ];
STOPPING_WINDOW_PRE_STIM = [ -0.5 0.0 ];
% STOPPING_THRESHOLD = 3.0; % total speed mm/s
YAW_STOPPING_THRESHOLD = 15; % deg/s

pre_stim_stopped_yaw = [];
pre_stim_moving_yaw = [];

post_stim_stopped_yaw = [];
post_stim_moving_yaw = [];

idx = 1;
for cond = [ 1 : size( full_bump, 1 ) ]

    for d = [ 1 : size( full_bump, 2 ) ]
                    
        cur_fwd_set = full_fwd{ cond, d };
        cur_yaw_set = full_yaw{ cond, d };
        cur_fwd_t_set = full_yaw_t{ cond, d };
                
        for s = [ 1 : size(cur_yaw_set,1) ]
            cur_fwd   = abs(convert_fwd_to_mms(cur_fwd_set(s,:) - RAW_FWD_ZERO_OFFSET));
            cut_fwd_t = cur_fwd_t_set( s, : );
            cur_yaw_tmp   = cur_yaw_set(s,:);
                                        
            cur_yaw = abs(correct_yaw_zero_offsets( cur_yaw_tmp, idx ));
            idx = idx + 1;
            
            pre_stim_win = find( cut_fwd_t >= STOPPING_WINDOW_PRE_STIM(1) & (cut_fwd_t <= STOPPING_WINDOW_PRE_STIM(2) )); 
            
            avg_fwd_pre_stim = mean( cur_fwd( pre_stim_win ) );
            
            if( avg_fwd_pre_stim <= STOPPING_THRESHOLD )
                pre_stim_stopped_yaw(end+1,:) = cur_yaw;                
            else
                pre_stim_moving_yaw(end+1,:) = cur_yaw;       
                
                post_stim_win = find( cut_fwd_t >= STOPPING_WINDOW_POST_STIM(1) & (cut_fwd_t <= STOPPING_WINDOW_POST_STIM(2) )); 
                    
                avg_fwd_post_stim = mean( cur_fwd( post_stim_win ) );
                
                if( avg_fwd_post_stim <= STOPPING_THRESHOLD )
                    post_stim_stopped_yaw(end+1, :) = cur_yaw;
                else
                    post_stim_moving_yaw(end+1, :) = cur_yaw;
                end
            end
        end
    end
end

f = figure;

YAW_SPEED_LIM = 250.0;

subplot(2,1,1);
hold on;
plot( cut_fwd_t, mean( pre_stim_stopped_yaw ) );
plot( cut_fwd_t, mean( pre_stim_moving_yaw ) );
legend([ {['pre stim not moving ' num2str(size(pre_stim_stopped_yaw,1)) ' trials']}, ...
         {['pre stim moving ' num2str(size(pre_stim_moving_yaw,1)) ' trials']}]);
ylabel('Yaw speed (deg/s)');
xlabel('Time (s)');
ylim([0.0 YAW_SPEED_LIM]);
title('All trials, use fwd to classify movement');

subplot(2,1,2);
hold on;
plot( cut_fwd_t, mean( post_stim_stopped_yaw ) );
plot( cut_fwd_t, mean( post_stim_moving_yaw ) );
legend([ {['Post stim not moving ' num2str( size( post_stim_stopped_yaw, 1 )) ' trials']}, ...
         {['Post stim moving ' num2str( size( post_stim_moving_yaw, 1 )) ' trials']}]);
ylabel('Yaw speed (deg/s)');
xlabel('Time (s)');
ylim([0.0 YAW_SPEED_LIM]);
title('Pre-stim moving trials only');

saveas(f, [savepath 'moving_vs_not_moving_pre_stim_fwd_to_yaw_analysis.fig']);
saveas(f, [savepath 'moving_vs_not_moving_pre_stim_fwd_to_yaw_analysis.png']);


if 0
% [ ] Quantify probability of behavioral stopping after the ATP injection, and also the across-trial mean (+/- SD) 
% duration of stopping after the ATP injection (Question 1B) @Sasha R 

STOPPING_WINDOW_POST_STIM = [0.5 1.5];
% STOPPING_THRESHOLD = 3.0; % total speed mm/s
STOPPING_THRESHOLD = 0.7; % mm/s

RAW_FWD_ZERO_OFFSET = 0.05;
RAW_YAW_ZERO_OFFSET = 24.0;

DEG_TO_MMS_FOR_YAW = 0.0553;

stopping_event_cnt = 0;
total_event_cnt = 0;

stop_durations = [];

total_speeds = [];

bump_jump_magnitudes_stopped = [];
bump_jump_magnitudes_moving = [];

idx = 1;
for cond = [ 1 : size( full_bump, 1 ) ]

    for d = [ 1 : size( full_bump, 2 ) ]
        
        cur_bump_set   = full_bump{ cond, d };
        cur_bump_t_set = full_bump_t{ cond, d };
            
        cur_fwd_set = full_fwd{ cond, d };
        cur_yaw_set = full_yaw{ cond, d };
                
        for s = [ 1 : size(cur_bump_set,1) ]
            cur_fwd = convert_fwd_to_mms(cur_fwd_set(s,:) - RAW_FWD_ZERO_OFFSET);
           
            
            cur_yaw = cur_yaw_set( s, :);
                
%             figure;
%             yyaxis right;
%             plot( cur_fwd );
%             yyaxis left;
%             plot( cur_yaw );
%             title( ['Idx: ' num2str(idx) ] );
            
            yaw_corrected_zero_offset = correct_yaw_zero_offsets( cur_yaw, idx );
            idx = idx + 1;
            
            cur_yaw_corr = DEG_TO_MMS_FOR_YAW * yaw_corrected_zero_offset;
            
            % total_speed = abs(cur_fwd) + abs(cur_yaw_corr);
            % total_speed = abs(cur_fwd);
            total_speed = cur_fwd;
                 
            stopping_win = find( cur_bump_t >= STOPPING_WINDOW_POST_STIM(1) & (cur_bump_t <= STOPPING_WINDOW_POST_STIM(2) ));        
            cur_fwd_vel_in_stopping_win = mean( total_speed( stopping_win ) );
            
            total_speeds(end+1) = cur_fwd_vel_in_stopping_win;
            
            cur_bump   = BUMP_TO_DEG_CONVERSION_FACTOR*cur_bump_set( s, : );
            cur_bump_t = cur_bump_t_set( s, : );
                        
            figure;
            hold on;
            yyaxis right;
            hold on; 
            plot( cur_bump_t, cur_bump );
            yyaxis left;
            plot( full_yaw_t{1,1}(1,:), cur_fwd );
            
            bump_jump_win = find( cur_bump_t >= BUMP_JUMP_WINDOW(1) & (cur_bump_t <= BUMP_JUMP_WINDOW(2) ));
                    
            % disp(size(bump_jump_win));
            cur_mean_bump_jump_mag = abs(mean( cur_bump( bump_jump_win ) ));
            
            if( cur_fwd_vel_in_stopping_win <= STOPPING_THRESHOLD )
                stopping_event_cnt = stopping_event_cnt + 1;
                
                % t=0 stim trigger
                stopping_win = find( cur_bump_t >= 0 );
                stop_durations(end+1) = compute_stop_duration( total_speed(stopping_win), cur_bump_t(2)-cur_bump_t(1), STOPPING_THRESHOLD );
                
                bump_jump_magnitudes_stopped( end+1 ) = cur_mean_bump_jump_mag;
            else
                bump_jump_magnitudes_moving( end+1 ) = cur_mean_bump_jump_mag;
            end
            
            total_event_cnt = total_event_cnt + 1;
        end
    end
end

duration_of_stop_avg = mean( stop_durations );
duration_of_stop_sd = std( stop_durations );

f = figure;
hold on;
text( 0, 0.2, ['Probability of stopping: ' num2str(stopping_event_cnt/total_event_cnt) ' N trials: ' num2str(total_event_cnt)]);
text( 0, 0.5, ['Avg duration of stopping: ' num2str(duration_of_stop_avg) ' (SD): ' num2str(duration_of_stop_sd) ' N trials: ' num2str(total_event_cnt)]);

avg_bump_mag_moving = mean( bump_jump_magnitudes_moving );
std_bump_mag_moving = std( bump_jump_magnitudes_moving );

avg_bump_mag_stopped = mean( bump_jump_magnitudes_stopped );
std_bump_mag_stopped = std( bump_jump_magnitudes_stopped );

text( 0, 0.8, ['Bump mag moving: ' num2str( avg_bump_mag_moving ) ' (SD): ' num2str( std_bump_mag_moving ) ' N trials: ' num2str(length(bump_jump_magnitudes_moving))]);
text( 0, 1.0, ['Bump mag stopped: ' num2str( avg_bump_mag_stopped ) ' (SD): ' num2str( std_bump_mag_stopped ) ' N trials: ' num2str(length(bump_jump_magnitudes_stopped))]);

ylim([0 1.4]);

saveas(f, [savepath 'stopping_analysis.fig']);
saveas(f, [savepath 'stopping_analysis.png']);
end

end

% ATTIC 
% [ ] Quantify probability of behavioral stopping after the ATP injection, and also the across-trial mean (+/- SD) 
% duration of stopping after the ATP injection (Question 1B) @Sasha R 
% 
% STOPPING_WINDOW_POST_STIM = [0.5 1.5];
% % STOPPING_THRESHOLD = 3.0; % total speed mm/s
% STOPPING_THRESHOLD = 0.7; % mm/s
% 
% RAW_FWD_ZERO_OFFSET = 0.05;
% RAW_YAW_ZERO_OFFSET = 24.0;
% 
% DEG_TO_MMS_FOR_YAW = 0.0553;
% 
% stopping_event_cnt = 0;
% total_event_cnt = 0;
% 
% stop_durations = [];
% 
% total_speeds = [];
% 
% bump_jump_magnitudes_stopped = [];
% bump_jump_magnitudes_moving = [];
% 
% idx = 1;
% for cond = [ 1 : size( full_bump, 1 ) ]
% 
%     for d = [ 1 : size( full_bump, 2 ) ]
%         
%         cur_bump_set   = full_bump{ cond, d };
%         cur_bump_t_set = full_bump_t{ cond, d };
%             
%         cur_fwd_set = full_fwd{ cond, d };
%         cur_yaw_set = full_yaw{ cond, d };
%                 
%         for s = [ 1 : size(cur_bump_set,1) ]
%             cur_fwd = convert_fwd_to_mms(cur_fwd_set(s,:) - RAW_FWD_ZERO_OFFSET);
%            
%             
%             cur_yaw = cur_yaw_set( s, :);
%                 
% %             figure;
% %             yyaxis right;
% %             plot( cur_fwd );
% %             yyaxis left;
% %             plot( cur_yaw );
% %             title( ['Idx: ' num2str(idx) ] );
%             
%             yaw_corrected_zero_offset = correct_yaw_zero_offsets( cur_yaw, idx );
%             idx = idx + 1;
%             
%             cur_yaw_corr = DEG_TO_MMS_FOR_YAW * yaw_corrected_zero_offset;
%             
%             % total_speed = abs(cur_fwd) + abs(cur_yaw_corr);
%             % total_speed = abs(cur_fwd);
%             total_speed = cur_fwd;
%                  
%             stopping_win = find( cur_bump_t >= STOPPING_WINDOW_POST_STIM(1) & (cur_bump_t <= STOPPING_WINDOW_POST_STIM(2) ));        
%             cur_fwd_vel_in_stopping_win = mean( total_speed( stopping_win ) );
%             
%             total_speeds(end+1) = cur_fwd_vel_in_stopping_win;
%             
%             cur_bump   = BUMP_TO_DEG_CONVERSION_FACTOR*cur_bump_set( s, : );
%             cur_bump_t = cur_bump_t_set( s, : );
%                         
%             figure;
%             hold on;
%             yyaxis right;
%             hold on; 
%             plot( cur_bump_t, cur_bump );
%             yyaxis left;
%             plot( full_yaw_t{1,1}(1,:), cur_fwd );
%             
%             bump_jump_win = find( cur_bump_t >= BUMP_JUMP_WINDOW(1) & (cur_bump_t <= BUMP_JUMP_WINDOW(2) ));
%                     
%             % disp(size(bump_jump_win));
%             cur_mean_bump_jump_mag = abs(mean( cur_bump( bump_jump_win ) ));
%             
%             if( cur_fwd_vel_in_stopping_win <= STOPPING_THRESHOLD )
%                 stopping_event_cnt = stopping_event_cnt + 1;
%                 
%                 % t=0 stim trigger
%                 stopping_win = find( cur_bump_t >= 0 );
%                 stop_durations(end+1) = compute_stop_duration( total_speed(stopping_win), cur_bump_t(2)-cur_bump_t(1), STOPPING_THRESHOLD );
%                 
%                 bump_jump_magnitudes_stopped( end+1 ) = cur_mean_bump_jump_mag;
%             else
%                 bump_jump_magnitudes_moving( end+1 ) = cur_mean_bump_jump_mag;
%             end
%             
%             total_event_cnt = total_event_cnt + 1;
%         end
%     end
% end
% 
% duration_of_stop_avg = mean( stop_durations );
% duration_of_stop_sd = std( stop_durations );
% 
% f = figure;
% hold on;
% text( 0, 0.2, ['Probability of stopping: ' num2str(stopping_event_cnt/total_event_cnt) ' N trials: ' num2str(total_event_cnt)]);
% text( 0, 0.5, ['Avg duration of stopping: ' num2str(duration_of_stop_avg) ' (SD): ' num2str(duration_of_stop_sd) ' N trials: ' num2str(total_event_cnt)]);
% 
% avg_bump_mag_moving = mean( bump_jump_magnitudes_moving );
% std_bump_mag_moving = std( bump_jump_magnitudes_moving );
% 
% avg_bump_mag_stopped = mean( bump_jump_magnitudes_stopped );
% std_bump_mag_stopped = std( bump_jump_magnitudes_stopped );
% 
% text( 0, 0.8, ['Bump mag moving: ' num2str( avg_bump_mag_moving ) ' (SD): ' num2str( std_bump_mag_moving ) ' N trials: ' num2str(length(bump_jump_magnitudes_moving))]);
% text( 0, 1.0, ['Bump mag stopped: ' num2str( avg_bump_mag_stopped ) ' (SD): ' num2str( std_bump_mag_stopped ) ' N trials: ' num2str(length(bump_jump_magnitudes_stopped))]);
% 
% ylim([0 1.4]);
% 
% saveas(f, [savepath 'stopping_analysis.fig']);
% saveas(f, [savepath 'stopping_analysis.png']);
