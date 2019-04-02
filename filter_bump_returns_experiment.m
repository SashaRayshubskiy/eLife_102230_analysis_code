function [ bump_jumps_up_returns_down, bump_jumps_down_returns_up, no_bump_jump ] = filter_bump_returns_experiment( basedir, exp_directories )

settings = sensor_settings;
BALL_FR = settings.sensorPollFreq;
EPHYS_FR = settings.sampRate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BEFORE_STIM_TIME = 2.0; % seconds
AFTER_STIM_TIME  = 10.0; % seconds

STABLE_BUMP_BEFORE_STIM_T = 1.0; % seconds
PRE_STIM_BUMP_STABILITY_THRESHOLD = 1.5; % units of EPG wedges

JUMP_LOOKAHEAD_WINDOW = 1.5; % seconds
JUMP_THRESHOLD = 0.5; % units of EPG wedge

BUMP_RETURN_MAX_TIME_CUTOFF = 7.0; % seconds
BUMP_RETURN_STABILITY_WINDOW = 1.0; % seconds

BUMP_RETURN_ZONE = 1.0; % EPG wedge

BUMP_RETURN_STATE_OUTSIDE = 22;
BUMP_RETURN_STATE_INSIDE = 23;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 = bump time cources
% 2 = stim ids
bump_jumps_up_returns_down = cell( length( exp_directories ), 3 );
bump_jumps_down_returns_up = cell( length( exp_directories ), 3 );
no_bump_jump = cell( length( exp_directories ), 2 );
    
for d = 1:length( exp_directories )
    
    cur_datapath = exp_directories{d}{1};
    cur_sid      = exp_directories{d}{2};
    
    datapath = [ basedir '/' cur_datapath ]; 
    analysis_path = [datapath '/analysis/'];
    stim_filepath = [analysis_path '/stim_window_data_sid_' num2str(cur_sid) '.mat'];
    
    cur_data = load( stim_filepath );
    % 'fwd_in_window', 'yaw_in_window', 'ephys_in_window', 'bump_in_window'
    
    raw_data_filepath = [analysis_path '/raw_glom_ball_ephys.mat'];
    VPS_tmp = load( raw_data_filepath );
    VPS = VPS_tmp.VPS;
    clear VPS_tmp;    
    
    bump_data   = cur_data.bump_in_window;
    yaw_data    = cur_data.yaw_in_window;
    ephys_data  = cur_data.ephys_in_window;
    
    bump_baseline_idx = [ floor(VPS * (BEFORE_STIM_TIME-0.5)) : floor(VPS * BEFORE_STIM_TIME) ];
    
    t_bump_w   = ([0:size(bump_data,3)-1] / VPS) - BEFORE_STIM_TIME;
    t_yaw_w    = ([0:size(yaw_data,2)-1] / BALL_FR) - BEFORE_STIM_TIME;
    t_ephys_w  = ([0:size(ephys_data,2)-1] / EPHYS_FR) - BEFORE_STIM_TIME;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filter trials with a stable bump for a period of time before stim.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    passed_bump_stability_check = [];
    failed_bump_stability_check = [];
    passed_bump_stability_tc = [];
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    smothed_bump_all_stims = [];
    vect_strength_check_all = [];
    
    for st = 1:size( bump_data, 1 )
        
        % [ smoothed_bump, cur_bump_tc ] = get_radial_weighted_avg_bump_pos_v3( squeeze( bump_data( st, :, : ) ) );
        [ smoothed_bump, cur_bump_tc, vect_strength_check ] = get_radial_weighted_avg_bump_pos_vect_strengh_check( squeeze( bump_data( st, :, : ) ) );
        
        vect_strength_check_all = horzcat( vect_strength_check_all, vect_strength_check );
        
        smothed_bump_all_stims(st,:,:) = smoothed_bump;
        
        % Rotate bump data so that the average pre-stim period is in the same
        % place for each trial.
        baseline_vals = cur_bump_tc( bump_baseline_idx );
        baseline_non_nan = baseline_vals(~isnan(baseline_vals));
        
        bump_delta_tc = medfilt1( cur_bump_tc - mean(baseline_non_nan), 5, 'truncate' );
        
        pre_stim_bump_t = find( (t_bump_w < 0) & ( t_bump_w > -1*STABLE_BUMP_BEFORE_STIM_T ));
        
        pre_stim_bump_tc = bump_delta_tc( pre_stim_bump_t );
        
        % Check that none of the values in the prestim period are NaNs, and that the bump was stable        
        if( ( sum( isnan( pre_stim_bump_tc ) ) == 0 ) && ( std( pre_stim_bump_tc ) < PRE_STIM_BUMP_STABILITY_THRESHOLD ) )
            passed_bump_stability_check(end+1) = st;
            passed_bump_stability_tc(end+1,:) = bump_delta_tc;
            
            subplot(1,2,1);
            hold on;
            plot( t_bump_w, bump_delta_tc );            
        else
            failed_bump_stability_check(end+1) = st;
            subplot(1,2,2);
            hold on;
            plot( t_bump_w, bump_delta_tc );
        end
    end
    
    DEBUG_BUMP_QUALITY = 0;
    if( DEBUG_BUMP_QUALITY )
        figure;
        
        histogram( vect_strength_check_all, 1000 );
        
        figure(f);
    end
    
    subplot(1,2,1)
    title(['Passed bump stability check: ' num2str(length(passed_bump_stability_check))]);
    xlabel('Time (s)');
    ylabel('EB bump position');
    
    subplot(1,2,2)
    title(['Failed bump stability check: ' num2str(length(failed_bump_stability_check))]);
    xlabel('Time (s)');
    ylabel('EB bump position');
    
    saveas(f, [analysis_path '/' cur_datapath '_post_stim_bump_stability_check.fig']);
    saveas(f, [analysis_path '/' cur_datapath '_post_stim_bump_stability_check.png']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Classify trials that passed bump stability check, by trials these
    % categories;
    % 1. Bump jumped up and returned
    % 2. Bump jumped down and returned
    % 3. Bump jumped up and did not return
    % 4. Bump jumped down and did not return
    % 5. No bump jump 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    jump_up_and_returned_tc = [];       jump_up_and_returned_ids = [];
    jump_down_and_returned_tc = [];     jump_down_and_returned_ids = [];
    jump_up_and_not_returned_tc = [];   jump_up_and_not_returned_ids = [];
    jump_down_and_not_returned_tc = []; jump_down_and_not_returned_ids = [];
    no_response_tc = [];                no_response_ids = [];

    jump_up_and_returned_time_idx = {};
    jump_down_and_returned_time_idx = {};
    
    BUMP_RETURN_STABILITY_WINDOW_IN_SAMPLES = floor(BUMP_RETURN_STABILITY_WINDOW * VPS);
    
    jump_lookahead_window = find( (t_bump_w < JUMP_LOOKAHEAD_WINDOW) & ( t_bump_w > 0 ));
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    for st = 1:length( passed_bump_stability_check )
        cur_stim = passed_bump_stability_check( st );
        cur_bump_tc = passed_bump_stability_tc( st, : );          
        
        cur_mean_bump_jump = mean( cur_bump_tc( jump_lookahead_window ) );
        
        % Criteria for bump return: bump crosses the initial position of
        % 0 and stays there for at least some time
        if( abs(cur_mean_bump_jump) > JUMP_THRESHOLD )
            % Jump up or down
            
            % Evaluate if the bump returns or not, with maximum amount of
            % time
            bump_return_cutoff = find( t_bump_w <  BUMP_RETURN_MAX_TIME_CUTOFF );
            bump_return_eval_period = [jump_lookahead_window(end) : bump_return_cutoff(end)];

            CUR_RETURN_STATE = BUMP_RETURN_STATE_OUTSIDE;
            bump_returns_flag = 0;

            returned_time_idx = -1;
            kk = 1;
            while( kk < length( bump_return_eval_period ) )
            % for j = bump_return_eval_period
                cur_j = bump_return_eval_period(kk);
                cur_bump_pos = cur_bump_tc( cur_j );
                
                if( CUR_RETURN_STATE == BUMP_RETURN_STATE_OUTSIDE )
                    
                    if( abs(cur_bump_pos) <= BUMP_RETURN_ZONE )
                        CUR_RETURN_STATE = BUMP_RETURN_STATE_INSIDE;
                        inside_return_zone_cnt = 1;
                    end
                
                elseif ( CUR_RETURN_STATE == BUMP_RETURN_STATE_INSIDE )
                    if( abs(cur_bump_pos) <= BUMP_RETURN_ZONE )
                        inside_return_zone_cnt = inside_return_zone_cnt + 1;
                    end
                    
                    if( inside_return_zone_cnt >= BUMP_RETURN_STABILITY_WINDOW_IN_SAMPLES ) 
                        % Victory, we found a bump return.
                        bump_returns_flag = 1;
                        returned_time_idx = cur_j;
                        
                        % Classify the direction of bump return, as the
                        % vector between the bump jump peak and the average
                        % bump return path
                        bump_return_direction = mean(cur_bump_tc(bump_return_eval_period(1:kk))) - cur_bump_tc(bump_return_eval_period(1));
                        break;
                    end
                end
                kk = kk + 1;
            end
            
            if( ( cur_mean_bump_jump > JUMP_THRESHOLD ) && (bump_returns_flag == 1))
                subplot(5,1,1);
                hold on;
                
                jump_up_and_returned_tc(end+1,:) = cur_bump_tc;       
                jump_up_and_returned_ids(end+1) = cur_stim; 
                jump_up_and_returned_time_idx{end+1} = [ jump_lookahead_window(end) : returned_time_idx ];
                
                
                
                
            elseif( ( cur_mean_bump_jump > JUMP_THRESHOLD ) && (bump_returns_flag == 0))
                subplot(5,1,2);
                hold on;
                jump_up_and_not_returned_tc(end+1,:) = cur_bump_tc;   jump_up_and_not_returned_ids(end+1) = cur_stim;
            
            elseif( ( cur_mean_bump_jump < -1.0*JUMP_THRESHOLD ) && (bump_returns_flag == 1))
                subplot(5,1,3);
                hold on;
                
                jump_down_and_returned_tc(end+1,:) = cur_bump_tc;     
                jump_down_and_returned_ids(end+1) = cur_stim; 
                jump_down_and_returned_time_idx{end+1} = [ jump_lookahead_window(end) : returned_time_idx ];
            
            elseif( ( cur_mean_bump_jump < -1.0*JUMP_THRESHOLD ) && (bump_returns_flag == 0))
                subplot(5,1,4);
                hold on;
                jump_down_and_not_returned_tc(end+1,:) = cur_bump_tc; jump_down_and_not_returned_ids(end+1) = cur_stim;            
            end
        else
            % No response
            no_response_tc(end+1,:) = cur_bump_tc;                    no_response_ids(end+1) = cur_stim;       
            subplot(5,1,5);
            hold on;            
        end
        
        plot( t_bump_w, cur_bump_tc );
    end
    
    subplot(5,1,1);
    title(['Jump up and returned down: ' num2str(length(jump_up_and_returned_ids))]);
    subplot(5,1,2);
    title(['Jump up and not returned: ' num2str(length(jump_up_and_not_returned_ids))]);
    subplot(5,1,3);
    title(['Jump down and returned up: ' num2str(length(jump_down_and_returned_ids))]);
    subplot(5,1,4);
    title(['Jump down and not returned: ' num2str(length(jump_down_and_not_returned_ids))]);
    subplot(5,1,5);
    title(['No response: ' num2str(length(no_response_ids))]);
    
    saveas(f, [analysis_path '/' cur_datapath '_bump_jump_classification.fig']);
    saveas(f, [analysis_path '/' cur_datapath '_bump_jump_classification.png']);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Return stims with succesful bump jump up/down with returns
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    bump_jumps_up_returns_down{d,1}   = jump_up_and_returned_tc;
    bump_jumps_up_returns_down{d,2}   = jump_up_and_returned_ids;
    bump_jumps_up_returns_down{d,3}   = jump_up_and_returned_time_idx;
    
    bump_jumps_down_returns_up{d,1} = jump_down_and_returned_tc;
    bump_jumps_down_returns_up{d,2} = jump_down_and_returned_ids;
    bump_jumps_down_returns_up{d,3} = jump_down_and_returned_time_idx;
    
    no_bump_jump{d,1} = no_response_tc;
    no_bump_jump{d,2} = no_response_ids;
end
end