function display_bump_jump_colorcode_returns( basedir, exp_directories, final_stims_passed_all_checks )



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
% JUMP_THRESHOLD = 0.5; % units of EPG wedge CHANGED ON 4/12/2019
JUMP_THRESHOLD = 0.25; % units of EPG wedge

BUMP_RETURN_MAX_TIME_CUTOFF = 7.0; % seconds

BUMP_RETURN_STABILITY_WINDOW = 0.5; % seconds CHANGED ON 4/12/2019
% BUMP_RETURN_STABILITY_WINDOW = 2.5; % seconds

BUMP_RETURN_ZONE = 45; % EPG wedge

BUMP_RETURN_STATE_OUTSIDE = 22;
BUMP_RETURN_STATE_INSIDE = 23;

DEBUG_SHOW_BUMP_TC = 11;
DEBUG_OFF = 12;
DEBUG_LEVEL = DEBUG_OFF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 = bump time cources
% 2 = stim ids
bump_jumps_up_returns_down = cell( length( exp_directories ), 3 );
bump_jumps_down_returns_up = cell( length( exp_directories ), 3 );
no_bump_jump = cell( length( exp_directories ), 2 );

for d = 1:length( exp_directories )
    
    JUMP_THRESHOLD               = exp_directories{d}{4};
    BUMP_RETURN_STABILITY_WINDOW = exp_directories{d}{5};
    
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
    
    % f = figure('units','normalized','outerposition',[0 0 1 1]);
    smothed_bump_all_stims = [];
    vect_strength_check_all = [];
    
    debug_stop = 1;
    
    all_bump_data = [];
    
    % Compute bump threshold percentile
    for st = 1:size( bump_data, 1 )        
        
        cur_bump_data = squeeze( bump_data( st, :, : ) );
                   
        [ smoothed_bump, cur_bump_tc, cur_bump_tc_unwrapped, vect_strength_check ] = get_radial_weighted_avg_bump_pos_vect_strengh_check_v2( cur_bump_data );
        
        vect_strength_check_all = horzcat( vect_strength_check_all, vect_strength_check );
    end
    
    PERCENTILE_CHECK = 7;    
    per = prctile( vect_strength_check_all, PERCENTILE_CHECK );

    BUMP_THREHOSLD_PERCENTILE = per;
   
    vect_strength_check_all = [];
    for st = 1:size( bump_data, 1 )
        
        cur_bump_data = squeeze( bump_data( st, :, : ) );
           
        all_bump_data = horzcat( all_bump_data, cur_bump_data );
        
        if(st == 48)
            stop_here = stop_here + 1;
        end
        
        [ smoothed_bump, cur_bump_tc, cur_bump_tc_unwrapped, bump_in_deg_unwrapped, vect_strength_check ] = get_radial_weighted_avg_bump_pos_vect_without_periodic_jumps( cur_bump_data, BUMP_THREHOSLD_PERCENTILE );
        
        vect_strength_check_all = horzcat( vect_strength_check_all, vect_strength_check );
        
        smothed_bump_all_stims(st,:,:) = smoothed_bump;
        
        % Rotate bump data so that the average pre-stim period is in the same
        % place for each trial.
        BUMP_TC_FILTER_SAMPLE_POINTS = 5;
        baseline_vals = cur_bump_tc( bump_baseline_idx );
        baseline_non_nan = baseline_vals(~isnan(baseline_vals));        
        bump_delta_tc = medfilt1( cur_bump_tc - mean(baseline_non_nan), BUMP_TC_FILTER_SAMPLE_POINTS, 'truncate' );      
        
        pre_stim_bump_t = find( (t_bump_w < 0) & ( t_bump_w > -1*STABLE_BUMP_BEFORE_STIM_T ));        
        pre_stim_bump_tc = bump_delta_tc( pre_stim_bump_t );

        baseline_vals_uw = cur_bump_tc_unwrapped( bump_baseline_idx );
        baseline_non_nan_uw = baseline_vals_uw(~isnan(baseline_vals_uw));        
        bump_delta_tc_uw = medfilt1( cur_bump_tc_unwrapped - mean(baseline_non_nan_uw), BUMP_TC_FILTER_SAMPLE_POINTS, 'truncate' );        
        pre_stim_bump_tc_uw = bump_delta_tc_uw( pre_stim_bump_t );
        
        baseline_vals_deg_uw = bump_in_deg_unwrapped( bump_baseline_idx );
        baseline_non_nan_deg_uw = baseline_vals_deg_uw(~isnan(baseline_vals_deg_uw));        
        bump_delta_tc_deg_uw = medfilt1( bump_in_deg_unwrapped - mean(baseline_non_nan_deg_uw), 4, 'truncate' );        
        % bump_delta_tc_deg_uw = bump_in_deg_unwrapped - mean(baseline_non_nan_deg_uw);
                
        stop_here = 0;
        
        % Check that none of the values in the prestim period are NaNs, and that the bump was stable        
        if( ( sum( isnan( pre_stim_bump_tc_uw ) ) == 0 ) && ( std( pre_stim_bump_tc_uw ) < PRE_STIM_BUMP_STABILITY_THRESHOLD  ) )
            passed_bump_stability_check(end+1) = st;
            
            bump_delta_tc_deg_uw(isnan(cur_bump_tc)) = NaN;           
            bump_delta_tc_deg_uw_fixed = bump_delta_tc_deg_uw;
            bump_delta_tc_deg_uw_fixed = assess_and_fix_bump_quality( bump_delta_tc_deg_uw, VPS );
         
            if(length(bump_delta_tc_deg_uw_fixed) == 0 )
                disp(['Trial: ' num2str(st) ' has been disqualified, too many time points with no clear bump.']);
                continue;
            end
            
if 0            
            WRAP_ANGLE = 360-45;
            
            for jj = 1:length(bump_delta_tc_deg_uw_fixed)
                
                cur_bump_deg = bump_delta_tc_deg_uw_fixed(jj);
                
                if( cur_bump_deg >= WRAP_ANGLE )
                    bump_delta_tc_deg_uw_fixed(jj) = bump_delta_tc_deg_uw_fixed(jj) - 360;
                elseif( cur_bump_deg <= -1.0*WRAP_ANGLE )
                    bump_delta_tc_deg_uw_fixed(jj) = bump_delta_tc_deg_uw_fixed(jj) + 360;                    
                end
            end
end

            passed_bump_stability_tc(end+1,:) =  bump_delta_tc_deg_uw_fixed;
            
            % subplot(1,2,1);
            % hold on;
            % plot( t_bump_w, bump_delta_tc, 'DisplayName', ['stim: ' num2str(st)] );  
            stop_here = stop_here + 1;
        else
            failed_bump_stability_check(end+1) = st;
            % subplot(1,2,2);
            % hold on;
            % plot( t_bump_w, bump_delta_tc, 'DisplayName', ['stim: ' num2str(st)] );
            stop_here = stop_here + 1;
        end
    end
    
%     DEBUG_BUMP_STATISTICS = 0;
%     if( DEBUG_BUMP_STATISTICS == 1 )
%         analyze_bump_statistics( all_bump_data );
%     end
%     
%     % Check distribution of bump values
%     DEBUG_BUMP_QUALITY = 0;
%     if( DEBUG_BUMP_QUALITY )
%         figure;
%         
% %         histogram( vect_strength_check_all, 50 );
% %         disp(['Avg bump strength vector: ' num2str( mean(vect_strength_check_all)) ]);
% %         disp(['Std bump strength vector: ' num2str( std(vect_strength_check_all)) ]);
% %         
% %         parmhat = evfit( vect_strength_check_all );
% %         parmhat
% 
%         OUTLIER_LIMIT = 5.0;
%         vect_strength_check_all( find(vect_strength_check_all > OUTLIER_LIMIT ) ) = OUTLIER_LIMIT;
%         avg_vect = mean( vect_strength_check_all );
%         std_vect = std( vect_strength_check_all );
%         
%         BUMP_QUALITY_THRESHOLD = 0.2;
%         ip = invprctile(vect_strength_check_all,BUMP_QUALITY_THRESHOLD);
%         
%         PERCENTILE_CHECK = 7;
%         per = prctile(vect_strength_check_all,PERCENTILE_CHECK);
%                 
%         disp(['Avg: ' num2str(avg_vect) ' Std: ' num2str(std_vect) ' Percentile (' num2str(PERCENTILE_CHECK) ') check: ' num2str(per) ]);
%         
%         figure(f);
%     end
    
%     
%     subplot(1,2,1)
%     title(['Passed bump stability check: ' num2str(length(passed_bump_stability_check))]);
%     xlabel('Time (s)');
%     ylabel('EB bump position');
%     
%     subplot(1,2,2)
%     title(['Failed bump stability check: ' num2str(length(failed_bump_stability_check))]);
%     xlabel('Time (s)');
%     ylabel('EB bump position');
%     
%     saveas(f, [analysis_path '/' cur_datapath '_post_stim_bump_stability_check.fig']);
%     saveas(f, [analysis_path '/' cur_datapath '_post_stim_bump_stability_check.png']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Classify trials that passed bump stability check, by trials these
    % categories;
    % 3. Bump returned up
    % 4. Bump returned down
    % 5. Bump jumped not return
    % 7. No bump jump 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    BUMP_RETURN_STABILITY_WINDOW_IN_SAMPLES = floor(BUMP_RETURN_STABILITY_WINDOW * VPS);
    
    jump_lookahead_window = find( (t_bump_w < JUMP_LOOKAHEAD_WINDOW) & ( t_bump_w > 0 ));
    
    % f = figure('units','normalized','outerposition',[0 0 1 1]);
    f = figure();
    debug_cnt = 1;
    total_cnt = 0;
    returned_cnt = 0;
    return_type_1_idx = 0;
    return_type_2_idx = 0;
    
    for st = 1:size( passed_bump_stability_tc,1 )
        
        if( st == 49 )
            continue;
        end
        cur_bump_tc = passed_bump_stability_tc( st, : );          
        
        cur_mean_bump_jump = mean( cur_bump_tc( jump_lookahead_window ) );
                    
        % Criteria for bump return: bump crosses the initial position of
        % 0 and stays there for at least some time
        if( abs(cur_mean_bump_jump) > JUMP_THRESHOLD )
            % Jump up or down
            
%             subplot(1,2,1);
%             hold on;
%             plot( t_bump_w, cur_bump_tc, 'DisplayName', ['stim: ' num2str(st)], 'color', 'b' );  
            
            % Evaluate if the bump returns or not, with maximum amount of
            % time
            bump_return_cutoff = find( t_bump_w <  BUMP_RETURN_MAX_TIME_CUTOFF );
            bump_return_eval_period = [jump_lookahead_window(end) : bump_return_cutoff(end)];
            
            bump_returns_flag = 0;
            returned_time_idx = -1;
            
            % Determine the first crossing back to the home zone.
            bump_home_return_crossings_0 = find( abs(cur_bump_tc(bump_return_eval_period)) < BUMP_RETURN_ZONE );
            
            % 
            bump_home_return_crossings_360 = find(( abs(cur_bump_tc(bump_return_eval_period)) < 360+BUMP_RETURN_ZONE ) & ( abs(cur_bump_tc(bump_return_eval_period)) > 360-BUMP_RETURN_ZONE) );
            
            % Check 0 
            if( length( bump_home_return_crossings_0 ) == 0 ) 
                bump_returns_flag = 0;
                returned_time_idx = -1;
            else
                % Check that the average bump position is within the home zone
                % for a certain window after the initial crossing.
                first_crossing_index = bump_home_return_crossings_0(1);
                
                if( (first_crossing_index+BUMP_RETURN_STABILITY_WINDOW_IN_SAMPLES) <= length(bump_return_eval_period) )
                    end_of_crossing_window = first_crossing_index + BUMP_RETURN_STABILITY_WINDOW_IN_SAMPLES;
                else
                    end_of_crossing_window = length( bump_return_eval_period );
                end
                
                avg_bump_in_home_window = mean( abs(cur_bump_tc( bump_return_eval_period( first_crossing_index : end_of_crossing_window ))));
                
                if( avg_bump_in_home_window < BUMP_RETURN_ZONE )
                    bump_returns_flag = 1;
                    returned_time_idx = bump_return_eval_period( end_of_crossing_window );
                    
                    % Classify the direction of bump return, as the
                    % vector between the bump jump peak and the average
                    % bump return path
                    post_jump_non_nan = cur_bump_tc(bump_return_eval_period(~isnan(cur_bump_tc(bump_return_eval_period(1:end_of_crossing_window))))); 
                    bump_return_direction = mean(post_jump_non_nan) - cur_bump_tc(bump_return_eval_period(1));
                    
                else
                    % Failed to find a home crossing.
                    bump_returns_flag = 0;
                    returned_time_idx = -1;
                end
            end

            % Check 360 
            if(bump_returns_flag == 0)
                if( length( bump_home_return_crossings_360 ) == 0 )
                    bump_returns_flag = 0;
                    returned_time_idx = -1;
                else
                    % Check that the average bump position is within the home zone
                    % for a certain window after the initial crossing.
                    first_crossing_index = bump_home_return_crossings_360(1);
                    
                    if( (first_crossing_index+BUMP_RETURN_STABILITY_WINDOW_IN_SAMPLES) <= length(bump_return_eval_period) )
                        end_of_crossing_window = first_crossing_index + BUMP_RETURN_STABILITY_WINDOW_IN_SAMPLES;
                    else
                        end_of_crossing_window = length( bump_return_eval_period );
                    end
                    
                    avg_bump_in_home_window = mean( abs(cur_bump_tc( bump_return_eval_period( first_crossing_index : end_of_crossing_window ))));
                    
                    if( ( avg_bump_in_home_window < (360 + BUMP_RETURN_ZONE) ) & ( avg_bump_in_home_window > (360 - BUMP_RETURN_ZONE) ) )
                        bump_returns_flag = 1;
                        returned_time_idx = bump_return_eval_period( end_of_crossing_window );
                        
                        % Classify the direction of bump return, as the
                        % vector between the bump jump peak and the average
                        % bump return path
                        post_jump_non_nan = cur_bump_tc(bump_return_eval_period(~isnan(cur_bump_tc(bump_return_eval_period(1:end_of_crossing_window)))));
                        bump_return_direction = mean(post_jump_non_nan) - cur_bump_tc(bump_return_eval_period(1));
                        
                    else
                        % Failed to find a home crossing.
                        bump_returns_flag = 0;
                        returned_time_idx = -1;
                    end
                end
            end
            
            return_type = -1;            
            subplot(1,1,1);
            hold on;
            if( bump_returns_flag ==  0)
                plot( t_bump_w, cur_bump_tc, 'DisplayName', [ 'stim: ' num2str( st ) ], 'color', rgb('Silver') );  
            else
                
                if( bump_return_direction > 0 )
                    return_type = 1;
                    return_type_1_idx = return_type_1_idx + 1;                    
                    cur_return_type_idx = return_type_1_idx;
                    disp('Got to return type 1');
                elseif( bump_return_direction < 0 )
                    return_type = 2;
                    return_type_2_idx = return_type_2_idx + 1;                                        
                    cur_return_type_idx = return_type_2_idx;
                    disp('Got to return type 2');
                end
                
                if( ismember( cur_return_type_idx, final_stims_passed_all_checks{ return_type } ) )                
                    plot( t_bump_w, cur_bump_tc, 'DisplayName', [ 'stim: ' num2str( st ) ], 'color', 'r' );    
                    returned_cnt = returned_cnt + 1;
                end
            end
            total_cnt = total_cnt + 1;
        end
    end
        
    xlim([-1 10]);
    % title(['Percent returned: ' num2str( 100.0*returned_cnt/total_cnt ) ' ( ' num2str( returned_cnt ) '/' num2str( total_cnt ) ' )']);
    title(['Returned count: ' num2str( returned_cnt ) ]);
    xlabel('Time from ATP stim (s)');
    xlabel('Bump position (deg)');
    set( gca(), 'TickDir', 'Out' );
    saveas(f, [analysis_path '/' cur_datapath '_bump_jump_colorcoded_return.fig']);
    saveas(f, [analysis_path '/' cur_datapath '_bump_jump_colorcoded_return.png']);
    saveas(f, [analysis_path '/' cur_datapath '_bump_jump_colorcoded_return.svg']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end
end