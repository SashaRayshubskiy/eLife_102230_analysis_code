function [ bump_jump_percent_per_fly ] = calculate_percent_of_bump_jumps( basedir, cur_dirs, experiment_type_str )

settings = sensor_settings;
BALL_FR = settings.sensorPollFreq;
EPHYS_FR = settings.sampRate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BEFORE_STIM_TIME = 2.0; % seconds

JUMP_LOOKAHEAD_WINDOW = 2.0; % seconds

STABLE_BUMP_BEFORE_STIM_T = 0.5; % seconds

PRE_STIM_BUMP_STABILITY_THRESHOLD = 1.5; % units of EPG wedges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d = 1:length( cur_dirs )
% for d = 2
    
    JUMP_THRESHOLD = cur_dirs{d}{4};
    
    cur_datapath   = cur_dirs{d}{1};
    cur_sid        = cur_dirs{d}{2};
    
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
    
    t_bump_w   = ([0:size(bump_data,3)-1] / VPS) - BEFORE_STIM_TIME;

    pre_stim_bump_t = find( (t_bump_w < 0) & ( t_bump_w > -1*STABLE_BUMP_BEFORE_STIM_T ));        
    
    jump_lookahead_window = find( (t_bump_w < JUMP_LOOKAHEAD_WINDOW) & ( t_bump_w > 0 ));
    
    passed_bump_stability_check = [];
    failed_bump_stability_check = [];
    passed_bump_stability_tc = [];
    
    smothed_bump_all_stims = [];
    vect_strength_check_all = [];
    
    f = figure();
    
    debug_stop = 1;
    for st = 1:size( bump_data, 1 )
        
        % [ smoothed_bump, cur_bump_tc ] = get_radial_weighted_avg_bump_pos_v3( squeeze( bump_data( st, :, : ) ) );
        % [ smoothed_bump, cur_bump_tc, vect_strength_check ] = get_radial_weighted_avg_bump_pos_vect_strengh_check( squeeze( bump_data( st, :, : ) ) );
        
%         if( st == 29 )
%             debug_stop = debug_stop + 1;
%             bump_debug = squeeze( bump_data( st, :, : ) ) ;
%             save('/tmp/bump_debug.mat', 'bump_debug', 'bump_baseline_idx');
%         end
        
        [ smoothed_bump, cur_bump_tc, cur_bump_tc_unwrapped, vect_strength_check ] = get_radial_weighted_avg_bump_pos_vect_strengh_check_v2( squeeze( bump_data( st, :, : ) ) );
        
        vect_strength_check_all = horzcat( vect_strength_check_all, vect_strength_check );
        
        smothed_bump_all_stims(st,:,:) = smoothed_bump;
        
        % Rotate bump data so that the average pre-stim period is in the same
        % place for each trial.
        BUMP_TC_FILTER_SAMPLE_POINTS = 10;
        baseline_vals = cur_bump_tc( pre_stim_bump_t );
        baseline_non_nan = baseline_vals(~isnan(baseline_vals));        
        bump_delta_tc = medfilt1( cur_bump_tc - mean(baseline_non_nan), BUMP_TC_FILTER_SAMPLE_POINTS, 'truncate' );             

        baseline_vals_uw = cur_bump_tc_unwrapped( pre_stim_bump_t );
        baseline_non_nan_uw = baseline_vals_uw(~isnan(baseline_vals_uw));        
        bump_delta_tc_uw = medfilt1( cur_bump_tc_unwrapped - mean(baseline_non_nan_uw), BUMP_TC_FILTER_SAMPLE_POINTS, 'truncate' );        
        pre_stim_bump_tc_uw = bump_delta_tc_uw( pre_stim_bump_t );
        
        stop_here = 0;
        
        % Check that none of the values in the prestim period are NaNs, and that the bump was stable        
        if( ( sum( isnan( pre_stim_bump_tc_uw ) ) == 0 ) && ( std( pre_stim_bump_tc_uw ) <= PRE_STIM_BUMP_STABILITY_THRESHOLD  ) )
            passed_bump_stability_check(end+1) = st;
            % passed_bump_stability_tc(end+1,:) = bump_delta_tc;
            passed_bump_stability_tc(end+1,:) = bump_delta_tc_uw;
            
            subplot(1,2,1);
            hold on;
            plot( t_bump_w, bump_delta_tc, 'DisplayName', ['stim: ' num2str(st)] );  
            % stop_here = stop_here + 1;
        else
            failed_bump_stability_check(end+1) = st;
            subplot(1,2,2);
            hold on;
            plot( t_bump_w, bump_delta_tc, 'DisplayName', ['stim: ' num2str(st)] );
            % stop_here = stop_here + 1;
        end
    end
    
    total_passed_stims = length( passed_bump_stability_check );
    jump_cnt = 0;

    f = figure;
    for st = 1:total_passed_stims
        cur_stim = passed_bump_stability_check( st );
        cur_bump_tc = passed_bump_stability_tc( st, : );             
              
        bump_in_win = cur_bump_tc(jump_lookahead_window);
        bump_in_win_non_nan = bump_in_win( ~isnan( bump_in_win ) );
        
        cur_mean_bump_jump = mean( bump_in_win_non_nan );
                
        if( abs(cur_mean_bump_jump) > JUMP_THRESHOLD )
            jump_cnt = jump_cnt + 1;
            subplot( 2, 1, 1 );
            hold on;
            
            % plot(cur_bump_tc( jump_lookahead_window ));
            plot( t_bump_w, cur_bump_tc );
            title(['Jump occured']);
        else
            subplot( 2, 1, 2 );
            hold on;
            
            % plot(cur_bump_tc( jump_lookahead_window ));            
            plot( t_bump_w, cur_bump_tc );
            title('Jump did not occured');
        end
    end
    
    cur_percentage_jump = jump_cnt / total_passed_stims;
    bump_jump_percent_per_fly( d ) = cur_percentage_jump;    
    title(['Fly: ' num2str(d) ' total passed stims: ' num2str( total_passed_stims ) ' jump cnt: ' num2str(jump_cnt)]);
end
end

