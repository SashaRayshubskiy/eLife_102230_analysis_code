function extract_stim_windows( basedir, directories )

DISPLAY_INDIVIDUAL_STIMS = 1;
DISPLAY_ALL_STIMS        = 1;

ac = get_analysis_constants;

BEFORE_STIM_TIME = 2.0; % seconds
AFTER_STIM_TIME  = 10.0; % seconds

STIM_TIME_WIN = 2.0; % KEEP THIS CONSTANT, OTHERWISE STIM IDs change

settings = sensor_settings;
BALL_FR = settings.sensorPollFreq;
EPHYS_FR = settings.sampRate;

EPHYS_BEFORE_STIM_FRAMES = EPHYS_FR*BEFORE_STIM_TIME;
EPHYS_AFTER_STIM_FRAMES = EPHYS_FR*AFTER_STIM_TIME;

for d = 1:length( directories )
    
    cur_datadir = directories{ d }{ 1 };
    cur_sid     = directories{ d }{ 2 };
    
    analysis_path         = [ basedir '/' cur_datadir '/analysis/' ];
    path_to_key_variables = [ analysis_path '/raw_glom_ball_ephys.mat' ];
    cur_data              = load( path_to_key_variables );
    
    cur_F                 = cur_data.glom_EPG_F_per_trial;
    ephys_time            = cur_data.ephys_time;
    ephys_data            = cur_data.ephys_data{1};
    
    ball_time             = cur_data.ball_time;
    ball_data             = cur_data.ball_data{1};
    
    stim_data             = cur_data.stim_data{1};
    
    VPS                   = cur_data.VPS;
    
    df_f_in_roi_per_trial = get_PB_dF_F_per_trial( cur_F );
    
    fwd_in_window = [];
    yaw_in_window = [];
    ephys_in_window = [];
    bump_in_window = [];
    
    global_stim_idx = 1;
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save all stims
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for tr = 1:size( stim_data, 1 )
        cur_stim = 10.0 * squeeze( stim_data( tr, : ) );
        
        % for each stim get a window of before and after stim
        stim_idx = find(diff(cur_stim) > 3.0 );
        
        % Throw out stims that are close together, was a result of a clogged
        % pipette
        stim_times = ephys_time(stim_idx);
        stim_idx_keep = [];
        prev_st = -999999;
        for st = 1:length( stim_times )-1
            cur_st = stim_times( st );
            next_st = stim_times( st+1 );
            
            if( ( (cur_st - prev_st) > STIM_TIME_WIN ) && ( (next_st - cur_st) > STIM_TIME_WIN ) )
                stim_idx_keep( end+1 ) = stim_idx( st );
            end
            
            prev_st = cur_st;
        end
        
        index_in_map = 1;
        for st = 1:length( stim_idx_keep )
            cur_stim_idx = stim_idx_keep( st );
            
            cur_stim_t = ephys_time( cur_stim_idx );
            
            before_t = cur_stim_t - BEFORE_STIM_TIME;
            after_t  = cur_stim_t + AFTER_STIM_TIME;
            
            if((before_t < 0) || (after_t > ephys_time(end) ) )
                continue;
            end
            
            begin_ephys_frame = cur_stim_idx - EPHYS_BEFORE_STIM_FRAMES;
            end_ephys_frame   = cur_stim_idx + EPHYS_AFTER_STIM_FRAMES;
            
            before_ball_frame = before_t * BALL_FR;
            after_ball_frame  = after_t  * BALL_FR;
            
            begin_bump_frame = before_t * VPS;
            end_bump_frame   = after_t  * VPS;
            
            fwd_in_window(end+1,:)    = squeeze( ball_data( tr, ac.VEL_FWD, before_ball_frame:after_ball_frame ) );
            yaw_in_window(end+1,:)    = squeeze( ball_data( tr, ac.VEL_YAW, before_ball_frame:after_ball_frame ) );
            ephys_in_window(end+1,:)  = squeeze( ephys_data(tr, begin_ephys_frame:end_ephys_frame) );
            bump_in_window(end+1,:,:) = squeeze( df_f_in_roi_per_trial( tr, :, begin_bump_frame:end_bump_frame ) );
            
            index_in_map = index_in_map + 1;
            global_stim_idx = global_stim_idx + 1;
        end
    end
    
    t_bump_w  = ([0:size(bump_in_window,3)-1] / VPS) - BEFORE_STIM_TIME;
    t_yaw_w  = ([0:size(yaw_in_window,2)-1] / BALL_FR) - BEFORE_STIM_TIME;
    t_ephys_w  = ([0:size(ephys_in_window,2)-1] / EPHYS_FR) - BEFORE_STIM_TIME;

    save([analysis_path '/stim_window_data_sid_' num2str( cur_sid ) '.mat'], ... 
        'fwd_in_window', 'yaw_in_window', 'ephys_in_window', 'bump_in_window', 'VPS', 't_bump_w', 't_yaw_w', 't_ephys_w' );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display all stim windows and take an average
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if( DISPLAY_ALL_STIMS == 1 )
        % Display averages in a window
        f = figure;
        bump_tc_all = [];
               
        bump_baseline_idx = [ 1 : VPS * BEFORE_STIM_TIME ];
        
        stim_ids_to_include = [1:size( yaw_in_window, 1 )];
        
        left_bump = [];
        left_yaw = [];
        left_bump = [];
        
        for i = stim_ids_to_include
            
            subplot(3,1,1);
            hold on;
            
            % [dummy, cur_bump_tc] = max( squeeze(bump_in_window(i,:,:)));
            cur_bump_tc = get_radial_weighted_avg_bump_pos_v3( squeeze(bump_in_window(i,:,:)) );
            
            % Rotate bump data so that the average pre-stim period is in the same
            % place for each trial.
            baseline_vals = cur_bump_tc(bump_baseline_idx);
            baseline_non_nan = baseline_vals(~isnan(baseline_vals));
            
            bump_delta_tc = medfilt1( cur_bump_tc - mean(baseline_non_nan), 5, 'truncate' );
            
            bump_tc_all(end+1,:) = bump_delta_tc;
            plot(t_bump_w, bump_delta_tc, 'LineWidth', 1 );
            
            subplot(3,1,2);
            hold on;
            plot(t_yaw_w, squeeze(yaw_in_window(i,:)), 'LineWidth', 1 );
            
            subplot(3,1,3);
            hold on;
            plot(t_ephys_w, squeeze(ephys_in_window(i,:)), 'LineWidth', 1 );
        end
        
        ax1(1) = subplot(3,1,1);
        hold on;
        plot( t_bump_w, squeeze(mean(bump_tc_all)), 'LineWidth', 2, 'color', 'b' );
        ylabel('EB bump location');
        
        ax1(2) = subplot(3,1,2);
        hold on;
        plot( t_yaw_w, squeeze(mean(yaw_in_window(stim_ids_to_include, :))), 'LineWidth', 2, 'color', 'b' );
        ylabel('Yaw (au)');
        
        ax1(3) = subplot(3,1,3);
        hold on;
        plot( t_ephys_w, squeeze(mean(ephys_in_window(stim_ids_to_include,:))), 'LineWidth', 2, 'color', 'b' );
        ylabel('Vm (mV)');
        
        linkaxes(ax1, 'x');
        
        xlabel('Time (s)');
        
        title(['Number of stim: ' num2str(size(yaw_in_window(stim_ids_to_include,:),1))]);
        
        saveas( f, [analysis_path '/avg_EB_yaw_ephys_triggered_on_stim_' num2str( cur_sid ) '.fig'] );
        saveas( f, [analysis_path '/avg_EB_yaw_ephys_triggered_on_stim_' num2str( cur_sid ) '.png'] );
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display individual stims in a separate file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    if( DISPLAY_INDIVIDUAL_STIMS == 1 )
        % Display all individual stim events
        for i = 1:size( yaw_in_window, 1 )
            f = figure;
            
            ax2(1) = subplot(4,1,1);
            hold on;
            imagesc( t_bump_w, [1:size(bump_in_window,2)], squeeze(bump_in_window(i,:,:)) );
            colormap(flipud(gray));
            axis tight;
            caxis([-0.5 4]);
            
            ax2(2) = subplot(4,1,2);
            hold on;
            plot(t_bump_w, squeeze(bump_tc_all(i,:)), 'LineWidth', 1 );
            
            ax2(3) = subplot(4,1,3);
            hold on;
            plot(t_yaw_w, squeeze(yaw_in_window(i,:)), 'LineWidth', 1 );
            
            ax2(4) = subplot(4,1,4);
            hold on;
            plot(t_ephys_w, squeeze(ephys_in_window(i,:)), 'LineWidth', 1 );
            ylim([-8 14]);
            
            linkaxes(ax2, 'x');
            
            saveas(f, [ analysis_path '/stim_' num2str(i) '_EB_yaw_ephys_sid_' num2str( cur_sid ) '.fig']);
            saveas(f, [ analysis_path '/stim_' num2str(i) '_EB_yaw_ephys_sid_' num2str( cur_sid ) '.png']);
            close(f);
        end
    end
end
end
