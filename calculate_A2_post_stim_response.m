function [ bump_response, fwd_response, yaw_response, ephys_response ] = calculate_A2_post_stim_response( basedir, data_dirs, experiment_type )

settings = sensor_settings;
BALL_FR = settings.sensorPollFreq;
EPHYS_FR = settings.sampRate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BEFORE_STIM_TIME = 2.0; % seconds
AFTER_STIM_TIME  = 10.0; % seconds

BASELINE_WINDOW_SIZE = 0.5; % seconds
BASELINE_WINDOW_MAX = 0.0;
BASELINE_WINDOW_MIN = BASELINE_WINDOW_MAX - BASELINE_WINDOW_SIZE;

POST_STIM_WINDOW_MIN = 0.0;
% POST_STIM_WINDOW_MAX = 2.0;
POST_STIM_WINDOW_MAX = 0.85;


DEBUG_SHOW_BUMP_TC = 11;
DEBUG_OFF = 12;
DEBUG_LEVEL = DEBUG_OFF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

color_per_fly = { rgb('Violet'), rgb('Red'), rgb('Black'), rgb('Blue') };
f = figure('units','normalized','outerposition',[0 0 1 1]);    

bump_response = [];
fwd_response = [];
yaw_response = [];
ephys_response = [];

for d = 1:length( data_dirs )
    
%    JUMP_THRESHOLD               = data_dirs{d}{4};
%    BUMP_RETURN_STABILITY_WINDOW = data_dirs{d}{5};
    
    cur_datapath = data_dirs{d}{1};
    cur_sid      = data_dirs{d}{2};
    
    datapath = [ basedir '/' cur_datapath ]; 
    analysis_path = [datapath '/analysis/'];
    stim_filepath = [analysis_path '/stim_window_data_sid_' num2str(cur_sid) '.mat'];
    
    cur_data = load( stim_filepath );
    
    raw_data_filepath = [analysis_path '/raw_glom_ball_ephys.mat'];
    VPS_tmp = load( raw_data_filepath );
    VPS = VPS_tmp.VPS;
    clear VPS_tmp;    
    
    bump_data   = cur_data.bump_in_window;
    fwd_data    = cur_data.fwd_in_window;
    yaw_data    = convert_yaw_to_degrees( cur_data.yaw_in_window );
    ephys_data  = cur_data.ephys_in_window;
            
    t_bump_w   = ([0:size(bump_data,3)-1] / VPS) - BEFORE_STIM_TIME;
    t_yaw_w    = ([0:size(yaw_data,2)-1] / BALL_FR) - BEFORE_STIM_TIME;
    t_ephys_w  = ([0:size(ephys_data,2)-1] / EPHYS_FR) - BEFORE_STIM_TIME;
    
    bump_baseline  = find( (t_bump_w <= BASELINE_WINDOW_MAX) & (t_bump_w > BASELINE_WINDOW_MIN) );
    yaw_baseline   = find( (t_yaw_w <= BASELINE_WINDOW_MAX) & (t_yaw_w > BASELINE_WINDOW_MIN) );
    ephys_baseline = find( (t_ephys_w <= BASELINE_WINDOW_MAX) & (t_ephys_w > BASELINE_WINDOW_MIN) );

    bump_stim_window   = find( (t_bump_w <= POST_STIM_WINDOW_MAX) & (t_bump_w > POST_STIM_WINDOW_MIN) );
    yaw_stim_window   = find( (t_yaw_w <= POST_STIM_WINDOW_MAX) & (t_yaw_w > POST_STIM_WINDOW_MIN) );
    ephys_stim_window = find( (t_ephys_w <= POST_STIM_WINDOW_MAX) & (t_ephys_w > POST_STIM_WINDOW_MIN) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collect an average of all trials post stim
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    bump_avg_per_stim = [];
    for st = 1:size( bump_data, 1 )
        
        [ smoothed_bump, cur_bump_tc, cur_bump_tc_unwrapped, vect_strength_check ] = get_radial_weighted_avg_bump_pos_vect_strengh_check_v2( squeeze( bump_data( st, :, : ) ) );
                       
        % Rotate bump data so that the average pre-stim period is in the same
        % place for each trial.
        BUMP_TC_FILTER_SAMPLE_POINTS = 5;
        baseline_vals = cur_bump_tc( bump_baseline );
        baseline_non_nan = baseline_vals(~isnan(baseline_vals));        
        bump_delta_tc = medfilt1( cur_bump_tc - mean(baseline_non_nan), BUMP_TC_FILTER_SAMPLE_POINTS, 'truncate' );
                
        bump_avg_per_stim(st, : ) = bump_delta_tc;
    end

    avg_bump = mean(bump_avg_per_stim);
    bump_response( d ) = mean( avg_bump( bump_stim_window ) ) - mean( avg_bump( bump_baseline ) );
    
    avg_fwd = mean(fwd_data);
    fwd_response( d ) = mean( avg_fwd( yaw_stim_window ) ) - mean( avg_fwd( yaw_baseline ) );

    avg_yaw = mean(yaw_data);
    yaw_response( d ) = mean( avg_yaw( yaw_stim_window ) ) - mean( avg_yaw( yaw_baseline ) );
    
    ephys_spikes_removed = [];
    for st = 1:size( bump_data, 1 )
        FILT_FACTOR = 0.04;
        cur_Vm_w_drift = medfilt1( squeeze(ephys_data(st,:)), FILT_FACTOR * EPHYS_FR, 'truncate' );
        cur_Vm = cur_Vm_w_drift - mean( cur_Vm_w_drift );

        ephys_spikes_removed(st,:) = cur_Vm;
    end
    avg_ephys = mean(ephys_spikes_removed);
    ephys_response(d) = mean( avg_ephys( ephys_stim_window ) ) - mean( avg_ephys( ephys_baseline ) );        
    

    % Display results
    subplot(4,1,1);
    hold on;
    plot( t_bump_w, avg_bump, 'color', color_per_fly{d} );
    plot( t_bump_w(bump_stim_window(end)), bump_response(d), 'X', 'color', color_per_fly{d} );

    subplot(4,1,2);
    hold on;
    plot( t_yaw_w, avg_fwd, 'color', color_per_fly{d} );
    plot( t_yaw_w( yaw_stim_window(end) ), fwd_response(d), 'X', 'color', color_per_fly{d} );

    subplot(4,1,3);
    hold on;
    plot( t_yaw_w, avg_yaw, 'color', color_per_fly{d} );
    plot( t_yaw_w( yaw_stim_window(end)), yaw_response(d), 'X', 'color', color_per_fly{d} );

    subplot(4,1,4);
    hold on;
    plot( t_ephys_w, avg_ephys, 'color', color_per_fly{d} );
    plot( t_ephys_w( ephys_stim_window(end)), ephys_response(d), 'X', 'color', color_per_fly{d} );        
end

saveas( f, [analysis_path '/' cur_datapath '_post_stim_avg.fig'] );
saveas( f, [analysis_path '/' cur_datapath '_post_stim_avg.png'] );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
