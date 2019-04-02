function [t_win, cx_vars] = analysis_of_bump_return_V2( exp_dirs, VPS, direction_of_return )

settings = sensor_settings;
ephys_FR = settings.sampRate;
ball_FR = settings.sensorPollFreq;
BEFORE_STIM_TIME = 2.0;
bump_baseline_idx = [ 1 : VPS * BEFORE_STIM_TIME ];

bump_win_all = [];
yaw_win_all = [];
ephys_win_all = [];

for d = 1:length(exp_dirs)
    cur_dir              = exp_dirs{ d }{ 1 };
    cur_sid              = exp_dirs{ d }{ 2 };
    cur_stims_to_include = exp_dirs{ d }{ 3 };
    bump_return_times    = exp_dirs{ d }{ 4 };
    
    analysis_path = cur_dir;
    cur_EB_yaw_ephys_data = load( [analysis_path '/bump_yaw_ephys_in_window_data_sid_' num2str(cur_sid) '.mat' ] );
        
    cur_yaw_in_window   = cur_EB_yaw_ephys_data.yaw_in_window;
    cur_ephys_in_window = cur_EB_yaw_ephys_data.ephys_in_window;
    cur_bump_in_window  = cur_EB_yaw_ephys_data.bump_in_window;
    
    t_bump_w   = ([0:size(cur_bump_in_window,3)-1] / VPS) - BEFORE_STIM_TIME;
    t_yaw_w    = ([0:size(cur_yaw_in_window,2)-1] / ball_FR) - BEFORE_STIM_TIME;
    t_ephys_w  = ([0:size(cur_ephys_in_window,2)-1] / ephys_FR) - BEFORE_STIM_TIME;
    
    assert( length(cur_stims_to_include) > 0 );
    
    f = figure;
    
    bump_to_avg = {};
    yaw_to_avg = {};
    ephys_to_avg = {};
    
    j = 1;
    bump_win  = [];
    ephys_win = [];
    yaw_win   = [];
    
    stim_ids = int32(squeeze(bump_return_times(:,1)));
    
    while( j <= length(stim_ids) )
        
        s = stim_ids(j);
        cur_stim_id = cur_stims_to_include( s, 1 );
        % cur_start_t = cur_stims_to_include( s, 2 ) + 2.0; % Roughly around the time of the bump turn after start
        cur_start_t = cur_stims_to_include( s, 2 );
        cur_end_t   = cur_stims_to_include( s, 3 );
        
        EB_bump_vel_align = bump_return_times( j, 2 ) + cur_start_t;
        j = j + 1;
        
        % disp([ 'About to process stim: ' num2str(s) ]);
        
        cur_bump_rois = squeeze(cur_bump_in_window(cur_stim_id,:,:));
        
        cur_EB_bump = get_radial_weighted_avg_bump_pos_v3( cur_bump_rois );
        if( length(cur_EB_bump) == 0 )
            % Skip this stim, because there was impossible to get a clean
            % bump trajectory, probably due to sporadic bump disappearance
            
            disp([ 'Rejecting stim: ' num2str(s) ]);
            continue;
        end
        
        % cur_EB_bump = medfilt1( cur_bump_tc - mean(cur_bump_tc(bump_baseline_idx)), 5, 'truncate' );
        cur_yaw     = cur_yaw_in_window( cur_stim_id, : );
        cur_ephys   = cur_ephys_in_window( cur_stim_id, : );
        
        EB_interval = find( (t_bump_w >= cur_start_t) & (t_bump_w <= cur_end_t));
        yaw_interval = find( (t_yaw_w >= cur_start_t) & (t_yaw_w <= cur_end_t));
        ephys_interval = find( (t_ephys_w >= cur_start_t) & (t_ephys_w <= cur_end_t));
        
        % cur_EB_bump_in_interval = cur_EB_bump( EB_interval );
        cur_EB_bump_in_interval = cur_EB_bump;
        
        %cur_yaw_in_interval     = cur_yaw( yaw_interval );
        cur_yaw_in_interval     = cur_yaw;
        
        % downsample yaw to be in the same time base as EB
        dT = floor(ball_FR / VPS);
        
        min_len = floor(length(cur_yaw_in_interval)/dT);
        dT_len = dT*min_len;
        
        cur_yaw_in_interval_down = mean(reshape( cur_yaw_in_interval(1:dT_len), [dT, min_len ]));
        
        %cur_ephys_in_interval   = cur_ephys( ephys_interval );
        cur_ephys_in_interval   = cur_ephys;
        
        bump_to_avg{end+1}  = cur_EB_bump_in_interval;
        yaw_to_avg{end+1}   = cur_yaw_in_interval;
        ephys_to_avg{end+1} = cur_ephys_in_interval;
        
        EB_bump_pos_filt = medfilt1( cur_EB_bump_in_interval, 10, 'truncate' );
        
        dt = 1 / VPS;
        EB_vel = diff(EB_bump_pos_filt) / dt;
        % EB_vel_filt = medfilt1( EB_vel, 10, 'truncate' );
        
        % Align by the EB vel
        xx = find( t_bump_w < EB_bump_vel_align );
        EB_bump_vel_align_idx = xx(end);
        
        TIME_BEFORE_EB_VEL_CHANGE = 1.0;
        TIME_AFTER_EB_VEL_CHANGE = 2.0;
        
        EB_FRAMES_BEFORE_EB_VEL_CHANGE = floor( TIME_BEFORE_EB_VEL_CHANGE * VPS );
        EB_FRAMES_AFTER_EB_VEL_CHANGE  = floor( TIME_AFTER_EB_VEL_CHANGE * VPS );
        cur_EB_bump_vel_win_start  = EB_bump_vel_align_idx - EB_FRAMES_BEFORE_EB_VEL_CHANGE;
        cur_EB_bump_vel_win_end  = EB_bump_vel_align_idx + EB_FRAMES_AFTER_EB_VEL_CHANGE;
        
        if( ( cur_EB_bump_vel_win_start < 1) || ( cur_EB_bump_vel_win_end > length(EB_vel) ) )
            continue;
        end
        
        cur_EB_vel_win = EB_vel( cur_EB_bump_vel_win_start:cur_EB_bump_vel_win_end );
        bump_win(end+1,:) = cur_EB_vel_win;
        t_bump_win = t_bump_w( cur_EB_bump_vel_win_start:cur_EB_bump_vel_win_end ) - t_bump_w(EB_bump_vel_align_idx);
        
        % Yaw
        % t_yaw_interval = t_yaw_w(yaw_interval);
        xx = find( t_yaw_w < EB_bump_vel_align );
        yaw_align_idx = xx(end);
        
        YAW_FRAMES_BEFORE_EB_VEL_CHANGE = floor( TIME_BEFORE_EB_VEL_CHANGE * ball_FR );
        YAW_FRAMES_AFTER_EB_VEL_CHANGE  = floor( TIME_AFTER_EB_VEL_CHANGE * ball_FR );
        cur_yaw_win_start  = yaw_align_idx - YAW_FRAMES_BEFORE_EB_VEL_CHANGE;
        cur_yaw_win_end  = yaw_align_idx + YAW_FRAMES_AFTER_EB_VEL_CHANGE;
        
        cur_yaw_win = cur_yaw_in_interval( cur_yaw_win_start:cur_yaw_win_end );
        yaw_win(end+1,:) = cur_yaw_win;
        t_yaw_win = t_yaw_w( cur_yaw_win_start:cur_yaw_win_end ) - t_yaw_w(yaw_align_idx);
        
        % Ephys
        % t_ephys_interval = t_ephys_w(ephys_interval);
        xx = find( t_ephys_w < EB_bump_vel_align );
        ephys_align_idx = xx(end);
        
        EPHYS_FRAMES_BEFORE_EB_VEL_CHANGE = floor( TIME_BEFORE_EB_VEL_CHANGE * ephys_FR );
        EPHYS_FRAMES_AFTER_EB_VEL_CHANGE  = floor( TIME_AFTER_EB_VEL_CHANGE * ephys_FR );
        cur_ephys_win_start  = ephys_align_idx - EPHYS_FRAMES_BEFORE_EB_VEL_CHANGE;
        cur_ephys_win_end  = ephys_align_idx + EPHYS_FRAMES_AFTER_EB_VEL_CHANGE;
        
        cur_ephys_win = cur_ephys_in_interval( cur_ephys_win_start:cur_ephys_win_end );
        ephys_win(end+1,:) = cur_ephys_win;
        t_ephys_win = t_ephys_w( cur_ephys_win_start:cur_ephys_win_end ) - t_ephys_w(ephys_align_idx);
        
        subplot(2,1,1);
        yyaxis left;
        hold on;
        plot( t_bump_win, cur_EB_vel_win, '-' );
        ylim([-15 15]);
        ylabel('EB Vel (au/s)');
        
        yyaxis right;
        hold on;
        plot( t_yaw_win, cur_yaw_win, '-' );
        ylabel('Yaw (au/s)');
        xlabel('Time (s)');
        
        subplot(2,1,2);
        yyaxis left;
        hold on;
        plot( t_bump_win, cur_EB_vel_win, '-' );
        ylim([-15 15]);
        ylabel('EB Vel (au/s)');
        
        yyaxis right;
        hold on;
        plot( t_ephys_win, cur_ephys_win, '-' );
        ylabel('Vm (mV)');
        xlabel('Time (s)');
        % waitforbuttonpress;
    end
    
    avg_bump_win = mean(bump_win);
    bump_win_all(d, :) = avg_bump_win;

    avg_yaw_win = mean(yaw_win);
    yaw_win_all(d, :) = avg_yaw_win;    

    avg_ephys_win = mean(ephys_win);
    ephys_win_all(d, :) = avg_ephys_win;    
    
    subplot(2,1,1);
    yyaxis left;
    hold on;
    plot( t_bump_win, avg_bump_win, '-', 'LineWidth', 3.0 );
    ylim([-15 15]);
    ylabel('EB Vel (au/s)');
    
    yyaxis right;
    hold on;
    plot( t_yaw_win, avg_yaw_win, '-', 'LineWidth', 3.0, 'color', rgb('Red') );
    ylabel('Yaw (au/s)');
    xlabel('Time (s)');
    
    axis tight;
    
    subplot(2,1,2);
    yyaxis left;
    hold on;
    plot( t_bump_win, avg_bump_win, '-', 'LineWidth', 3.0 );
    ylim([-15 15]);
    ylabel('EB Vel (au/s)');
    
    yyaxis right;
    hold on;
    plot( t_ephys_win, avg_ephys_win, '-', 'LineWidth', 3.0, 'color', rgb('Red') );
    ylabel('Vm (mV)');
    xlabel('Time (s)');
    axis tight;
    
    saveas( f, [analysis_path '/bump_return_vs_yaw_vs_ephys_' direction_of_return '.fig'] );
    saveas( f, [analysis_path '/bump_return_vs_yaw_vs_ephys_' direction_of_return '.png'] ); 
end

figure_path = '/data/drive2/sasha/CX_summary/';

f = figure;
% Plots for all flies
for d = 1:length(exp_dirs)
    ax1(1) = subplot(3,1,1);
    hold on;
    plot( t_bump_win, squeeze( bump_win_all(d,:) ), 'color', rgb('Gray'));    
    ylabel('EB Vel (au/s)');

    ax1(2) = subplot(3,1,2);
    hold on;
    plot( t_yaw_win, squeeze( yaw_win_all(d,:) ), 'color', rgb('Gray'));    
    ylabel('Yaw (au/s)');

    ax1(3) = subplot(3,1,3);
    hold on;
    plot( t_ephys_win, squeeze( ephys_win_all(d,:) ), 'color', rgb('Gray'));    
    ylabel('Yaw (au/s)');
    xlabel('Time (s)');
end

ax1(1) = subplot(3,1,1);
plt_hdl1 = plot( t_bump_win, mean( bump_win_all ), 'color', rgb('Red'), 'LineWidth', 3, 'DisplayName', ['n = ' num2str(length(exp_dirs))] );
ylabel('EB Vel (au/s)');
title( [ 'Bump direction: ' direction_of_return ] );

ax1(2) = subplot(3,1,2);
plot( t_yaw_win, mean( yaw_win_all ), 'color', rgb('Red'), 'LineWidth', 3 );
ylabel('Yaw (au/s)');

ax1(3) = subplot(3,1,3);
plot( t_ephys_win, mean( ephys_win_all ), 'color', rgb('Red'), 'LineWidth', 3 );
ylabel('Yaw (au/s)');
xlabel('Time (s)');

legend(plt_hdl1, ['n = ' num2str(length(exp_dirs))]);
linkaxes(ax1, 'x');
saveas( f, [figure_path '/bump_return_' direction_of_return '_all_flies.fig'] );
saveas( f, [figure_path '/bump_return_' direction_of_return '_all_flies.png'] );

t_win = cell(1,3);
cx_vars = cell(1,3);

t_win{1} = t_bump_win;
t_win{2} = t_yaw_win;
t_win{3} = t_ephys_win;

cx_vars{1} = bump_win_all;
cx_vars{2} = yaw_win_all;
cx_vars{3} = ephys_win_all;

end