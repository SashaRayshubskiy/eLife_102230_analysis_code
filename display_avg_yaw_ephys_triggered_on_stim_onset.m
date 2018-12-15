function [stim_events] = display_avg_yaw_ephys_triggered_on_stim_onset( df_f_in_roi_per_trial, pico_stim_data, ephys_time, ephys_data, bdata_vel_time, bdata_vel, VPS, analysis_path, sid, stims_to_include )

BEFORE_STIM_TIME = 2.0; % seconds
AFTER_STIM_TIME  = 10.0; % seconds

STIM_TIME_WIN = 2.0; % KEEP THIS CONSTANT, OTHERWISE STIM IDs change

EPHYS_FR = 10000;
EPHYS_BEFORE_STIM_FRAMES = EPHYS_FR*BEFORE_STIM_TIME;
EPHYS_AFTER_STIM_FRAMES = EPHYS_FR*AFTER_STIM_TIME;

BALL_FR = 100;
BALL_BEFORE_STIM_FRAMES = BALL_FR*BEFORE_STIM_TIME;
BALL_AFTER_STIM_FRAMES  = BALL_FR*AFTER_STIM_TIME;

yaw_in_window = [];
ephys_in_window = [];
bump_in_window = [];
ac = get_analysis_constants;

for tr = 1:size(pico_stim_data{1},1)
    cur_stim = 10.0*squeeze(pico_stim_data{1}( tr, : ));

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
            stim_idx_keep(end+1) = stim_idx(st);
        end
        
        prev_st = cur_st;
    end
    
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
              
        yaw_in_window(end+1,:)    = squeeze( bdata_vel{ 1 }( tr, ac.VEL_YAW, before_ball_frame:after_ball_frame ) );
        ephys_in_window(end+1,:)  = squeeze( ephys_data(tr, begin_ephys_frame:end_ephys_frame) );
        bump_in_window(end+1,:,:) = squeeze( df_f_in_roi_per_trial( tr, :, begin_bump_frame:end_bump_frame ) );
    end
end

save([analysis_path '/bump_yaw_ephys_in_window_data_sid_' num2str(sid) '.mat'], 'yaw_in_window', 'ephys_in_window', 'bump_in_window' );

% Display averages in a window
f = figure;
bump_tc_all = [];

t_bump_w  = ([0:size(bump_in_window,3)-1] / VPS) - BEFORE_STIM_TIME;
t_yaw_w  = ([0:size(yaw_in_window,2)-1] / BALL_FR) - BEFORE_STIM_TIME;
t_ephys_w  = ([0:size(ephys_in_window,2)-1] / EPHYS_FR) - BEFORE_STIM_TIME;

bump_baseline_idx = [ 1 : VPS * BEFORE_STIM_TIME ];

if( isempty( stims_to_include ) == 0 )
    stim_ids_to_include = stims_to_include(:,1)';
else
    stim_ids_to_include = [1:size( yaw_in_window, 1 )];
end


for i = stim_ids_to_include
    
    subplot(3,1,1);
    hold on;
    
    % [dummy, cur_bump_tc] = max( squeeze(bump_in_window(i,:,:)));      
    cur_bump_tc = get_radial_weighted_avg_bump_pos( squeeze(bump_in_window(i,:,:)) );
    
    % Rotate bump data so that the average pre-stim period is in the same
    % place for each trial.
    bump_delta_tc = medfilt1( cur_bump_tc - mean(cur_bump_tc(bump_baseline_idx)), 5, 'truncate' );

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

saveas( f, [analysis_path '/avg_EB_yaw_ephys_triggered_on_stim_' num2str(sid) '.fig'] );
saveas( f, [analysis_path '/avg_EB_yaw_ephys_triggered_on_stim_' num2str(sid) '.png'] );

if 1
    % Display all individual stim events
    for i = 1:size( yaw_in_window, 1 )
%    for i = 1
    
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

        saveas(f, [ analysis_path '/stim_' num2str(i) '_EB_yaw_ephys_sid_' num2str(sid) '.fig']);
        saveas(f, [ analysis_path '/stim_' num2str(i) '_EB_yaw_ephys_sid_' num2str(sid) '.png']);
        close(f);
    end   
end

% Save a matrix of EB bump location, yaw and ephys for piecewise time warping
if( isempty(stims_to_include) == 0 )
    
    % Create an array for the largest span of time
    max_t_range = max(stims_to_include(:,3) - stims_to_include(:,2));
    stim_event_length = (max_t_range + 1) * EPHYS_FR;
    
    stim_events.EB    = cell( size(stims_to_include, 1) );
    stim_events.yaw   = cell( size(stims_to_include, 1) );
    stim_events.ephys = cell( size(stims_to_include, 1) );
        
    for i = 1:size(stims_to_include,1)
        cur_stim_id = stims_to_include( i, 1 );
        cur_start_t = stims_to_include( i, 2 );
        cur_end_t   = stims_to_include( i, 3 );

        cur_EB_bump = bump_tc_all( i, : );
        cur_yaw     = yaw_in_window( cur_stim_id, : );
        cur_ephys   = ephys_in_window( cur_stim_id, : );
        
        % Extract the time interval
        EB_interval = find( (t_bump_w >= cur_start_t) & (t_bump_w <= cur_end_t));
        yaw_interval = find( (t_yaw_w >= cur_start_t) & (t_yaw_w <= cur_end_t));
        ephys_interval = find( (t_ephys_w >= cur_start_t) & (t_ephys_w <= cur_end_t));
        
        cur_EP_bump_in_interval = cur_EB_bump( EB_interval );
        cur_yaw_in_interval     = cur_yaw( yaw_interval );
        cur_ephys_in_interval   = cur_ephys( ephys_interval );
        
        % Upsample EB and yaw to bring in the same timebase as ephys
        upscale_interval_EB = length(cur_EP_bump_in_interval)/length(cur_ephys_in_interval);
        cur_EB_bump_upsample = interp1( [1:length(cur_EP_bump_in_interval)],  cur_EP_bump_in_interval, [ 1 : upscale_interval_EB : length(cur_EP_bump_in_interval) ] );

        upscale_interval_yaw = length(cur_yaw_in_interval)/length(cur_ephys_in_interval);
        cur_yaw_upsample = interp1( [1:length(cur_yaw_in_interval)],  cur_yaw_in_interval, [ 1 : upscale_interval_yaw : length(cur_yaw_in_interval) ] );
        
        save_len = min( length(cur_EB_bump_upsample), min( length(cur_yaw_upsample), length(cur_ephys_in_interval)) );
        
        stim_events.EB{i}(1:save_len)    = cur_EB_bump_upsample(1:save_len) - mean(cur_EB_bump_upsample(1:save_len));
        stim_events.yaw{i}(1:save_len)   = cur_yaw_upsample(1:save_len);
        stim_events.ephys{i}(1:save_len) = cur_ephys_in_interval(1:save_len);
        
        % Not using upsample
%         stim_events.EB{i}    = cur_EP_bump_in_interval - mean(cur_EP_bump_in_interval);
%         stim_events.yaw{i}   = cur_yaw_in_interval;        
%         stim_events.ephys{i} = cur_ephys_in_interval;
    end
    
    
    % save([analysis_path '/stim_event_data.mat'], 'stim_events.EB', 'stim_events_yaw', 'stim_events_ephys');
else
    stim_events.EB = {};
    stim_events.yaw = {};
    stim_events.ephys = {};
end
end

