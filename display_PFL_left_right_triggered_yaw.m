function [] = display_PFL_left_right_triggered_yaw( PFL_LAL_dF_F_per_trial, bdata_vel, imaging_data_1, VPS, analysis_path )

ac = get_analysis_constants;
settings = sensor_settings;
ballSR = settings.sensorPollFreq;

num_trials = size(imaging_data_1, 1);
nframes = size(imaging_data_1, 5);

t = [0:nframes-1]./VPS;

good_trials_02 = [1:12];

good_trials = good_trials_02;

PFL_ASYMMETRY_CUTOFF = 0.28;

EXAMINE_LEFT_BIAS = 11;
EXAMINE_RIGHT_BIAS = 12;

WINDOW_BACK = 0.75;
WINDOW_FRONT = 0.75;

WINDOW_IMAGING_SAMPLES_BACK  = floor(WINDOW_BACK * VPS);
WINDOW_IMAGING_SAMPLES_FRONT = floor(WINDOW_FRONT * VPS);

WINDOW_BALL_SAMPLES_BACK  = WINDOW_BACK * ballSR;
WINDOW_BALL_SAMPLES_FRONT = WINDOW_FRONT * ballSR;

f = figure;
cnt = 1;

stop_here = 1;

for examine_cond = [EXAMINE_RIGHT_BIAS EXAMINE_LEFT_BIAS ]
    PFL_asymmetrical_activity_events = [];
    yaw_events = [];
    
    EXAMINE_BIAS = examine_cond;
    
    if( EXAMINE_BIAS == EXAMINE_RIGHT_BIAS )
        examine_bias_str = 'right';
    else
        examine_bias_str = 'left';
    end
            
    for tr = good_trials
            
        cur_PFL_LAL = squeeze(PFL_LAL_dF_F_per_trial(tr, :,:));
        
        left_PFLm_LAL  = squeeze( cur_PFL_LAL( 2, : ) );
        right_PFLm_LAL = squeeze( cur_PFL_LAL( 4, : ) );
        
        FFT_FILTER_CUTOFF = 1;
        left_PFLm_LAL_filt  = fft_filter( left_PFLm_LAL',  FFT_FILTER_CUTOFF, VPS );
        right_PFLm_LAL_filt = fft_filter( right_PFLm_LAL', FFT_FILTER_CUTOFF, VPS );
        
        left_right_PB_LAL_delta = left_PFLm_LAL_filt - right_PFLm_LAL_filt;
        left_right_PB_LAL_delta = left_right_PB_LAL_delta - mean(left_right_PB_LAL_delta);
        
        if( EXAMINE_BIAS == EXAMINE_RIGHT_BIAS )
            idx = find(left_right_PB_LAL_delta > -1.0 * PFL_ASYMMETRY_CUTOFF );
            left_right_PB_LAL_delta(idx) = 0;
            
            [~, locs] = findpeaks( -1.0*left_right_PB_LAL_delta, 'Annotate', 'extents' );
        else
            idx = find(left_right_PB_LAL_delta < PFL_ASYMMETRY_CUTOFF );
            left_right_PB_LAL_delta(idx) = 0;
            
            [~, locs] = findpeaks( left_right_PB_LAL_delta, 'Annotate', 'extents' );
        end
        
%         figure;
%         hold on;
%         plot(left_right_PB_LAL_delta);
%         plot(locs, left_right_PB_LAL_delta(locs),  'X', 'color', 'g' );
        
        for i = 1:length( locs )
            cur_loc = locs(i);
            
            imaging_win_before = cur_loc - WINDOW_IMAGING_SAMPLES_BACK;
            imaging_win_after = cur_loc + WINDOW_IMAGING_SAMPLES_FRONT;
            
            if( ( imaging_win_before >= 1 ) && ( imaging_win_after <= length( left_right_PB_LAL_delta ) ) )
                PFL_asymmetrical_activity_events( end+1, : ) = left_right_PB_LAL_delta( imaging_win_before : imaging_win_after );
            else
                continue;
                disp('Got here');
            end
            
            cur_t = t(cur_loc);
            
            cur_ball_idx = cur_t * ballSR;
            
            ball_win_before = cur_ball_idx - WINDOW_BALL_SAMPLES_BACK;
            ball_win_after = cur_ball_idx + WINDOW_BALL_SAMPLES_FRONT;
            
            cur_yaw = squeeze(bdata_vel{ 1 }( tr, ac.VEL_YAW, : ));
            cur_yaw = cur_yaw - mean(cur_yaw);
            
            if( ( ball_win_before >= 1 ) && ( ball_win_after <= length( cur_yaw ) ) )
                yaw_events( end+1, : ) = cur_yaw( ball_win_before : ball_win_after );
            end
        end
        
        if(i==length(locs))
            stop_here = stop_here + 1;
        end
    end
    
    ball_window_t = ([0:size(yaw_events,2)-1] / ballSR) - WINDOW_BACK;
    imaging_window_t = ([0:size(PFL_asymmetrical_activity_events,2)-1] / VPS) - WINDOW_BACK;
    
    subplot(1,2,cnt);
    yyaxis left;
    hold on;
    
    PFL_avg = mean(PFL_asymmetrical_activity_events);
    PFL_sem = get_sem(PFL_asymmetrical_activity_events,1);
    
    fh = fill( [imaging_window_t, fliplr(imaging_window_t)], ...
        [(PFL_avg-PFL_sem) fliplr((PFL_avg+PFL_sem))], ...
        rgb('Violet'));
    set( fh, 'EdgeColor', 'None' );
    
    plot( imaging_window_t, PFL_avg, 'color', rgb('DarkMagenta') );
    xlabel('Time from PFL peak(s)');
    ylabel('PFL left/right');
    
    yyaxis right;
    hold on;
    
    yaw_avg = mean(yaw_events);
    yaw_sem = get_sem(yaw_events,1);
    
    fh = fill( [ball_window_t, fliplr(ball_window_t)], ...
        [(yaw_avg-yaw_sem) fliplr((yaw_avg+yaw_sem))], ...
        rgb('PaleGreen'));
    set( fh, 'EdgeColor', 'None' );
    
    plot( ball_window_t, yaw_avg, 'color', rgb('SeaGreen') );
    xlabel('Time from PFL peak(s)');
    ylabel('Yaw vel (au/s)');
    title(['Number of events: ' num2str(size(PFL_asymmetrical_activity_events,1))] );
    grid on;
    
    cnt = cnt + 1;
end

saveas(f, [analysis_path '/PFLm.LAL_asymmetric_events_vs_yaw.fig']);
saveas(f, [analysis_path '/PFLm.LAL_asymmetric_events_vs_yaw.png']);
end