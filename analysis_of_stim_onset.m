function analysis_of_stim_onset( exp_dirs, VPS, analysis_path )

% exp_dirs = { {[ pad '181022_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15/analysis/'], 0 }, ...
%              {[ pad '181022_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15/analysis/'], 1 }, ...
%              {[ pad '181203_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16/analysis/'], 0 }, ...
%              {[ pad '181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17/analysis/'], 0 }, ...
%              {[ pad '181211_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18/analysis/'], 1 } };

settings = sensor_settings;
ephys_FR = settings.sampRate;
ball_FR = settings.sensorPollFreq;
BEFORE_STIM_TIME = 2.0;
bump_baseline_idx = [ 1 : VPS * BEFORE_STIM_TIME ];

for d = 1:length(exp_dirs)
    cur_dir              = exp_dirs{d}{1};
    cur_sid              = exp_dirs{d}{2};
    cur_stims_to_include = exp_dirs{d}{3};
    
    cur_EB_yaw_ephys_data = load( [cur_dir '/bump_yaw_ephys_in_window_data_sid_' num2str(cur_sid) '.mat' ] );
    
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
    
    bump_win  = [];
    ephys_win = [];
    yaw_win   = [];
    
    for s = 1:length( cur_stims_to_include )
        cur_stim_id = cur_stims_to_include( s, 1 );
        % cur_start_t = cur_stims_to_include( s, 2 ) + 2.0; % Roughly around the time of the bump turn after start
        cur_start_t = cur_stims_to_include( s, 2 ); 
        cur_end_t   = cur_stims_to_include( s, 3 );
             
        disp([ 'About to process stim: ' num2str(s) ]);

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
        
        cur_EB_bump_in_interval = cur_EB_bump( EB_interval );        
       
        cur_yaw_in_interval     = cur_yaw( yaw_interval );
        
        % downsample yaw to be in the same time base as EB
        dT = floor(ball_FR / VPS);

        min_len = floor(length(cur_yaw_in_interval)/dT);
        dT_len = dT*min_len;
        
        cur_yaw_in_interval_down = mean(reshape( cur_yaw_in_interval(1:dT_len), [dT, min_len ]));
        
        cur_ephys_in_interval   = cur_ephys( ephys_interval );        
      
        bump_to_avg{end+1}  = cur_EB_bump_in_interval; 
        yaw_to_avg{end+1}   = cur_yaw_in_interval;
        ephys_to_avg{end+1} = cur_ephys_in_interval; 
        
        EB_bump_pos_filt = medfilt1( cur_EB_bump_in_interval, 10, 'truncate' );
                
        dt = 1 / VPS;
        EB_vel = diff(EB_bump_pos_filt) / dt;  
        
        if 1

            subplot(2,1,1);
            hold on;
            imagesc( t_bump_w(EB_interval), [1:size(cur_bump_rois,1)], cur_bump_rois(:,EB_interval) );
            colormap(flipud(gray));
            axis tight;
            ylim([0 9]);
            caxis([-0.5 4]);
            title(['Stim: ' num2str(s)]);
            plot( t_bump_w(EB_interval), EB_bump_pos_filt );    
            
            subplot(2,1,2);
            hold on;
            plot( t_bump_w(EB_interval(1:end-1)), EB_vel );
             
            xlabel('Time (s)');
            ylabel('EB vel');
            
            saveas(f, [analysis_path '/EB_bump_location_and_velocity_stimID_' num2str( s ) '.fig']);
            saveas(f, [analysis_path '/EB_bump_location_and_velocity_stimID_' num2str( s ) '.png']);
            % close(f);
            waitforbuttonpress;            
            subplot(2,1,1); cla();
            subplot(2,1,2); cla();
        end
        
        if 0
            hold on;
            EB_speed = abs(EB_vel_filt(1:min_len));
            yaw_speed = abs( cur_yaw_in_interval_down(1:min_len) );
            scatter( EB_speed, yaw_speed, [], rgb('DarkGray'), 'LineWidth', 0.8 );
            xlabel('EB bump velocity');
            ylabel('Yaw velocity');
        end
        
        if 0
            yyaxis left;
            plot( t_bump_w(EB_interval(1:end-1)), EB_vel_filt );
            ylabel('EB Vel (au/s)');
            
            yyaxis right;
            plot( t_yaw_w(yaw_interval), cur_yaw_in_interval );
            ylabel('Yaw (au/s)');
            waitforbuttonpress;
        end
        
        if 0
            % Calculate bump speed
            subplot(3,1,1);
            hold on;
            plot( cur_EB_bump_in_interval );
            ylabel('EB position');
            
            subplot(3,1,2);
            hold on;
            
            dt = 1 / VPS;
            EB_vel = diff(cur_EB_bump_in_interval) / dt;
            %EB_vel_filt = medfilt1( EB_vel, 10, 'truncate' );
            
            plot( EB_vel );
            ylabel('EB Vel (au/s)');
            % ylim([-100 100]);
            
            subplot(3,1,3);
            hold on;
            plot( cur_yaw_in_interval );
            ylabel('Yaw (au/s)');
        end
    end           
    
    if 0
        % Plot average responses
        bump_max_len = -1;
        for i = 1:length( bump_to_avg )
            if( length(bump_to_avg{i}) > bump_max_len )
                bump_max_len = length( bump_to_avg{i} );
            end
        end
        
        yaw_max_len = -1;
        for i = 1:length( yaw_to_avg )
            if( length(yaw_to_avg{i}) > yaw_max_len )
                yaw_max_len = length(yaw_to_avg{i});
            end
        end
        
        ephys_max_len = -1;
        for i = 1:length( ephys_to_avg )
            if( length(ephys_to_avg{i}) > ephys_max_len )
                ephys_max_len = length(ephys_to_avg{i});
            end
        end
        
        bump_all = zeros( length( bump_to_avg ), bump_max_len );
        yaw_all = zeros( length( yaw_to_avg ), yaw_max_len );
        ephys_all = zeros( length( ephys_to_avg ), ephys_max_len );
        
        f = figure;
        for i = 1:length( yaw_to_avg )
            bump_all( i, 1:length(bump_to_avg{i}) ) = bump_to_avg{i};
            yaw_all( i, 1:length(yaw_to_avg{i}) ) = yaw_to_avg{i};
            ephys_all( i, 1:length(ephys_to_avg{i}) ) = ephys_to_avg{i};
            
            subplot(3,1,1);
            hold on;
            plot( bump_all( i, 1:length(bump_to_avg{i}) ) );
            ylabel('EB position');
            
            subplot(3,1,2);
            hold on;
            plot( yaw_all( i, 1:length(yaw_to_avg{i}) ) );
            ylabel('Yaw (au/s)');
            
            subplot(3,1,3);
            hold on;
            plot( ephys_all( i, 1:length(ephys_to_avg{i}) ) );
            ylabel('Vm (mV)');
        end
        
        subplot(3,1,1);
        hold on;
        plot( mean(bump_all), 'LineWidth', 3.0 );
        
        subplot(3,1,2);
        hold on;
        plot( mean(yaw_all), 'LineWidth', 3.0 );
        
        subplot(3,1,3);
        hold on;
        plot( mean(ephys_all), 'LineWidth', 3.0 );
    end
end

end