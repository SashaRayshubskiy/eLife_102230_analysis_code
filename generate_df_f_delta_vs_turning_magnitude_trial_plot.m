function generate_df_f_delta_vs_turning_magnitude_trial_plot( sid, rois, PLANE, cdata_raw, turn_metadata, frame_start_offsets_per_plane, VPS, btrial_meta, analysis_path )

ac = get_analysis_constants;
settings = sensor_settings;

prestim = settings.pre_stim;
stim    = settings.stim;

base_begin = 1;
baseline_time = prestim;
base_end = floor(baseline_time*VPS);

nframes = size(cdata_raw{ 1 }(1,:,:,PLANE,:),5);
t = (([0:nframes-1]))./VPS + frame_start_offsets_per_plane( PLANE );

stim_frames = find((t >= prestim) & (t<=(prestim+stim)));

f = figure;
for trial_type = 1:2
    trial_cnt = size( turn_metadata{ trial_type }, 1 );
    
    x = [];
    y = [];
    
    for trial_ord_idx = 1:trial_cnt
        
        % Collect calcium data
        trial_cdata = squeeze( cdata_raw{ trial_type }( trial_ord_idx, :, :, PLANE, : ) );
        
        left_df_f  = get_df_f_in_roi( trial_cdata, base_begin, base_end, rois{ 1 } );
        right_df_f = get_df_f_in_roi( trial_cdata, base_begin, base_end, rois{ 2 } );
        
        left_max_stim_df_f = mean( left_df_f( stim_frames ) );
        right_max_stim_df_f = mean( right_df_f( stim_frames ) );
        
        turn_magnitude = turn_metadata{ trial_type }( trial_ord_idx, 2 );        
        
        subplot(1,2,trial_type)
        hold on;
        if (trial_type == ac.LEFT )
            yy = left_max_stim_df_f-right_max_stim_df_f;
        else
            yy = right_max_stim_df_f-left_max_stim_df_f;            
        end
        
        if( abs(turn_magnitude) < 0.03 )
            continue;
        end
        
        y(end+1) = yy;
        x(end+1) = turn_magnitude;
        
        plot( turn_magnitude, yy, 'X' );
        drawnow();

        disp(['Processed: ' ac.task_str{ trial_type } ' trial: ' num2str( trial_ord_idx ) ]);
    end
    
    subplot(1,2,trial_type)    
    ylabel('dF/F');
    xlabel('Turn magnitude');
    
    p = polyfit(x,y,1);
    xl = linspace(min(x), max(x), length(y));
    yfit = polyval(p,xl);
    plot(xl, yfit, 'r--', 'LineWidth', 3.0);
    
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    
    rsq = 1 - SSresid/SStotal;
    
    title( [ ac.task_str{ trial_type } ': turn magnitude vs. delta dF/F -- R^2: ' num2str(rsq)]);
    
    pause(0.01); % This is so I can cancel if the program runs for too long.
end

savepath = [ analysis_path '/df_f_delta_vs_turn_magnitude_sid_' num2str( sid ) ];
saveas( f, [savepath '.png'] );
saveas( f, [savepath '.fig'] );

end

