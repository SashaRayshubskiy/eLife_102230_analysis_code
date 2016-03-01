function generate_df_f_vs_turning_magnitude_trial_plot( sid, rois, PLANE, cdata_raw, turn_metadata, frame_start_offsets_per_plane, VPS, btrial_meta, analysis_path )

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

    for trial_ord_idx = 1:size( turn_metadata{ trial_type }, 1 )        
        
        % Collect calcium data
        trial_cdata = squeeze( cdata_raw{ trial_type }( trial_ord_idx, :, :, PLANE, : ) );
        
        left_df_f  = get_df_f_in_roi( trial_cdata, base_begin, base_end, rois{ 1 } );
        right_df_f = get_df_f_in_roi( trial_cdata, base_begin, base_end, rois{ 2 } );
        
        left_max_stim_df_f = mean( left_df_f( stim_frames ) );
        right_max_stim_df_f = mean( right_df_f( stim_frames ) );
        
        turn_magnitude = turn_metadata{ trial_type }( trial_ord_idx, 2 );        
        
        subplot(2,2,trial_type)
        hold on;
        plot( turn_magnitude, left_max_stim_df_f, 'X' );
        drawnow();

        subplot(2,2,trial_type+2)
        hold on;
        plot( turn_magnitude, right_max_stim_df_f, 'X' );              
        drawnow();
        disp(['Processed: ' ac.task_str{ trial_type } ' trial: ' num2str( trial_ord_idx ) ]);
    end
    
    subplot(2,2,trial_type)
    title( [ ac.task_str{ trial_type } ': turn magnitude vs. left hemisphere dF/F']);
    ylabel('dF/F');
    xlabel('Turn magnitude');
    yy = ylim;
    ylim([0 yy(2)]);
    
    subplot(2,2,trial_type+2)
    title( [ ac.task_str{ trial_type } ': turn magnitude vs. right hemisphere dF/F']);
    ylabel('dF/F');
    xlabel('Turn magnitude');
    yy = ylim;
    ylim([0 yy(2)]); 
    
    pause(0.01); % This is so I can cancel if the program runs for too long.
end

savepath = [ analysis_path '/df_f_vs_turn_magnitude_sid_' num2str( sid ) ];
saveas( f, [savepath '.png'] );
saveas( f, [savepath '.fig'] );

end

