function analyze_left_right_PFLm_diff_and_fwd( imaging_data_1, PFL_FB_dF_F_per_trial, PFL_LAL_dF_F_per_trial, bdata_vel_time, bdata_vel, analysis_path, VPS, sid )

ac = get_analysis_constants;

num_trials = size(imaging_data_1, 1);
nframes = size(imaging_data_1, 5);

t = [0:nframes-1]./VPS;

FFT_FILTER_CUTOFF = 1;

for tr = 1:num_trials
    f = figure('units','normalized','outerposition',[0 0 1 1]);
            
    ax( 1 ) = subplot( 4, 1, 1 );
    cur_PFL_FB = squeeze(PFL_FB_dF_F_per_trial(tr, :,:));        
    imagesc( t, [1:size(PFL_FB_dF_F_per_trial,2)], cur_PFL_FB );
    colormap(flipud(gray));
    caxis([-0.5 4]);
    
    ax( 2 ) = subplot( 4, 1, 2 );
    cur_PFL_LAL = squeeze(PFL_LAL_dF_F_per_trial(tr, :,:));        
%     imagesc( t, [1:size(PFL_LAL_dF_F_per_trial,2)], cur_PFL_LAL );
%     colormap(flipud(gray));
%     caxis([-0.5 2]);
    
    left_PFLm_LAL_axon_end    = squeeze( cur_PFL_LAL(1,:));
    right_PFLm_LAL_axon_end    = squeeze( cur_PFL_LAL(3,:));
    
    left_PFLm_LAL_axon_end_filt     = fft_filter( left_PFLm_LAL_axon_end', FFT_FILTER_CUTOFF, VPS );    
    right_PFLm_LAL_axon_end_filt    = fft_filter( right_PFLm_LAL_axon_end', FFT_FILTER_CUTOFF, VPS );    


    left_PFLm_LAL_axon_end_filt = left_PFLm_LAL_axon_end_filt - mean( left_PFLm_LAL_axon_end_filt );
    right_PFLm_LAL_axon_end_filt = right_PFLm_LAL_axon_end_filt - mean( right_PFLm_LAL_axon_end_filt );

    hold on;
        
    plt1 = plot( t, left_PFLm_LAL_axon_end_filt, 'color', 'r', 'LineWidth', 2 );            
    plt2 = plot( t, right_PFLm_LAL_axon_end_filt, 'g', 'LineWidth', 2 );        
    ylim([-0.8 0.8]);

    ax( 3 ) = subplot( 4, 1, 3 );
    left_right_PB_LAL_sum_end = left_PFLm_LAL_axon_end_filt + right_PFLm_LAL_axon_end_filt;
    left_right_PB_LAL_sum_end = left_right_PB_LAL_sum_end - mean(left_right_PB_LAL_sum_end);
    
    hold on;
    plot( t, left_right_PB_LAL_sum_end,   'LineWidth', 2, 'color', 'b', 'DisplayName', 'axon end' );      
    ylim([-1.5 1.5]);
    legend();
    
    ax( 4 ) = subplot( 4, 1, 4 );
    plot( bdata_vel_time,  squeeze(bdata_vel{ 1 }( tr, ac.VEL_FWD, : )) );
    
    ylabel('Fwd vel');
    xlabel('Time (s)');    
    linkaxes(ax,'x');
        
    legend( [plt1, plt2], {'Left LAL', 'Right LAL'} );
    
    saveas(f,[analysis_path '/fwd_roi_traces_sid_' num2str(sid) '_tid_' num2str(tr) '.fig']);
    saveas(f,[analysis_path '/fwd_roi_traces_sid_' num2str(sid) '_tid_' num2str(tr) '.png']);
    % close(f);          
end
end

