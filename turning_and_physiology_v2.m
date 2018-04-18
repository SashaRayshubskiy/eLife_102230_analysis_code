%% Load all data from all the animals to be included in this analysis

clear all;
close all;

working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/descending_neurons/';

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Configuration for this analysis
% path to dir, sid, cutoff_time in seconds (this is the time when the
% usable recording cut out). If this value is -1, then a plot of the ephys
% will appear and the user will decide this value by eye.
% Junk:  % { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_03', 1, 1968 }, ...
%          { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_04', 2, 1500 }, ...                           
directories_to_analyze = { { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_03', 0, 1900 }, ...                          
                            { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_04', 1, 2346 }, ...                           
                            { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_05', 0, 1740 }, ... 
                            { '161210_lexAOpCsChrimson_gfp_83blexA_ss730_06', 0, 1882 } };

%directories_to_analyze = { { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_05', 0, 1740 } };                        
%directories_to_analyze = { { '161218_lexAOpCsChrimson_gfp_83blexA_ss731_02', 0, 3000 } }; 

% directories_to_analyze = { { '161218_lexAOpCsChrimson_gfp_83blexA_ss731_02', 0, 3000 }, ...                          
%                            { '170721_lexAOpCsChrimson_gfp_83blexA_ss731_05', 2, 2829 }, ...                           
%                            { '170726_lexAOpCsChrimson_gfp_83blexA_ss731_08', 0, 1656 }, ...
%                            { '170727_lexAOpCsChrimson_gfp_83blexA_ss731_12', 0, 2800 }};

pre_stim_t = 3.0;
stim_t =  3.5;
post_stim_t = 3.0;
inter_trial_t = 5.0;

FIRING_RATE_THRESHOLD = 40;
YAW_THRESHOLD_LEFT = -600;
YAW_THRESHOLD_RIGHT = 600;

ALIGN_ON_FIRING_RATE = 10;
ALIGN_ON_YAW_LEFT = 11;
ALIGN_ON_YAW_RIGHT = 12;
WINDOW_BACK = 0.6;
WINDOW_FRONT = 0.4;

ANALYSIS_PERIOD_TYPE_ALL         = 20;
ANALYSIS_PERIOD_TYPE_ODOR_EVOKED = 21;
ANALYSIS_PERIOD_TYPE_SPONTANEOUS = 22;

%%%% KEY VARIABLES %%%%
ANALYSIS_PERIOD_TYPE = ANALYSIS_PERIOD_TYPE_ALL;
ANALYSIS_PERIOD_TYPE_TEXT = 'all';

ALIGNMENT_METHOD = ALIGN_ON_FIRING_RATE;
REPORT_THRESHOLD = FIRING_RATE_THRESHOLD;
save_name = 'FR_predicting_turning';

% ALIGNMENT_METHOD = ALIGN_ON_YAW_LEFT;
% REPORT_THRESHOLD = YAW_THRESHOLD_LEFT;
% save_name = 'turning_predicting_FR';

% ALIGNMENT_METHOD = ALIGN_ON_YAW_RIGHT;
% REPORT_THRESHOLD = YAW_THRESHOLD_RIGHT;
% save_name = 'turning_predicting_FR';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timenow_str = datestr(datetime, 'yymmdd_HHMMSS');
summary_analysis_path = [working_dir '/summary_analysis/ephys_vs_yaw_analysis/triggered_analysis_' timenow_str];

if(~exist(summary_analysis_path, 'dir'))
    mkdir(summary_analysis_path);
end

[t_all, t_vel_all, yaw_all, fwd_all, ephys_all_A, ephys_all_B] = load_LAL_DN_data( working_dir, directories_to_analyze, ephys_SR, ball_SR );

EVENT_FIND_PEAKS = 10;
EVENT_FR_THRESHOLD_CROSSING = 11;

EVENT_FINDING_METHOD = EVENT_FR_THRESHOLD_CROSSING;

%% 
% 1. Does firing rate predict a subsequent turn? 
% 2. Does the turning predict the firing rate?

WINDOW_BACK_SR = WINDOW_BACK * ball_SR;
WINDOW_FRONT_SR = WINDOW_FRONT * ball_SR;

WINDOW_BACK_EPHYS_SR = WINDOW_BACK * ephys_SR;
WINDOW_FRONT_EPHYS_SR = WINDOW_FRONT * ephys_SR;

psth_win_all_data = {};
yaw_win_all_data = {};

yaw_mins = {};
psth_max = {};

dataset_cnt = length(directories_to_analyze);

SPIKE_THRESHOLD = 0.4;
BIN_SIZE = ephys_SR/ball_SR;

for d = 1:dataset_cnt
    
    data_dirname = directories_to_analyze{ d }{ 1 };
    sid = directories_to_analyze{ d }{ 2 };
    
    psth1 =  calculate_psth( t_all{d}, t_vel_all{d}, ephys_all_A{d}, ephys_SR, SPIKE_THRESHOLD, BIN_SIZE);
    
    cur_yaw = yaw_all{d};
        
    if( ANALYSIS_PERIOD_TYPE == ANALYSIS_PERIOD_TYPE_ODOR_EVOKED )
        
    elseif( ANALYSIS_PERIOD_TYPE == ANALYSIS_PERIOD_TYPE_SPONTANEOUS )
    
    end
    
    if 0
        f = figure;        
        [r, lags] = xcorr( psth1 );        
        plot(lags,r, 'DisplayName', 'xcorr on PSTH');  
        legend('show');
        saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_xcorr_on_psdh_' ANALYSIS_PERIOD_TYPE_TEXT '.fig' ]);
        saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_xcorr_on_psdh_' ANALYSIS_PERIOD_TYPE_TEXT '.png' ]);
        close(f);

        f = figure;        
        [r, lags] = xcorr( cur_yaw );        
        plot(lags,r, 'DisplayName', 'xcorr on yaw');        
        legend('show');
        saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_xcorr_on_yaw_' ANALYSIS_PERIOD_TYPE_TEXT '.fig' ]);
        saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_xcorr_on_yaw_' ANALYSIS_PERIOD_TYPE_TEXT '.png' ]);
        close(f);
    end
    
    if ( d == 3 )
        f = figure;
        hold on;
        yyaxis left;
        hold on;
        %plot( t_vel_all{d}, psth1 );
        plot( t_all{d}, ephys_all_A{d} );
        ylabel( 'Voltage (mV)' );
        
        yyaxis right;
        hold on;
        plot( t_vel_all{d}, cur_yaw );
        ylabel( 'Yaw velocity ( deg/s )' );
        xlabel('Time (s)');
        title( ['Firing rate vs. yaw velocity' ] );
        %close(f);        
    end
    
    if 0
        f = figure;
        hold on;
        yyaxis left;
        hold on;
        %plot( t_vel_all{d}, psth1 );
        plot( t_vel_all{d}, psth1 );
        ylabel( 'Firing rate ( spikes/s )' );
        
        yyaxis right;
        hold on;
        plot( t_vel_all{d}, cur_yaw );
        ylabel( 'Yaw velocity ( deg/s )' );
        xlabel('Time (s)');
        title( ['Firing rate vs. yaw velocity' ] );
        saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_FR_and_yaw_velocity_' ANALYSIS_PERIOD_TYPE_TEXT '.fig' ]);
        saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_FR_and_yaw_velocity_' ANALYSIS_PERIOD_TYPE_TEXT '.png' ]);
        close(f);        
    end
    
    if( ALIGNMENT_METHOD == ALIGN_ON_FIRING_RATE )
        
        if ( EVENT_FINDING_METHOD == EVENT_FIND_PEAKS )
            psth1(find(psth1 < FIRING_RATE_THRESHOLD) ) = 0;
            
            % [~, locs] = findpeaks(psth1, 'Threshold', 20, 'MinPeakProminence', 60, 'MinPeakDistance', 20, 'Annotate','extents');
            [~, locs] = findpeaks(psth1, 'MinPeakDistance', 50, 'Annotate','extents');
            
            ff = figure;
            
            hold on;
            plot(t_vel_all{d}, psth1);
            plot(t_vel_all{d}(locs), psth1(locs), 'x', 'color', 'g');
            %         plot(t_all{d}, psth1);
            %         plot(t_all{d}(locs), psth1(locs), 'x', 'color', 'g');
            
            % waitforbuttonpress;
            %close(ff);
            % break;
        elseif( EVENT_FINDING_METHOD == EVENT_FR_THRESHOLD_CROSSING )
            
            FR_THRESHOLD = 80;
            
            locs = find_thresholds_in_PSTH( t_vel_all{d}, psth1, FR_THRESHOLD );
            
            ff = figure;
            hold on;
            plot(t_vel_all{d}, psth1);
            plot(t_vel_all{d}(locs), psth1(locs), 'x', 'color', 'g');            
        end               
        
    elseif( ALIGNMENT_METHOD == ALIGN_ON_YAW_LEFT )
        cur_yaw(find(cur_yaw > YAW_THRESHOLD_LEFT) ) = 0; % for LEFT turns
        [~, locs] = findpeaks(-1.0*cur_yaw, 'Annotate','extents');
    elseif( ALIGNMENT_METHOD == ALIGN_ON_YAW_RIGHT )
        cur_yaw(find(cur_yaw < YAW_THRESHOLD_RIGHT) ) = 0; % for RIGHT turns
        [~, locs] = findpeaks(cur_yaw, 'Annotate','extents');
    end
        
    %f = figure;
    latencies = [];
    for i=1:length(locs)        
        
        cur_loc = locs(i);
%        cur_loc = ceil(locs(i)/40.0);
        
        %subplot(2,1,1)
        %hold on;
        
        beg_offset = cur_loc - WINDOW_BACK_SR;
        end_offset = cur_loc + WINDOW_FRONT_SR;
        
%         beg_offset_e = cur_loc_e - WINDOW_BACK_EPHYS_SR;
%         end_offset_e = cur_loc_e + WINDOW_FRONT_EPHYS_SR;
                
        if(beg_offset < 1) continue; end
        if(end_offset > length(psth1)) continue; end
        
        t_plot   = t_vel_all{ d }( beg_offset:end_offset) - t_vel_all{d}( cur_loc );
        % t_plot_e = t_all{ d }( beg_offset_e:end_offset_e ) - t_all{d}( cur_loc_e );
        cur_psth_plot = psth1( beg_offset:end_offset);

        psth_win_all_data{d}(i, :) = cur_psth_plot;
        
        %plot( t_plot, cur_psth_plot, 'color', rgb('PaleGreen') );
        %ylabel('Firing rate (spikes/s)');
        %xlabel('Time (s)');
        %xlim([-1*WINDOW_BACK WINDOW_FRONT]);
        
        %subplot(2,1,2)
        %hold on;
        
        yaw_plot = cur_yaw( beg_offset: end_offset );
        
        yaw_win_all_data{d}(i, :) = yaw_plot;
        
%         plot( t_plot, yaw_plot, 'color', rgb('Violet') );
%         ylabel('Yaw velocity (deg/s)');
%         xlabel('Time (s)');
%         xlim([-1*WINDOW_BACK WINDOW_FRONT]);

        % [~, locs1] = findpeaks(cur_psth_plot, 'NPeaks', 1, 'SortStr', 'descend', 'Annotate','extents');   
        
        %[~, locs0] = find( t_plot > 0 );        
        
        [~, locs2] = findpeaks(-1*yaw_plot(find( t_plot > 0 )), 'NPeaks', 1, 'SortStr', 'descend', 'Annotate','extents');
        
        
        %if( (length(locs1) > 0) & (length(locs2) > 0) )      
        %if( (length(locs2) > 0) )
        if( 0 )

            yaw_offset = locs2(1) + WINDOW_BACK_SR;
            
            % FR_peak_time = t_plot(locs1(1));
            FR_peak_time = 0;
            yaw_min_time = t_plot( yaw_offset );
            cur_latency = yaw_min_time - FR_peak_time;
            latencies(end+1) = cur_latency;    
            
            if 1
            fff = figure;

            yyaxis left;
            hold on;
            plot(t_plot, cur_psth_plot);
            plot(t_plot(locs1(1)), cur_psth_plot(WINDOW_BACK_SR), 'x');
            
            yyaxis right;
            hold on;
            plot(t_plot, yaw_plot);
            plot(t_plot(yaw_offset), yaw_plot(yaw_offset), 'x');

            waitforbuttonpress;
            close(fff);
            end
        end
    end
    
    f = figure;
    histogram(latencies, 20);
    xlabel('Time (s)');
    ylabel('Turn event count');
    set(gca, 'FontSize', 16);
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_turn_peak_latency_hist_' ANALYSIS_PERIOD_TYPE_TEXT '.fig' ]);
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_turn_peak_latency_hist_' ANALYSIS_PERIOD_TYPE_TEXT '.png' ]);
    close(f);

    
%     subplot(2,1,1);
%     hold on;
    psth_avg = squeeze( mean( psth_win_all_data{d}, 1 ) );
    psth_sem = get_sem( psth_win_all_data{d}, 1);
    
%     plot(t_plot, psth_avg, 'color', rgb('SeaGreen'), 'LineWidth', 2.0);
%     xlim([-1*WINDOW_BACK WINDOW_FRONT]);
    
%     subplot(2,1,2);
%     hold on;
    yaw_avg = squeeze( mean( yaw_win_all_data{d} ,1 ) );
    yaw_sem = get_sem( yaw_win_all_data{d}, 1 );
    
%     plot(t_plot, yaw_avg, 'color', rgb('DarkMagenta'), 'LineWidth', 2.0);
%     xlim([-1*WINDOW_BACK WINDOW_FRONT]);

%     if( ALIGNMENT_METHOD == ALIGN_ON_FIRING_RATE )
%         title(['Firing rate peak triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(psth_win_all_data{d},1)) ]);
%     else
%         title(['Turning triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(psth_win_all_data{d},1)) ]);        
%     end
%     
%     saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_all_thresh_' num2str(REPORT_THRESHOLD) '.fig' ]);
%     saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_all_thresh_' num2str(REPORT_THRESHOLD) '.png' ]);
%     close(f);
    
    f = figure;
    
    left_color = rgb('SeaGreen');
    right_color = rgb('DarkMagenta');
    set(f,'defaultAxesColorOrder',[right_color; left_color]);

    yyaxis left;
    hold on;

    fh = fill( [t_plot', fliplr(t_plot')], ...
        [(yaw_avg-yaw_sem) fliplr((yaw_avg+yaw_sem))], ...
        rgb('Violet'));
    set(fh, 'EdgeColor', 'None');
    
    plot(t_plot, yaw_avg, 'color', rgb('DarkMagenta'), 'LineStyle', '-', 'LineWidth', 2.0 );
    xlim([-1*WINDOW_BACK WINDOW_FRONT]);    
    
    ylabel( 'Yaw velocity ( deg/s )');
    xlabel('Time (s)');
    
    yyaxis right;
    hold on;
        
    fh = fill( [t_plot' fliplr(t_plot')], ...
        [(psth_avg-psth_sem) fliplr(psth_avg+psth_sem)], ...
        rgb('PaleGreen'));
    set(fh, 'EdgeColor', 'None');

    plot(t_plot, psth_avg, 'color', rgb('SeaGreen'), 'LineStyle', '-', 'LineWidth', 2.0 );
    xlim([-1*WINDOW_BACK WINDOW_FRONT]);

    ylabel( 'Firing rate( spikes/s )');
    
    if( ALIGNMENT_METHOD == ALIGN_ON_FIRING_RATE )
        title(['Firing rate peak triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(psth_win_all_data{d},1)) ]);
    else
        title(['Turning triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(psth_win_all_data{d},1)) ]);
    end

    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_avg_thresh_' num2str(REPORT_THRESHOLD) '.fig' ]);
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_avg_thresh_' num2str(REPORT_THRESHOLD) '.png' ]);
    close(f);
    
    if 0
    % Plot a scatter plot of FR peak vs. turning magnitude, right after the
    % peak
    PEAK_WINDOW_AHEAD = 0.3 * ball_SR;
    PEAK_WINDOW_BACK = 0.3 * ball_SR;
    for i=1:length(locs)
        
        cur_loc = locs(i);
        
        beg_offset = cur_loc - PEAK_WINDOW_BACK;
        end_offset = cur_loc + PEAK_WINDOW_AHEAD;
        if( beg_offset < 1 ) continue; end
        if( end_offset > length(psth1)) continue; end
        
        if ( ALIGNMENT_METHOD == ALIGN_ON_FIRING_RATE )
            psth_max{d}(i) = psth1( cur_loc );
            yaw_mins{d}(i) = min( cur_yaw( beg_offset : end_offset ));
        else
            psth_max{d}(i) = max( psth1( beg_offset : end_offset ));
            yaw_mins{d}(i) = cur_yaw(cur_loc);        
        end
    end
    
    f = figure;
    p = polyfit(psth_max{d}, yaw_mins{d}, 1 );
    
    xx = linspace( 0, 140, 14000 );
    
    yy = polyval(p,xx);
    hold on;
    scatter( psth_max{d}, yaw_mins{d} );
    hold on;
    plot(xx,yy);
    
    yresid = yaw_mins{d} - polyval(p,psth_max{d});
    
    SSresid = sum(yresid.^2);
    SStotal = (length(yaw_mins{d})-1) * var(yaw_mins{d});
    rsq = 1 - SSresid/SStotal;
    
    xlabel('Firing rate (spikes/s)');
    ylabel('Yaw velocity (deg/s)');
    title(['Firing rate peak vs. yaw magnitude threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(psth_win_all_data{d},1)) ' R^2: ' num2str(rsq) ]);
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_FR_peak_vs_turning_mag_thresh_' num2str(REPORT_THRESHOLD) '.fig' ]);
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_FR_peak_vs_turning_mag_thresh_' num2str(REPORT_THRESHOLD) '.png' ]);
    close(f);
    end
end

%%
% Plot data for all fly datasets specied. This code needs the above to run first to generate all the data. 

% FR Peak triggered turning magnitude, all trials.
psth_win_all_data_flat = [];
yaw_win_all_data_flat = [];
for d = 1:length(psth_win_all_data)
    
    cur_psths = psth_win_all_data{d};
    cur_yaws = yaw_win_all_data{d};
    
    psth_win_all_data_flat = vertcat( psth_win_all_data_flat, cur_psths );
    yaw_win_all_data_flat = vertcat( yaw_win_all_data_flat, cur_yaws );    
end

% f = figure;
% for i=1:size(psth_win_all_data_flat,1)
%     
%     cur_psth_plot = psth_win_all_data_flat(i,:);
%     cur_yaw_plot = yaw_win_all_data_flat(i,:);    
%     
%     subplot(2,1,1)
%     hold on;
%             
%     plot( t_plot, cur_psth_plot, 'color', rgb('PaleGreen') );
%     ylabel('Firing rate (spikes/s)');
%     xlabel('Time (s)');
%     xlim([-1*WINDOW_BACK WINDOW_FRONT]);
%     
%     subplot(2,1,2)
%     hold on;
%         
%     plot( t_plot, cur_yaw_plot, 'color', rgb('Violet') );
%     ylabel('Yaw velocity (deg/s)');
%     xlabel('Time (s)');
%     xlim([-1*WINDOW_BACK WINDOW_FRONT]);
% end
% 
% subplot(2,1,1);
% hold on;

psth_avg = squeeze( mean( psth_win_all_data_flat, 1 ) );
psth_sem = get_sem( psth_win_all_data_flat, 1);

% plot(t_plot, psth_avg, 'color', rgb('SeaGreen'), 'LineWidth', 2.0);
% 
% subplot(2,1,2);
% hold on;
yaw_avg = squeeze( mean( yaw_win_all_data_flat ,1 ) );
yaw_sem = get_sem( yaw_win_all_data_flat, 1 );

% plot(t_plot, yaw_avg, 'color', rgb('DarkMagenta'), 'LineWidth', 2.0);
% if( ALIGNMENT_METHOD == ALIGN_ON_FIRING_RATE )
%     title(['Firing rate peak triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(yaw_win_all_data_flat,1)) ]);
% else
%     title(['Turning triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(yaw_win_all_data_flat,1)) ]);
% end

timenow_str = datestr(datetime, 'yymmdd_HHMMSS');
save([summary_analysis_path '/' timenow_str '_directories_to_analyze.mat'],'directories_to_analyze');
% saveas(f, [summary_analysis_path '/' timenow_str '_' save_name '_all_thresh_' num2str(REPORT_THRESHOLD) '.fig' ]);
% saveas(f, [summary_analysis_path '/' timenow_str '_' save_name '_all_thresh_' num2str(REPORT_THRESHOLD) '.png' ]);
% close(f);

% FR Peak triggered turning magnitude, avg.
f = figure;

left_color = rgb('SeaGreen');
right_color = rgb('DarkMagenta');
set(f,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left;
hold on;

fh = fill( [t_plot', fliplr(t_plot')], ...
    [(psth_avg-psth_sem) fliplr((psth_avg+psth_sem))], ...
    rgb('PaleGreen'));
set(fh, 'EdgeColor', 'None');

plot(t_plot, psth_avg, 'color', rgb('SeaGreen'), 'LineStyle', '-', 'LineWidth', 2.0 );

ylabel( 'Firing rate( spikes/s )');

yyaxis right;
hold on;

fh = fill( [t_plot', fliplr(t_plot')], ...
    [(yaw_avg-yaw_sem) fliplr((yaw_avg+yaw_sem))], ...
    rgb('Violet'));
set(fh, 'EdgeColor', 'None');

plot(t_plot, yaw_avg, 'color', rgb('DarkMagenta'), 'LineStyle', '-', 'LineWidth', 2.0 );

ylabel( 'Yaw velocity ( deg/s )');
xlabel('Time (s)');
if( ALIGNMENT_METHOD == ALIGN_ON_FIRING_RATE )
    title(['Firing rate peak triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(yaw_win_all_data_flat,1)) ]);
else
    title(['Turning triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(yaw_win_all_data_flat,1)) ]);
end

saveas(f, [summary_analysis_path '/' timenow_str '_' save_name '_avg_thresh_' num2str(REPORT_THRESHOLD) '.fig' ]);
saveas(f, [summary_analysis_path '/' timenow_str '_' save_name '_avg_thresh_' num2str(REPORT_THRESHOLD) '.png' ]);
close(f);

% Scatter plot of FP peak vs. turning mag
psth_max_all = [];
yaw_mins_all = [];

for d=1:length(psth_max)
    psth_max_all = horzcat(psth_max_all, psth_max{d});
    yaw_mins_all = horzcat(yaw_mins_all, yaw_mins{d});
end

f = figure;
p = polyfit(psth_max_all, yaw_mins_all, 1 );

xx = linspace( 0, 140, 14000 );

yy = polyval(p,xx);
hold on;
scatter( psth_max_all, yaw_mins_all );
hold on;
plot(xx,yy);

yresid = yaw_mins_all - polyval(p,psth_max_all);

SSresid = sum(yresid.^2);
SStotal = (length(yaw_mins_all)-1) * var(yaw_mins_all);
rsq = 1 - SSresid/SStotal;

xlabel('Firing rate (spikes/s)');
ylabel('Yaw velocity (deg/s)');
title(['Firing rate peak vs. yaw magnitude threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(length(yaw_mins_all)) ' R^2: ' num2str(rsq) ]);
saveas(f, [summary_analysis_path '/' timenow_str '_FR_peak_vs_turning_mag_thresh_' num2str(REPORT_THRESHOLD) '.fig' ]);
saveas(f, [summary_analysis_path '/' timenow_str '_FR_peak_vs_turning_mag_thresh_' num2str(REPORT_THRESHOLD) '.png' ]);
close(f);

%% Future questions:
% 3. Compare odor evoked and spontaneous firing rates and turns.
% 3. Compare derivative of FR with turning velocity, derivative of turning velocity.
% 4. Do contraleteral turns correspond to a decrease (or firing rate of
% zero)?

