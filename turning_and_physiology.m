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
%directories_to_analyze = { { '161218_lexAOpCsChrimson_gfp_83blexA_ss731_02', 0, 3000 } }; 

pre_stim_t = 3.0;
stim_t =  3.5;
post_stim_t = 3.0;
inter_trial_t = 5.0;

FIRING_RATE_THRESHOLD = 40;
YAW_THRESHOLD_LEFT = -600;
YAW_THRESHOLD_RIGHT = 600;

SPIKE_THRESHOLD = 0.4;
BIN_SIZE = ephys_SR/ball_SR;

ALIGN_ON_FIRING_RATE = 10;
ALIGN_ON_YAW_LEFT = 11;
ALIGN_ON_YAW_RIGHT = 12;
WINDOW_BACK = 1;
WINDOW_FRONT = 1;

ANALYSIS_PERIOD_TYPE_ALL         = 20;
ANALYSIS_PERIOD_TYPE_ODOR_EVOKED = 21;
ANALYSIS_PERIOD_TYPE_SPONTANEOUS = 22;

%%%% KEY VARIABLES %%%%
ANALYSIS_PERIOD_TYPE = ANALYSIS_PERIOD_TYPE_ALL;
ANALYSIS_PERIOD_TYPE_TEXT = 'all';

%ALIGNMENT_METHOD = ALIGN_ON_FIRING_RATE;
%REPORT_THRESHOLD = FIRING_RATE_THRESHOLD;
%save_name = 'FR_predicting_turning';

% ALIGNMENT_METHOD = ALIGN_ON_YAW_LEFT;
% REPORT_THRESHOLD = YAW_THRESHOLD_LEFT;
% save_name = 'turning_predicting_FR';

ALIGNMENT_METHOD = ALIGN_ON_YAW_RIGHT;
REPORT_THRESHOLD = YAW_THRESHOLD_RIGHT;
save_name = 'turning_predicting_FR';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
summary_analysis_path = [working_dir '/summary_analysis/ephys_vs_yaw_analysis/'];

if(~exist(summary_analysis_path, 'dir'))
    mkdir(summary_analysis_path);
end

dataset_cnt = length(directories_to_analyze);

t_all = cell(1,dataset_cnt);
t_vel_all = cell(1,dataset_cnt);
yaw_all = cell(1,dataset_cnt);
fwd_all = cell(1,dataset_cnt);
ephys_all_A = cell(1,dataset_cnt);
ephys_all_B = cell(1,dataset_cnt);

for d = 1:dataset_cnt
    
    datapath = [ working_dir directories_to_analyze{ d }{ 1 } ];
    
    sid = [ directories_to_analyze{ d }{ 2 }];
    
    disp(['Loading directory: ' datapath '   sid:' num2str(sid) ]);
    
    TIME_CUTOFF_SECONDS = directories_to_analyze{ d }{ 3 }; % Analysis for: 161208_lexAOpCsChrimson_gfp_83blexA_ss730_05, sid=0
    
    dd = {};
    for s = 1:length(sid)
        
        dd_tmp = dir([datapath '/*_sid_' num2str(sid(s)) '_*.mat']);
        
        % Make sure that the data is in sorted order
        num_trials = length(dd_tmp);
        sort_vector1 = zeros(num_trials, 3);
        
        for i=1:num_trials
            cur_name = dd_tmp(i).name;
            
            cur_name_token = strsplit(cur_name, '_');
            sid = str2double(cur_name_token{6});
            tid_str = strsplit(cur_name_token{8}, '.');
            
            tid = str2double(tid_str{1});
            
            sort_vector1(i, 1) = sid;
            sort_vector1(i, 2) = tid;
            sort_vector1(i, 3) = i;
        end
        
        sort_result_1 = sortrows(sort_vector1);
        
        dd_sorted = cell(1,num_trials);
        
        for i=1:num_trials
            cur_idx = sort_result_1(i,3);
            dd_sorted{i} = dd_tmp(cur_idx).name;
        end
        
        dd = vertcat(dd, dd_sorted);
    end       
    
    tmp_t_all = [];
    tmp_t_vel_all = [];
    tmp_yaw_all = [];
    tmp_fwd_all = [];
    tmp_ephys_all_A = [];
    tmp_ephys_all_B = [];
    
    % Load in all the DAQ data
    for i=1:length(dd)
        filepath = [datapath '/' dd{i}];
        
        cur_data = load(filepath);
        
        D = cur_data.trial_bdata;
        
        [currentA, voltageA, currentB, voltageB] = get_dual_scaled_voltage_and_current( D );
        
        t = cur_data.trial_time;
        
        [ t_vel, vel_forward, vel_side, vel_yaw ] = get_velocity_ephysrig(t, D, datapath, 1);
        
        tmp_t_all = vertcat( tmp_t_all, t );
        tmp_t_vel_all = vertcat( tmp_t_vel_all, t_vel' );
        tmp_yaw_all = vertcat( tmp_yaw_all, vel_yaw' );
        tmp_fwd_all = vertcat( tmp_fwd_all, vel_forward' );
        tmp_ephys_all_A = vertcat( tmp_ephys_all_A, voltageA );
        tmp_ephys_all_B = vertcat( tmp_ephys_all_B, voltageB );
    end
        
    if( TIME_CUTOFF_SECONDS == -1 )
        figure;
        subplot(2,1,1);
        plot(tmp_t_all, tmp_ephys_all_A);
        xlabel('Time (s)');
        ylabel('Voltage (mV)');
        title('Patch A');

        subplot(2,1,2);
        plot(tmp_t_all, tmp_ephys_all_B);
        xlabel('Time (s)');
        ylabel('Voltage (mV)');
        title('Patch B');
        waitforbuttonpress;
    else        
        TIME_CUTOFF = TIME_CUTOFF_SECONDS * ephys_SR;
        TIME_CUTOFF_VEL = TIME_CUTOFF_SECONDS * ball_SR;
        
        t_all{d}       = tmp_t_all(1:TIME_CUTOFF);
        ephys_all_A{d} = tmp_ephys_all_A(1:TIME_CUTOFF);
        ephys_all_B{d} = tmp_ephys_all_B(1:TIME_CUTOFF);
        
        yaw_all{d}     = tmp_yaw_all(1:TIME_CUTOFF_VEL);
        fwd_all{d}     = tmp_fwd_all(1:TIME_CUTOFF_VEL);
        t_vel_all{d}   = tmp_t_vel_all(1:TIME_CUTOFF_VEL);
    end
end

%% 
% 1. Does firing rate predict a subsequent turn? 
% 2. Does the turning predict the firing rate?

WINDOW_BACK_SR = WINDOW_BACK * ball_SR;
WINDOW_FRONT_SR = WINDOW_FRONT * ball_SR;

psth_win_all_data = {};
yaw_win_all_data = {};

yaw_mins = {};
psth_max = {};

for d = 1:dataset_cnt
    
    data_dirname = directories_to_analyze{ d }{ 1 };
    sid = directories_to_analyze{ d }{ 2 };
    
    psth1 =  calculate_psth( t_all{d}, t_vel_all{d}, ephys_all_A{d}, ephys_SR, SPIKE_THRESHOLD, BIN_SIZE);
    
    cur_yaw = yaw_all{d};
        
    if( ANALYSIS_PERIOD_TYPE == ANALYSIS_PERIOD_TYPE_ODOR_EVOKED )
        
    elseif( ANALYSIS_PERIOD_TYPE == ANALYSIS_PERIOD_TYPE_SPONTANEOUS )
    
    end
    
    f = figure;
    hold on;
    yyaxis left; 
    hold on;
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
    
    if( ALIGNMENT_METHOD == ALIGN_ON_FIRING_RATE )
        psth1(find(psth1 < FIRING_RATE_THRESHOLD) ) = 0;
        [~, locs] = findpeaks(psth1, 'Annotate','extents');
    elseif( ALIGNMENT_METHOD == ALIGN_ON_YAW_LEFT )
        cur_yaw(find(cur_yaw > YAW_THRESHOLD_LEFT) ) = 0; % for LEFT turns
        [~, locs] = findpeaks(-1.0*cur_yaw, 'Annotate','extents');
    elseif( ALIGNMENT_METHOD == ALIGN_ON_YAW_RIGHT )
        cur_yaw(find(cur_yaw < YAW_THRESHOLD_RIGHT) ) = 0; % for RIGHT turns
        [~, locs] = findpeaks(cur_yaw, 'Annotate','extents');
    end
        
    f = figure;
    for i=1:length(locs)
        
        cur_loc = locs(i);
        
        subplot(2,1,1)
        hold on;
        
        beg_offset = cur_loc - WINDOW_BACK_SR;
        end_offset = cur_loc + WINDOW_FRONT_SR;
        
        if(beg_offset < 1) continue; end
        if(end_offset > length(psth1)) continue; end
        
        t_plot = t_vel_all{d}( beg_offset:end_offset) - t_vel_all{d}(cur_loc);
        cur_psth_plot = psth1( beg_offset:end_offset);

        psth_win_all_data{d}(i, :) = cur_psth_plot;
        
        plot( t_plot, cur_psth_plot, 'color', rgb('PaleGreen') );
        ylabel('Firing rate (spikes/s)');
        xlabel('Time (s)');
        xlim([-1*WINDOW_BACK WINDOW_FRONT]);
        
        subplot(2,1,2)
        hold on;
        
        yaw_plot = cur_yaw( beg_offset: end_offset );
        
        yaw_win_all_data{d}(i, :) = yaw_plot;
        
        plot( t_plot, yaw_plot, 'color', rgb('Violet') );
        ylabel('Yaw velocity (deg/s)');
        xlabel('Time (s)');
        xlim([-1*WINDOW_BACK WINDOW_FRONT]);
    end
    
    subplot(2,1,1);
    hold on;
    psth_avg = squeeze( mean( psth_win_all_data{d}, 1 ) );
    psth_sem = get_sem( psth_win_all_data{d}, 1);
    
    plot(t_plot, psth_avg, 'color', rgb('SeaGreen'), 'LineWidth', 2.0);
    xlim([-1*WINDOW_BACK WINDOW_FRONT]);
    
    subplot(2,1,2);
    hold on;
    yaw_avg = squeeze( mean( yaw_win_all_data{d} ,1 ) );
    yaw_sem = get_sem( yaw_win_all_data{d}, 1 );
    
    plot(t_plot, yaw_avg, 'color', rgb('DarkMagenta'), 'LineWidth', 2.0);
    xlim([-1*WINDOW_BACK WINDOW_FRONT]);

    if( ALIGNMENT_METHOD == ALIGN_ON_FIRING_RATE )
        title(['Firing rate peak triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(psth_win_all_data{d},1)) ]);
    else
        title(['Turning triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(psth_win_all_data{d},1)) ]);        
    end
    
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_all_thresh_' num2str(REPORT_THRESHOLD) '.fig' ]);
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_all_thresh_' num2str(REPORT_THRESHOLD) '.png' ]);
    close(f);
    
    f = figure;
    
    left_color = rgb('SeaGreen');
    right_color = rgb('DarkMagenta');
    set(f,'defaultAxesColorOrder',[left_color; right_color]);

    yyaxis left;
    hold on;
        
    fh = fill( [t_plot' fliplr(t_plot')], ...
        [(psth_avg-psth_sem) fliplr(psth_avg+psth_sem)], ...
        rgb('PaleGreen'));
    set(fh, 'EdgeColor', 'None');

    plot(t_plot, psth_avg, 'color', rgb('SeaGreen'), 'LineStyle', '-', 'LineWidth', 2.0 );
    xlim([-1*WINDOW_BACK WINDOW_FRONT]);

    ylabel( 'Firing rate( spikes/s )');

    yyaxis right;
    hold on;

    fh = fill( [t_plot', fliplr(t_plot')], ...
        [(yaw_avg-yaw_sem) fliplr((yaw_avg+yaw_sem))], ...
        rgb('Violet'));
    set(fh, 'EdgeColor', 'None');
    
    plot(t_plot, yaw_avg, 'color', rgb('DarkMagenta'), 'LineStyle', '-', 'LineWidth', 2.0 );
    xlim([-1*WINDOW_BACK WINDOW_FRONT]);

    ylabel( 'Yaw velocity ( deg/s )');
    xlabel('Time (s)');
    
    if( ALIGNMENT_METHOD == ALIGN_ON_FIRING_RATE )
        title(['Firing rate peak triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(psth_win_all_data{d},1)) ]);
    else
        title(['Turning triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(psth_win_all_data{d},1)) ]);
    end

    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_avg_thresh_' num2str(REPORT_THRESHOLD) '.fig' ]);
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_avg_thresh_' num2str(REPORT_THRESHOLD) '.png' ]);
    close(f);
    
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

f = figure;
for i=1:size(psth_win_all_data_flat,1)
    
    cur_psth_plot = psth_win_all_data_flat(i,:);
    cur_yaw_plot = yaw_win_all_data_flat(i,:);    
    
    subplot(2,1,1)
    hold on;
            
    plot( t_plot, cur_psth_plot, 'color', rgb('PaleGreen') );
    ylabel('Firing rate (spikes/s)');
    xlabel('Time (s)');
    xlim([-1*WINDOW_BACK WINDOW_FRONT]);
    
    subplot(2,1,2)
    hold on;
        
    plot( t_plot, cur_yaw_plot, 'color', rgb('Violet') );
    ylabel('Yaw velocity (deg/s)');
    xlabel('Time (s)');
    xlim([-1*WINDOW_BACK WINDOW_FRONT]);
end

subplot(2,1,1);
hold on;
psth_avg = squeeze( mean( psth_win_all_data_flat, 1 ) );
psth_sem = get_sem( psth_win_all_data_flat, 1);

plot(t_plot, psth_avg, 'color', rgb('SeaGreen'), 'LineWidth', 2.0);

subplot(2,1,2);
hold on;
yaw_avg = squeeze( mean( yaw_win_all_data_flat ,1 ) );
yaw_sem = get_sem( yaw_win_all_data_flat, 1 );

plot(t_plot, yaw_avg, 'color', rgb('DarkMagenta'), 'LineWidth', 2.0);
if( ALIGNMENT_METHOD == ALIGN_ON_FIRING_RATE )
    title(['Firing rate peak triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(yaw_win_all_data_flat,1)) ]);
else
    title(['Turning triggered average threshold: ' num2str(REPORT_THRESHOLD) ' n:' num2str(size(yaw_win_all_data_flat,1)) ]);
end

timenow_str = datestr(datetime, 'yymmdd_HHMMSS');
save([summary_analysis_path '/' timenow_str '_directories_to_analyze.mat'],'directories_to_analyze');
saveas(f, [summary_analysis_path '/' timenow_str '_' save_name '_all_thresh_' num2str(REPORT_THRESHOLD) '.fig' ]);
saveas(f, [summary_analysis_path '/' timenow_str '_' save_name '_all_thresh_' num2str(REPORT_THRESHOLD) '.png' ]);
close(f);

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

%% ATTIC
%%

all_yaw_splined = spline(t_vel_all, yaw_all, t_all);
all_yaw_splined = all_yaw_splined(1:TIME_CUTOFF);

figure;
plotyy( t_all, ephys_all_A, t_all, all_yaw_splined );

