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
% directories_to_analyze = { { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_03', 0, 1900 }, ...                          
%                             { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_04', 1, 2346 }, ...                           
%                             { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_05', 0, 1740 }, ... 
%                             { '161210_lexAOpCsChrimson_gfp_83blexA_ss730_06', 0, 1882 } };

directories_to_analyze = { { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_05', 0, 1740 } };                        
%directories_to_analyze = { { '161218_lexAOpCsChrimson_gfp_83blexA_ss731_02', 0, 3000 } }; 

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
            
    thresh_idx = 1;
    MAX_FR = 150;
    MIN_FR = 5;
    INC_FR = 5;
    
    WINDOW_TOTAL_SIZE = WINDOW_BACK_SR+WINDOW_FRONT_SR+1;
    FR_RATE_INC_SIZE = (MAX_FR-MIN_FR) / INC_FR;
    yaw_image = zeros( FR_RATE_INC_SIZE, WINDOW_TOTAL_SIZE ); 
    fr_image = zeros( FR_RATE_INC_SIZE, WINDOW_TOTAL_SIZE ); 
    latencies = [];
    
    for thresh = MIN_FR:INC_FR:MAX_FR
    
        FR_THRESHOLD = thresh;
    
        locs = find_thresholds_in_PSTH( t_vel_all{d}, psth1, FR_THRESHOLD );
    
%         ff = figure;
%         hold on;
%         plot(t_vel_all{d}, psth1);
%         plot(t_vel_all{d}(locs), psth1(locs), 'x', 'color', 'g');
    
        for i=1:length(locs)

            cur_loc = locs(i);

            beg_offset = cur_loc - WINDOW_BACK_SR;
            end_offset = cur_loc + WINDOW_FRONT_SR;

            if(beg_offset < 1) continue; end
            if(end_offset > length(psth1)) continue; end

            t_plot   = t_vel_all{ d }( beg_offset:end_offset) - t_vel_all{d}( cur_loc );
            cur_psth_plot = psth1( beg_offset:end_offset);

            psth_win_all_data{thresh_idx}(i, :) = cur_psth_plot;

            yaw_plot = cur_yaw( beg_offset: end_offset );
            yaw_win_all_data{thresh_idx}(i, :) = yaw_plot;
        end
        
        psth_avg = squeeze( mean( psth_win_all_data{thresh_idx}, 1 ) );
        fr_image(thresh_idx ,:) = psth_avg;
        psth_sem = get_sem( psth_win_all_data{thresh_idx}, 1);
        ipt_fr = findchangepts(psth_avg, 'Statistic', 'std');
        
        yaw_avg = squeeze( mean( yaw_win_all_data{thresh_idx} ,1 ) );
        yaw_image(thresh_idx ,:) = yaw_avg;
        yaw_sem = get_sem( yaw_win_all_data{thresh_idx}, 1 );   
        ipt_yaw = findchangepts(yaw_avg, 'Statistic', 'std');

        cur_latency = t_plot(ipt_yaw) - t_plot(ipt_fr);
        latencies(thresh_idx) = cur_latency;
        
        %if( thresh == 70 )
        if( 0 )
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
            pp = plot(t_plot(ipt_yaw), yaw_avg(ipt_yaw) );
            pp.Marker = '*';

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
            pp = plot(t_plot(ipt_fr), psth_avg(ipt_fr) );
            pp.Marker = '*';

            xlim([-1*WINDOW_BACK WINDOW_FRONT]);
            
            ylabel( 'Firing rate( spikes/s )');        
        end
        
        thresh_idx = thresh_idx + 1;
    end
    
    f = figure;
    plot([MIN_FR:INC_FR:MAX_FR], latencies );
    
    f = figure;
    
    surf(t_plot, [MIN_FR:INC_FR:MAX_FR], yaw_image );
    title('Yaw vs. threshold');
    xlabel('Time (s)');
    ylabel('Firing rate threshold (spikes/s)');
    zlabel('Triggered yaw (deg/s)');
    set(gca, 'FontSize', 14);
    axis tight;
    
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_yaw_img_thresholds.fig' ]);
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_yaw_img_thresholds.png' ]);
    
    f = figure;
    
    surf(t_plot, [MIN_FR:INC_FR:MAX_FR], fr_image );
    title('FR vs. threshold');
    xlabel('Time (s)');
    ylabel('Firing rate threshold (spikes/s)');
    zlabel('Avg FR');
    set(gca, 'FontSize', 14);
    axis tight;
    
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_fr_img_thresholds.fig' ]);
    saveas(f, [summary_analysis_path '/' data_dirname '_sid_' num2str(sid) '_' save_name '_fr_img_thresholds.png' ]);
end

