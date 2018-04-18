%% Ephys vs. yaw scatter plots
% Load data
clear all;

global slash;
slash = '/';

working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/descending_neurons/';

% A2
analysis_path_type = 'A2';
directories_to_analyze = { { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_03', 0, 1900, 0 }, ...                          
                           { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_04', 1, 2346, 0 }, ...                           
                           { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_05', 0, 1740, 0 }, ... 
                           { '161210_lexAOpCsChrimson_gfp_83blexA_ss730_06', 0, 1882, 0 } };

% A1
% analysis_path_type = 'A1';
% directories_to_analyze = { { '161218_lexAOpCsChrimson_gfp_83blexA_ss731_02', 0, 3000, 0 }, ...                          
%                            { '170721_lexAOpCsChrimson_gfp_83blexA_ss731_05', 2, 2829, 0 }, ...                           
%                            { '170726_lexAOpCsChrimson_gfp_83blexA_ss731_08', 0, 1656, 5.0 }, ...
%                            { '170727_lexAOpCsChrimson_gfp_83blexA_ss731_12', 0, 2800, 0 }};

analysis_path_base = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/descending_neurons/summary_analysis/';

analysis_path = [analysis_path_base '/' analysis_path_type '_ephys_vs_yaw_analysis/'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

trial_type_cnt = 3;

% Transform raw data into velocity
t_volts = {};
volts_all = {};

t_vel_all = {};
yaw_vel_all = {};
forward_vel_all = {};  

trial_idx_cntrs = ones(trial_type_cnt,1);

dataset_cnt = length(directories_to_analyze);

for d = 1:dataset_cnt
    
    datapath = [ working_dir directories_to_analyze{ d }{ 1 } ];
    
    sid = [ directories_to_analyze{ d }{ 2 }];
    
    disp(['Loading directory: ' datapath '   sid:' num2str(sid) ]);
    
    TIME_CUTOFF_SECONDS = directories_to_analyze{ d }{ 3 };
    FWD_VEL_OFFSET = directories_to_analyze{ d }{ 4 };
    
    cur_TRIAL_CNT_MAX = TIME_CUTOFF_SECONDS / 11.5;
    
    [ bdata_raw, bdata_time, trial_metadata ] = load_behavioral_data( sid, datapath, trial_type_cnt );

    TIME_OFFSET = bdata_time(1);
    
    for tt = 1:trial_type_cnt
        for trial = 1:size(bdata_raw{tt},1)
            
            tid = trial_metadata{tt}(trial,2);
            if(tid > cur_TRIAL_CNT_MAX)
                continue;
            end
            
            t = bdata_time;
            D = squeeze(bdata_raw{tt}(trial,:,:));
            
            [currentA, voltageA, currentB, voltageB] = get_dual_scaled_voltage_and_current( D );
                        
            cur_idx = trial_idx_cntrs(tt);
            
            volts_all{tt}(cur_idx, :) = voltageA;
            current_all{tt}(cur_idx, :) = currentA;

            t_volts{tt}(cur_idx,:) = t-TIME_OFFSET;
            
            [ t_vel, vel_forward, vel_side, vel_yaw ] = get_velocity_ephysrig(t, D, datapath, 1 );            
            
            t_vel_all{tt}( cur_idx, : ) = t_vel-TIME_OFFSET;            
            yaw_vel_all{tt}( cur_idx, : ) = vel_yaw;
            forward_vel_all{tt}( cur_idx, : ) = vel_forward + FWD_VEL_OFFSET;
            trial_idx_cntrs(tt) = trial_idx_cntrs(tt) + 1;
        end
    end
end

%% Calculate the optimal shift between yaw and firing rate.
CALC_PSTH = 1;
generate_optimal_shift(sid, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH, analysis_path_type);


%% generate Vm vs. yaw plot for the stimulation period
CALC_PSTH = 1;

generate_odor_evoked_Vm_vs_yaw_plot_LAL_DNs(sid, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH);

%% generate a Vm vs. yaw trajectory plot
CALC_PSTH = 1;
image_cnts = generate_Vm_vs_yaw_trajectory_LAL_DNs_v2(sid, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH, analysis_path_type);

%% generate a Vm vs. yaw trajectory plot with fwd vel overlay 
CALC_PSTH = 1;
generate_Vm_vs_yaw_trajectory_LAL_DNs_v3_fwd_overlay(sid, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH, analysis_path_type);

%% generate a fwd vs. yaw trajectory plot with FR overlay
CALC_PSTH = 1;
generate_fwd_vs_yaw_trajectory_LAL_DNs(sid, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH, analysis_path_type);

%% Calculate FR threshold for left turning
DOWN_YAW = 25;
DOWN_FR = 5;
image_cnts_reshaped = reshape(image_cnts, [DOWN_FR, size(image_cnts,1)/DOWN_FR, DOWN_YAW, size(image_cnts,2)/DOWN_YAW]);
image_cnts_sum_d = squeeze(sum(squeeze(sum(image_cnts_reshaped, 1)), 2));

%figure;
%histogram( image_cnts_sum_d(8,:) );
%hold on;

FR_s = 40;

FR_start = FR_s/DOWN_FR;
Y = image_cnts_sum_d(FR_start,:);

p_val_table = zeros(size(image_cnts_sum_d,1));

for FR_s = 5:5:250
    FR_start = FR_s/DOWN_FR;
    Y = image_cnts_sum_d(FR_start,:);

    for i=1:size(image_cnts_sum_d,1)
        Y1 = image_cnts_sum_d(i,:);
        %X1 = [1:length(Y1)];
        %f1 = fit( X1.', Y1.', 'gauss1');
        % plot(f1, X1, Y1);
                
        %disp(['ranksum( ' num2str(i) ' ): ' num2str( ranksum(Y, Y1) ) ' firing rate: ' num2str(DOWN_FR*i)]);
        [h,p] = ttest2(Y, Y1, 'Tail', 'both');
        %disp(['FR_start:  ' num2str(FR_s) '  ttest2( ' num2str(i) ' ): (h: ' num2str( h ) '   p: ' num2str(p) ' ) firing rate: ' num2str(DOWN_FR*i)]);
        p_val_table( FR_s, i ) = p;
    end
end

p_val_table(find(p_val_table <= 0.05)) = 0.0;

f = figure;
imagesc([1:250], [1:250], p_val_table);
axis image;
ylabel('Fixed FR (spikes/s)');
xlabel('Checked FR (spikes/s)');

caxis([0 0.1])
colorbar;

saveas(f, [analysis_path '/fixed_checked_p_value_plot.fig']);
saveas(f, [analysis_path '/fixed_checked_p_value_plot.png']);

% Say ~ 50 spikes/s ???

%% Calculate the yaw
WINDOW_CNT = 10000;
WINDOW_SIZE = 100; % samples @ 100 Hz

ITER_CNT = 1000;

yaw_rand_NF = zeros(1,ITER_CNT);
fwd_rand_NF = zeros(1,ITER_CNT);

for ii = 1:ITER_CNT
    
    yaw_wins = zeros(WINDOW_CNT, WINDOW_SIZE);
    fwd_wins = zeros(WINDOW_CNT, WINDOW_SIZE);
    
    for i = 1:WINDOW_CNT
        tt = ceil(rand * trial_type_cnt);
        
        trials_in_tt = size(yaw_vel_all{tt},1);
        num_samples_in_trial = size(yaw_vel_all{tt},2);
        
        rand_trial_idx = ceil(rand * trials_in_tt);
        
        % Get a random number in the range
        a = WINDOW_SIZE/2;
        b = num_samples_in_trial-WINDOW_SIZE/2;
        r = ceil((b-a).*rand) + a;
        
        rand_range = [(r-WINDOW_SIZE/2):((r+WINDOW_SIZE/2)-1)];
        
        cur_yaw_win = yaw_vel_all{tt}( rand_trial_idx, rand_range );
        cur_fwd_win = forward_vel_all{tt}( rand_trial_idx, rand_range );
        
        yaw_wins(i,:) = cur_yaw_win;
        fwd_wins(i,:) = cur_fwd_win;
    end
    
%     fwd_rand_avg(ii) = mean(squeeze(mean(fwd_wins)));
%     yaw_rand_avg(ii) = mean(squeeze(mean(yaw_wins)));

    fwd_rand_NF(ii) = mean( std(fwd_wins,1,1) .* 1.96 );
    yaw_rand_NF(ii) = mean( std(yaw_wins,1,1) .* 1.96 );

end

f = figure;
NBins = 20;
subplot(2,1,1);
histogram( fwd_rand_NF, NBins );
xlabel('fwd vel (mm/s)');
ylabel('Count');

subplot(2,1,2);
histogram( yaw_rand_NF, NBins );
xlabel('yaw vel (deg/s)');
ylabel('Count');

title(['Yaw threshold for turning: ' num2str(mean(yaw_rand_NF))]);

saveas(f, [analysis_path '/running_baseline_hist_windows_' num2str(WINDOW_CNT) '_iter_' num2str(ITER_CNT) ' .fig']);
saveas(f, [analysis_path '/running_baseline_hist_windows_' num2str(WINDOW_CNT) '_iter_' num2str(ITER_CNT) ' .png']);


% Yaw threshold for turning: 425.8076
%%
f1 = figure;
t = [0:WINDOW_SIZE-1]/100;

subplot(2,1,1);
hold on;

avg_fwd = mean( fwd_wins );
NF_fwd = std( fwd_wins,1 ) * 1.96;
fh = fill( [t, fliplr(t)], ... 
        [(avg_fwd+NF_fwd) fliplr((avg_fwd-NF_fwd))], ...
        rgb('Salmon'));

plot( t, avg_fwd, 'b' );

ylabel('Fwd (mm/s)');

title(['Avg of random ' num2str(WINDOW_CNT) ' windows ']);
set(fh, 'EdgeColor', 'None');


subplot(2,1,2);
hold on;

avg_yaw = mean( yaw_wins );
NF_yaw = std( fwd_wins,1 ) * 1.96;

fh = fill( [t, fliplr(t)], ... 
        [(avg_yaw+NF_yaw) fliplr((avg_yaw-NF_yaw))], ...
        rgb('Salmon'));
set(fh, 'EdgeColor', 'None');

plot( t, avg_yaw, 'b' );

ylabel('Yaw (deg/s)');
xlabel('Time (s)');

saveas(f1, [analysis_path '/fwd_yaw_rand_windows_' num2str(WINDOW_CNT) '.fig' ]);
saveas(f1, [analysis_path '/fwd_yaw_rand_windows_' num2str(WINDOW_CNT) '.png' ]);


%% Plot binned YAW and FR image, also plot the FR vs. left turn bias plot.

DOWN_YAW = 25;
DOWN_FR = 10;
image_cnts_reshaped = reshape(image_cnts, [DOWN_FR, size(image_cnts,1)/DOWN_FR, DOWN_YAW, size(image_cnts,2)/DOWN_YAW]);
image_cnts_sum_d = squeeze(sum(squeeze(sum(image_cnts_reshaped, 1)), 2));

f = figure; 
imagesc([-2000 2000], [0 250], flipdim(image_cnts_sum_d,1));
set(gca,'YTickLabel',flipud(get(gca,'YTickLabel')));
colormap jet;
caxis([0 100]);
xlabel('Yaw (deg/s)');
ylabel('Firing rate (spikes/s)');
title(['FR bin size: ' num2str(DOWN_FR) '   Yaw bin size: ' num2str(DOWN_YAW) ]);
set(gca, 'FontSize', 16);

saveas(f, [analysis_path '/trajectory_hist_img.fig' ]);
saveas(f, [analysis_path '/trajectory_hist_img.png' ]);


left_right_bias = zeros(1, size(image_cnts_sum_d,1));

BIAS_PADDING_FROM_ZERO = 4;
zero_yaw_line = size(image_cnts_sum_d,2) / 2;
for i = 1:size(image_cnts_sum_d,1)
    
    cur_left_sum = sum( squeeze( image_cnts_sum_d( i, 1 : zero_yaw_line-BIAS_PADDING_FROM_ZERO ) ) );
    cur_right_sum = sum( squeeze( image_cnts_sum_d( i, zero_yaw_line+BIAS_PADDING_FROM_ZERO : end ) ) );
    
    cur_left_right_bias = (cur_left_sum-cur_right_sum) / (cur_left_sum+cur_right_sum);
    left_right_bias(i) = cur_left_right_bias;   
end

f1 = figure; 
%plot(left_right_bias, DOWN_FR * [1:size(image_cnts_sum_d,1)]);
plot( DOWN_FR * [1:size(image_cnts_sum_d,1)], left_right_bias, 'LineWidth', 2.0 );
ylabel('Left Turn Index [(Left-Right)/(Left+Right)]');
xlabel('Firing rate (spikes/s)');
title(['FR bin size: ' num2str(DOWN_FR) ]);

set(gca, 'FontSize', 16);

saveas(f1, [analysis_path '/trajectory_FR_vs_bias_plot.fig' ]);
saveas(f1, [analysis_path '/trajectory_FR_vs_bias_plot.png' ]);
