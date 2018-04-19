% Load all data from all the animals to be included in this analysis

clear all;
close all;

working_dir = '/data/drive1/sasha/';

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directories_to_analyze_A2 =  { { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_03', 0, 1900 }, ...                          
                               { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_04', 1, 2346 }, ...                           
                               { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_05', 0, 1740 }, ... 
                               { '161210_lexAOpCsChrimson_gfp_83blexA_ss730_06', 0, 1882 } };

%directories_to_analyze = { { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_05', 0, 1740 } };                        
%directories_to_analyze = { { '161218_lexAOpCsChrimson_gfp_83blexA_ss731_02', 0, 3000 } }; 

%                              { '170721_lexAOpCsChrimson_gfp_83blexA_ss731_05', 2, 2829 }, ...                           

CELL_TYPE = 'A1';
directories_to_analyze_A1 = { { '161218_lexAOpCsChrimson_gfp_83blexA_ss731_02', 0, 3000 }, ...                          
                              { '170726_lexAOpCsChrimson_gfp_83blexA_ss731_08', 0, 1656 }, ...
                              { '170727_lexAOpCsChrimson_gfp_83blexA_ss731_12', 0, 2800 }};

pre_stim_t = 3.0;
stim_t =  3.5;
post_stim_t = 3.0;
inter_trial_t = 5.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timenow_str = datestr(datetime, 'yymmdd_HHMMSS');
analysis_path = [working_dir '/summary_analysis/A1_A2_multi_fly_analysis/'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

[t_all_A2, t_vel_all_A2, yaw_all_A2, fwd_all_A2, ephys_all_A2, dummy] = load_LAL_DN_data( working_dir, directories_to_analyze_A2, ephys_SR, ball_SR, 0 );
[t_all_A1, t_vel_all_A1, yaw_all_A1, fwd_all_A1, ephys_all_A1, dummy] = load_LAL_DN_data( working_dir, directories_to_analyze_A1, ephys_SR, ball_SR, 0 );

SAVE_SUFFIX_STR = '';

%% Calculate PSTH and downsample for scatter plots

num_flies_A1 = length( directories_to_analyze_A1 );
num_flies_A2 = length( directories_to_analyze_A2 );
psth_dt_samples = ephys_SR/ball_SR;

psth_A1         = cell(1, num_flies_A1);
A1_psth_down    = cell(1, num_flies_A1);
fwd_all_down_A1 = cell(1, num_flies_A1);
yaw_all_down_A1 = cell(1, num_flies_A1);


psth_A2         = cell(1, num_flies_A2);
A2_psth_down    = cell(1, num_flies_A2);
fwd_all_down_A2 = cell(1, num_flies_A2);
yaw_all_down_A2 = cell(1, num_flies_A2);

BIN_SIZE = 0.05;
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

for f = 1:num_flies_A1
    SPIKE_THRESHOLD_LAL_DN = 2.0;
    tic; psth_A1{f} = calculate_psth_A1( t_all_A1{f}, t_vel_all_A1{f}, ephys_all_A1{f}, ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples ); toc;

    A1_psth_down{f} = squeeze(mean(reshape( psth_A1{f}, [ DT_YAW, length( psth_A1{f} ) / DT_YAW ] ), 1));
    fwd_all_down_A1{f} = squeeze(mean(reshape( fwd_all_A1{f}, [ DT_YAW, length( yaw_all_A1{f} )/DT_YAW ]), 1));
    yaw_all_down_A1{f} = squeeze(mean(reshape( yaw_all_A1{f}, [ DT_YAW, length( yaw_all_A1{f} ) / DT_YAW ]), 1));
end

for f = 1:num_flies_A2
    SPIKE_THRESHOLD_LAL_DN = 0.3;
    tic; psth_A2{f} = calculate_psth_A2( t_all_A2{f}, t_vel_all_A2{f}, ephys_all_A2{f}, ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples ); toc;

    A2_psth_down{f} = squeeze(mean(reshape( psth_A2{f}, [ DT_YAW, length( psth_A2{f})/DT_YAW ] ), 1));
    fwd_all_down_A2{f} = squeeze(mean(reshape( fwd_all_A2{f}, [ DT_YAW, length( yaw_all_A2{f} )/DT_YAW ]), 1));
    yaw_all_down_A2{f} = squeeze(mean(reshape( yaw_all_A2{f}, [ DT_YAW, length( yaw_all_A2{f} )/DT_YAW ]), 1));
end

A1_psth_all = [];
A1_yaw_all  = [];
A1_fwd_all  = [];

A2_psth_all = [];
A2_yaw_all  = [];
A2_fwd_all  = [];

SHIFT_FACTOR = 3;

for f = 1:num_flies_A1
    A1_psth_all = horzcat( A1_psth_all, A1_psth_down{f}(1:end-SHIFT_FACTOR+1));
    A1_yaw_all  = horzcat( A1_yaw_all,  yaw_all_down_A1{f}(SHIFT_FACTOR:end));
    A1_fwd_all  = horzcat( A1_fwd_all,  fwd_all_down_A1{f}(SHIFT_FACTOR:end));        
end

for f = 1:num_flies_A2
    A2_psth_all = horzcat( A2_psth_all, A2_psth_down{f}(1:end-SHIFT_FACTOR+1));
    A2_yaw_all  = horzcat( A2_yaw_all,  yaw_all_down_A2{f}(SHIFT_FACTOR:end));
    A2_fwd_all  = horzcat( A2_fwd_all,  fwd_all_down_A2{f}(SHIFT_FACTOR:end));        
end


%% Exclude data where the fly is standing still
SAVE_SUFFIX_STR = 'exclude_standing';

% Copy all to a temporary variable
A1_psth_down_tmp    = A1_psth_all;
fwd_all_down_A1_tmp = A1_fwd_all;
yaw_all_down_A1_tmp = A1_yaw_all;

A2_psth_down_tmp    = A2_psth_all;
fwd_all_down_A2_tmp = A2_fwd_all;
yaw_all_down_A2_tmp = A2_yaw_all;


FWD_THRESHOLD = 1.0;
A1_psth_all = [];
A1_fwd_all = [];
A1_yaw_all = [];

for i = 1:length(fwd_all_down_A1_tmp)
    
    cur_fwd = fwd_all_down_A1_tmp(i);
    if( (cur_fwd < -1.0*FWD_THRESHOLD) || (cur_fwd > FWD_THRESHOLD))
        A1_psth_all(end+1) = A1_psth_down_tmp(i);
        A1_fwd_all(end+1) = fwd_all_down_A1_tmp(i);
        A1_yaw_all(end+1) = yaw_all_down_A1_tmp(i);
    end    
end

A2_psth_all = [];
A2_fwd_all = [];
A2_yaw_all = [];

for i = 1:length(fwd_all_down_A2_tmp)
    
    cur_fwd = fwd_all_down_A2_tmp(i);
    if( (cur_fwd < -1.0*FWD_THRESHOLD) || (cur_fwd > FWD_THRESHOLD))
        A2_psth_all(end+1) = A2_psth_down_tmp(i);
        A2_fwd_all(end+1) = fwd_all_down_A2_tmp(i);
        A2_yaw_all(end+1) = yaw_all_down_A2_tmp(i);
    end    
end


%% Include running straight only 
SAVE_SUFFIX_STR = 'include_running_straight_only';

% Copy all to a temporary variable
A1_psth_down_tmp    = A1_psth_all;
fwd_all_down_A1_tmp = A1_fwd_all;
yaw_all_down_A1_tmp = A1_yaw_all;

A2_psth_down_tmp    = A2_psth_all;
fwd_all_down_A2_tmp = A2_fwd_all;
yaw_all_down_A2_tmp = A2_yaw_all;


FWD_THRESHOLD = 1.0;
YAW_THRESHOLD = 350;
A1_psth_all = [];
A1_fwd_all = [];
A1_yaw_all = [];

for i = 1:length(fwd_all_down_A1_tmp)    
    cur_fwd = fwd_all_down_A1_tmp(i);
    cur_yaw = yaw_all_down_A1_tmp(i);
    
    if( ((cur_yaw > -1.0*YAW_THRESHOLD) && (cur_yaw < YAW_THRESHOLD)) && (cur_fwd > FWD_THRESHOLD) )
        A1_psth_all(end+1) = A1_psth_down_tmp( i );
        A1_fwd_all(end+1) = fwd_all_down_A1_tmp( i );
        A1_yaw_all(end+1) = yaw_all_down_A1_tmp( i );
    end    
end

A2_psth_all = [];
A2_fwd_all = [];
A2_yaw_all = [];

for i = 1:length(fwd_all_down_A2_tmp)    
    cur_fwd = fwd_all_down_A2_tmp(i);
    cur_yaw = yaw_all_down_A2_tmp(i);
    
    if( ((cur_yaw > -1.0*YAW_THRESHOLD) && (cur_yaw < YAW_THRESHOLD)) && (cur_fwd > FWD_THRESHOLD) )
        A2_psth_all(end+1) = A2_psth_down_tmp( i );
        A2_fwd_all(end+1) = fwd_all_down_A2_tmp( i );
        A2_yaw_all(end+1) = yaw_all_down_A2_tmp( i );
    end    
end



%% Plot 2D histograms 
%f = figure;
f = figure('units','normalized','outerposition',[0 0 1 1]);

CAXIS_LIMIT_YAW = 50;
CAXIS_LIMIT_FWD = CAXIS_LIMIT_YAW;

YAW_LIM = 1000;

subplot(2,2,1);
hist3([A1_psth_all', A1_yaw_all'], 'NBins', [50 100], 'CDataMode','auto','FaceColor','interp' );
xlabel('A1 PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
zlabel('count');
colormap gray;
colorbar
caxis([0 CAXIS_LIMIT_YAW]);
view(2)
axis tight;
ylim([-1.0*YAW_LIM YAW_LIM]);
title(['Yaw (yaw shift factor: ' num2str(SHIFT_FACTOR_YAW*BIN_SIZE) ' s)']);

subplot(2,2,2);
hist3([A1_psth_all', A1_fwd_all'], 'NBins', [50 100], 'CDataMode','auto','FaceColor','interp' );
xlabel('A1 PSTH (spikes/s)');
ylabel('Fwd (mm/s)');
zlabel('count');
colorbar
caxis([0 CAXIS_LIMIT_FWD]);
view(2)
axis tight;
ylim([-40 80]);
zlabel('count');
title(['Fwd (fwd shift factor: ' num2str(SHIFT_FACTOR_FWD*BIN_SIZE) ' s)']);

subplot(2,2,3);
hist3([A2_psth_all', A2_yaw_all'], 'NBins', [50 100], 'CDataMode','auto','FaceColor','interp' );
xlabel('A2 PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
zlabel('count');
colorbar
caxis([0 CAXIS_LIMIT_YAW]);
view(2)
axis tight;
ylim([-1.0*YAW_LIM YAW_LIM]);
zlabel('count');

subplot(2,2,4);
hist3([A2_psth_all', A2_fwd_all'], 'NBins', [50 100], 'CDataMode','auto','FaceColor','interp' );
xlabel('A2 PSTH (spikes/s)');
ylabel('Fwd (mm/s)');
zlabel('count');
colorbar
caxis([0 CAXIS_LIMIT_FWD]);
view(2)
axis tight;
ylim([-40 80]);
zlabel('count');

savestr = '/A1_A2_psth_vs_yaw_and_fwd_hist2_all_flies';
savestr = [ savestr '_' SAVE_SUFFIX_STR ];

saveas(f, [analysis_path savestr '.fig']);
saveas(f, [analysis_path savestr '.png']);


%% Plot x,y weighted averages from 2D histograms above

f = figure('units','normalized','outerposition',[0 0 1 1]);

% A1
[N, c] = hist3( [A1_psth_all', A1_yaw_all'], 'NBins', [50 100], 'CDataMode', 'auto', 'FaceColor', 'interp' );

psth_N_weighted = N .* repmat(c{1}', [1 size(N,2)]);
yaw_N_weighted = N .* repmat(c{2}, [size(N,1) 1]);

weight_sum_yaw = sum(N,1);
weight_sum_psth = sum(N,2);

A1_psth_all_wavg_yaw_x = squeeze(sum(psth_N_weighted, 1)) ./ weight_sum_yaw;
yaw_all_wavg_psth_x = sum(yaw_N_weighted, 2) ./ weight_sum_psth;

subplot(2,4,1);
plot(c{2}, A1_psth_all_wavg_yaw_x);
ylim([0 100]);
xlim([-1000 1000]);
xlabel('Yaw (deg/s)');
ylabel('A1 PSTH (spikes/s)');

subplot(2,4,2);
plot(c{1}, yaw_all_wavg_psth_x');
xlim([0 100]);
ylim([-500 500]);
ylabel('Yaw (deg/s)');
xlabel('A1 PSTH (spikes/s)');

[N,c] = hist3([A1_psth_all', A1_fwd_all'], 'NBins', [50 100], 'CDataMode','auto','FaceColor','interp' );

psth_N_weighted = N .* repmat(c{1}', [1 size(N,2)]);
fwd_N_weighted = N .* repmat(c{2}, [size(N,1) 1]);

weight_sum_fwd = sum(N,1);
weight_sum_psth = sum(N,2);

A1_psth_all_wavg_fwd_x = squeeze(sum(psth_N_weighted, 1)) ./ weight_sum_fwd;
fwd_all_wavg_psth_x = sum(fwd_N_weighted, 2) ./ weight_sum_psth;

subplot(2,4,3);
plot(c{2}, A1_psth_all_wavg_fwd_x);
ylim([0 100]);
xlim([-10 20]);
xlabel('Fwd (mm/s)');
ylabel('A1 PSTH (spikes/s)');

subplot(2,4,4);
plot(c{1}, fwd_all_wavg_psth_x');
xlim([0 100]);
ylim([-10 20]);
ylabel('Fwd (mm/s)');
xlabel('A1 PSTH (spikes/s)');

% A2
[N, c] = hist3( [A2_psth_all', A2_yaw_all'], 'NBins', [50 100], 'CDataMode', 'auto', 'FaceColor', 'interp' );

psth_N_weighted = N .* repmat(c{1}', [1 size(N,2)]);
yaw_N_weighted = N .* repmat(c{2}, [size(N,1) 1]);

weight_sum_yaw = sum(N,1);
weight_sum_psth = sum(N,2);

A2_psth_all_wavg_yaw_x = squeeze(sum(psth_N_weighted, 1)) ./ weight_sum_yaw;
yaw_all_wavg_psth_x = sum(yaw_N_weighted, 2) ./ weight_sum_psth;

subplot(2,4,5);
plot(c{2}, A2_psth_all_wavg_yaw_x);
ylim([0 100]);
xlim([-1000 1000]);
xlabel('Yaw (deg/s)');
ylabel('A2 PSTH (spikes/s)');

subplot(2,4,6);
plot(c{1}, yaw_all_wavg_psth_x');
xlim([0 100]);
ylim([-500 500]);
ylabel('Yaw (deg/s)');
xlabel('A2 PSTH (spikes/s)');

[N,c] = hist3([A2_psth_all', A2_fwd_all'], 'NBins', [50 100], 'CDataMode','auto','FaceColor','interp' );

psth_N_weighted = N .* repmat(c{1}', [1 size(N,2)]);
fwd_N_weighted = N .* repmat(c{2}, [size(N,1) 1]);

weight_sum_fwd = sum(N,1);
weight_sum_psth = sum(N,2);

A2_psth_all_wavg_fwd_x = squeeze(sum(psth_N_weighted, 1)) ./ weight_sum_fwd;
fwd_all_wavg_psth_x = sum(fwd_N_weighted, 2) ./ weight_sum_psth;

subplot(2,4,7);
plot(c{2}, A2_psth_all_wavg_fwd_x);
ylim([0 100]);
xlim([-10 20]);
xlabel('Fwd (mm/s)');
ylabel('A2 PSTH (spikes/s)');

subplot(2,4,8);
plot(c{1}, fwd_all_wavg_psth_x');
xlim([0 100]);
ylim([-10 20]);
ylabel('Fwd (mm/s)');
xlabel('A2 PSTH (spikes/s)');

savestr = '/A1_A2_psth_vs_yaw_and_fwd_weighted_avg_all_flies';
savestr = [ savestr '_' SAVE_SUFFIX_STR ];

saveas(f, [analysis_path savestr '.fig']);
saveas(f, [analysis_path savestr '.png']);


