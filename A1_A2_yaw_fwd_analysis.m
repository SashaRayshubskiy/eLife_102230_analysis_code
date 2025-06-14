% Load all data from all the animals to be included in this analysis

clear all;
close all;

working_dir = '/data/drive1/sasha/';

%%% WARNING: MAKE SURE EPHYS SAMPLE RATE IS CORRECT  4000 Hz
%%% WARNING: MAKE SURE EPHYS SAMPLE RATE IS CORRECT  4000 Hz
%%% WARNING: MAKE SURE EPHYS SAMPLE RATE IS CORRECT  4000 Hz

settings = sensor_settings;
%settings.sampRate = 4000;
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

%% Plot yaw vs. fwd vel for each fly

num_flies_A1 = length( directories_to_analyze_A1 );
num_flies_A2 = length( directories_to_analyze_A2 );

for i = 1:num_flies_A1
    f = figure;
    c=[1 10 100];
    
    [N,cen] = hist3( [ yaw_all_A1{i}, fwd_all_A1{i} ], 'NBins', [50 50] );
    
    %hold on;
    imagesc(cen{1}([1 end]), cen{2}([1 end]), log(N'));
    
    colormap( gray );    
    caxis(log([c(1) c(length(c))]));
    h = colorbar('FontSize',11,'YTick',log(c),'YTickLabel',c);
    ylabel(h, 'Counts');
    
    xlabel('Yaw (deg/s)');
    ylabel('Fwd (mm/s)');
    zlabel('count');
    set(gca,'YDir','normal');
    grid on;
    
    saveas(f,[analysis_path '/A1_yaw_vs_fwd_scatter_logcolor_v2_gray_fly_' num2str(i) '.fig']);
    saveas(f,[analysis_path '/A1_yaw_vs_fwd_scatter_logcolor_v2_gray_fly_' num2str(i) '.png']);
end

for i = 1:num_flies_A2
    f = figure;
    c=[1 10 100];
    
    [N,cen] = hist3( [ yaw_all_A2{i}, fwd_all_A2{i} ], 'NBins', [50 50] );
    
    %hold on;
    imagesc(cen{1}([1 end]), cen{2}([1 end]), log(N'));
    
    colormap( gray );    
    caxis(log([c(1) c(length(c))]));
    h = colorbar('FontSize',11,'YTick',log(c),'YTickLabel',c);
    ylabel(h, 'Counts');
    
    xlabel('Yaw (deg/s)');
    ylabel('Fwd (mm/s)');
    zlabel('count');
    set(gca,'YDir','normal');
    grid on;
    
    saveas(f,[analysis_path '/A2_yaw_vs_fwd_scatter_logcolor_v2_gray_fly_' num2str(i) '.fig']);
    saveas(f,[analysis_path '/A2_yaw_vs_fwd_scatter_logcolor_v2_gray_fly_' num2str(i) '.png']);
end



%% Calculate PSTH and downsample for scatter plots

num_flies_A1 = length( directories_to_analyze_A1 );
num_flies_A2 = length( directories_to_analyze_A2 );
psth_dt_samples = ephys_SR/ball_SR;

psth_A1         = cell(1, num_flies_A1);
A1_psth_down    = cell(1, num_flies_A1);
fwd_all_down_A1 = cell(1, num_flies_A1);
fwd_std_all_down_A1 = cell(1, num_flies_A1);
yaw_all_down_A1 = cell(1, num_flies_A1);

psth_A2         = cell(1, num_flies_A2);
A2_psth_down    = cell(1, num_flies_A2);
fwd_all_down_A2 = cell(1, num_flies_A2);
fwd_std_all_down_A2 = cell(1, num_flies_A2);
yaw_all_down_A2 = cell(1, num_flies_A2);

BIN_SIZE = 0.05;
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

for f = 1:num_flies_A1
    SPIKE_THRESHOLD_LAL_DN = 2.0;
    tic; psth_A1{f} = calculate_psth_A1( t_all_A1{f}, t_vel_all_A1{f}, ephys_all_A1{f}, ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples ); toc;

    A1_psth_down{f} = squeeze(mean(reshape( psth_A1{f}, [ DT_YAW, length( psth_A1{f} ) / DT_YAW ] ), 1));
    fwd_all_down_A1{f} = squeeze(mean(reshape( fwd_all_A1{f}, [ DT_YAW, length( yaw_all_A1{f} )/DT_YAW ]), 1));
    fwd_std_all_down_A1{f} = squeeze(std(reshape( fwd_all_A1{f}, [ DT_YAW, length( yaw_all_A1{f} )/DT_YAW ]), 0, 1));
    yaw_all_down_A1{f} = squeeze(mean(reshape( yaw_all_A1{f}, [ DT_YAW, length( yaw_all_A1{f} ) / DT_YAW ]), 1));
end

for f = 1:num_flies_A2
    SPIKE_THRESHOLD_LAL_DN = 0.3;
    tic; psth_A2{f} = calculate_psth_A2( t_all_A2{f}, t_vel_all_A2{f}, ephys_all_A2{f}, ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples ); toc;

    A2_psth_down{f} = squeeze(mean(reshape( psth_A2{f}, [ DT_YAW, length( psth_A2{f})/DT_YAW ] ), 1));
    fwd_all_down_A2{f} = squeeze(mean(reshape( fwd_all_A2{f}, [ DT_YAW, length( yaw_all_A2{f} )/DT_YAW ]), 1));
    fwd_std_all_down_A2{f} = squeeze(std(reshape( fwd_all_A2{f}, [ DT_YAW, length( yaw_all_A2{f} )/DT_YAW ]), 0, 1));
    yaw_all_down_A2{f} = squeeze(mean(reshape( yaw_all_A2{f}, [ DT_YAW, length( yaw_all_A2{f} )/DT_YAW ]), 1));
end

A1_psth_all = [];
A1_yaw_all  = [];
A1_fwd_all  = [];
A1_fwd_std_all  = [];

A2_psth_all = [];
A2_yaw_all  = [];
A2_fwd_all  = [];
A2_fwd_std_all  = [];

SHIFT_FACTOR = 3;

for f = 1:num_flies_A1
    A1_psth_all = horzcat( A1_psth_all, A1_psth_down{f}(1:end-SHIFT_FACTOR+1));
    A1_yaw_all  = horzcat( A1_yaw_all,  yaw_all_down_A1{f}(SHIFT_FACTOR:end));
    A1_fwd_all  = horzcat( A1_fwd_all,  fwd_all_down_A1{f}(SHIFT_FACTOR:end));        
    A1_fwd_std_all  = horzcat( A1_fwd_std_all,  fwd_std_all_down_A1{f}(1:end-SHIFT_FACTOR+1));
end


for f = 1:num_flies_A2
    A2_psth_all = horzcat( A2_psth_all, A2_psth_down{f}(1:end-SHIFT_FACTOR+1));
    A2_yaw_all  = horzcat( A2_yaw_all,  yaw_all_down_A2{f}(SHIFT_FACTOR:end));
    A2_fwd_all  = horzcat( A2_fwd_all,  fwd_all_down_A2{f}(SHIFT_FACTOR:end));        
    A2_fwd_std_all  = horzcat( A2_fwd_std_all,  fwd_std_all_down_A2{f}(1:end-SHIFT_FACTOR+1));
end


% Exclude data where the fly is standing still
SAVE_SUFFIX_STR = 'exclude_standing';

% Copy all to a temporary variable
A1_psth_down_tmp        = A1_psth_all;
fwd_all_down_A1_tmp     = A1_fwd_all;
fwd_std_all_down_A1_tmp = A1_fwd_std_all;
yaw_all_down_A1_tmp     = A1_yaw_all;

A2_psth_down_tmp        = A2_psth_all;
fwd_all_down_A2_tmp     = A2_fwd_all;
fwd_std_all_down_A2_tmp = A2_fwd_std_all;
yaw_all_down_A2_tmp     = A2_yaw_all;

FWD_STD_THRESHOLD = 0.01;
A1_psth_all = [];
A1_fwd_all = [];
A1_std_fwd_all = [];
A1_yaw_all = [];

i=1;
while( i < length(fwd_all_down_A1_tmp) )
    
    cur_fwd_std = fwd_std_all_down_A1_tmp(i);
    
    if( cur_fwd_std > FWD_STD_THRESHOLD )
        A1_psth_all(end+1) = A1_psth_down_tmp(i);
        A1_fwd_all(end+1) = fwd_all_down_A1_tmp(i);
        A1_yaw_all(end+1) = yaw_all_down_A1_tmp(i);
    end
    
    i = i + 1;
end

A2_psth_all = [];
A2_fwd_all = [];
A2_fwd_std_all = [];
A2_yaw_all = [];

i=1;
while( i < length(fwd_all_down_A2_tmp) )
    
    cur_fwd_std = fwd_std_all_down_A2_tmp(i);
    
    if( cur_fwd_std > FWD_STD_THRESHOLD )
        A2_psth_all(end+1) = A2_psth_down_tmp(i);
        A2_fwd_all(end+1) = fwd_all_down_A2_tmp(i);
        A2_yaw_all(end+1) = yaw_all_down_A2_tmp(i);
    end
    
    i = i + 1;
end


% Exclude standing for each fly
N_FLIES = 3;
FWD_STD_THRESHOLD = 0.01;
A1_psth_per_fly = cell(1,N_FLIES);
A1_fwd_per_fly = cell(1,N_FLIES);
A1_std_fwd_per_fly = cell(1,N_FLIES);
A1_yaw_per_fly = cell(1,N_FLIES);

A2_psth_per_fly = cell(1,N_FLIES);
A2_fwd_per_fly = cell(1,N_FLIES);
A2_fwd_std_per_fly = cell(1,N_FLIES);
A2_yaw_per_fly = cell(1,N_FLIES);

for fl = 1:3
    A1_psth_per_fly{fl} = [];
    A1_fwd_per_fly{fl} = [];
    A1_std_fwd_per_fly{fl} = [];
    A1_yaw_per_fly{fl} = [];
    
    A2_psth_per_fly{fl} = [];
    A2_fwd_per_fly{fl} = [];
    A2_fwd_std_per_fly{fl} = [];
    A2_yaw_per_fly{fl} = [];
end

for fl = 1:3
    
    i=1;
    while( i < length(fwd_std_all_down_A1{fl}) )
        
        cur_fwd_std = fwd_std_all_down_A1{fl}(i);
        
        if( cur_fwd_std > FWD_STD_THRESHOLD )
            A1_psth_per_fly{fl}(end+1) = A1_psth_down{fl}(i);
            A1_fwd_per_fly{fl}(end+1) = fwd_all_down_A1{fl}(i);
            A1_yaw_per_fly{fl}(end+1) = yaw_all_down_A1{fl}(i);
        end
        
        i = i + 1;
    end
    
    j=1;
    while( j < length(fwd_std_all_down_A2{fl}) )
        
        cur_fwd_std = fwd_std_all_down_A2{fl}(j);
        
        if( cur_fwd_std > FWD_STD_THRESHOLD )
            A2_psth_per_fly{fl}(end+1) = A2_psth_down{fl}(j);
            A2_fwd_per_fly{fl}(end+1) = fwd_all_down_A2{fl}(j);
            A2_yaw_per_fly{fl}(end+1) = yaw_all_down_A2{fl}(j);
        end
        
        j = j + 1;
    end
end
% RUN THIS TO RESTART THE ANALYSIS

%% Scatter plots for each fly

%f = figure;
f = figure('units','normalized','outerposition',[0 0 1 1]);

SHIFT_FACTOR = 3;

CAXIS_LIMIT_YAW = 30;
CAXIS_LIMIT_FWD = CAXIS_LIMIT_YAW;

if(strcmp(SAVE_SUFFIX_STR,'exclude_standing') == 1)
    YAW_LIM = 1500;
else
    YAW_LIM = YAW_THRESHOLD;
end

for fl = 1:3
    
    cur_psth_A1 = A1_psth_per_fly{fl}(1:end-SHIFT_FACTOR+1)';
    cur_yaw_A1 = A1_yaw_per_fly{fl}(SHIFT_FACTOR:end)';
    cur_fwd_A1 = A1_fwd_per_fly{fl}(SHIFT_FACTOR:end)';

    cur_psth_A2 = A2_psth_per_fly{fl}(1:end-SHIFT_FACTOR+1)';
    cur_yaw_A2 = A2_yaw_per_fly{fl}(SHIFT_FACTOR:end)';
    cur_fwd_A2 = A2_fwd_per_fly{fl}(SHIFT_FACTOR:end)';

    
    
    subplot(2,2,1);
    hold on;
    scatter( cur_psth_A1, cur_yaw_A1, 1);
    xlabel('A1 PSTH (spikes/s)');
    ylabel('Yaw (deg/s)');
    grid on;
    axis tight;
    ylim([-1.0*YAW_LIM YAW_LIM]);
    title(['Yaw (shift: ' num2str(SHIFT_FACTOR*BIN_SIZE) ' s) Avg PSTH: ' num2str(mean(A1_psth_all)) ' Avg yaw: ' num2str(mean(A1_yaw_all))]);
    
    [N, c] = hist3( [cur_psth_A1, cur_yaw_A1], 'NBins', [50 100], 'CDataMode', 'auto', 'FaceColor', 'interp' );
    
    psth_N_weighted = N .* repmat(c{1}', [1 size(N,2)]);
    yaw_N_weighted = N .* repmat(c{2}, [size(N,1) 1]);
    
    weight_sum_yaw = sum(N,1);
    weight_sum_psth = sum(N,2);
    
    A1_psth_all_wavg_yaw_x = squeeze(sum(psth_N_weighted, 1)) ./ weight_sum_yaw;
    yaw_all_wavg_psth_x = sum(yaw_N_weighted, 2) ./ weight_sum_psth;
     
    A1_FR_yaw(fl,:,:) = [c{1}; yaw_all_wavg_psth_x'];
    
    subplot(2,2,2);
    hold on;
    scatter( cur_psth_A1, cur_fwd_A1, 1 );
    xlabel('A1 PSTH (spikes/s)');
    ylabel('Fwd (mm/s)');
    grid on;
    axis tight;
    ylim([-40 80]);
    title(['Fwd (fwd shift factor: ' num2str(SHIFT_FACTOR*BIN_SIZE) ' s) Avg PSTH: ' num2str(mean(A1_psth_all)) ' Avg fwd: ' num2str(mean(A1_fwd_all))]);
    
    [N,c] = hist3([cur_psth_A1, cur_fwd_A1], 'NBins', [50 100], 'CDataMode','auto','FaceColor','interp' );
    
    psth_N_weighted = N .* repmat(c{1}', [1 size(N,2)]);
    fwd_N_weighted = N .* repmat(c{2}, [size(N,1) 1]);
    
    weight_sum_fwd = sum(N,1);
    weight_sum_psth = sum(N,2);
    
    A1_psth_all_wavg_fwd_x = squeeze(sum(psth_N_weighted, 1)) ./ weight_sum_fwd;
    fwd_all_wavg_psth_x = sum(fwd_N_weighted, 2) ./ weight_sum_psth;
    
    A1_FR_fwd(fl,:,:) = [c{1}; fwd_all_wavg_psth_x'];
    
    subplot(2,2,3);
    hold on;
    scatter( cur_psth_A2, cur_yaw_A2, 1);
    xlabel('A2 PSTH (spikes/s)');
    ylabel('Yaw (deg/s)');
    grid on;
    axis tight;
    ylim([-1.0*YAW_LIM YAW_LIM]);
    title(['Avg PSTH: ' num2str(mean(A2_psth_all)) ' Avg yaw: ' num2str(mean(A2_yaw_all))]);
    
    [N, c] = hist3( [cur_psth_A2, cur_yaw_A2], 'NBins', [50 100], 'CDataMode', 'auto', 'FaceColor', 'interp' );
    
    psth_N_weighted = N .* repmat(c{1}', [1 size(N,2)]);
    yaw_N_weighted = N .* repmat(c{2}, [size(N,1) 1]);
    
    weight_sum_yaw = sum(N,1);
    weight_sum_psth = sum(N,2);
    
    A2_psth_all_wavg_yaw_x = squeeze(sum(psth_N_weighted, 1)) ./ weight_sum_yaw;
    yaw_all_wavg_psth_x = sum(yaw_N_weighted, 2) ./ weight_sum_psth;

    A2_FR_yaw(fl,:,:) = [c{1}; yaw_all_wavg_psth_x'];
    
    subplot(2,2,4);
    hold on;
    scatter( cur_psth_A2, cur_fwd_A2, 1 );
    xlabel('A2 PSTH (spikes/s)');
    ylabel('Fwd (mm/s)');
    grid on;
    axis tight;
    ylim([-40 80]);
    title(['Avg PSTH: ' num2str(mean(A2_psth_all)) ' Avg fwd: ' num2str(mean(A2_fwd_all))]);
   
    [N, c] = hist3( [cur_psth_A2, cur_fwd_A2], 'NBins', [50 100], 'CDataMode', 'auto', 'FaceColor', 'interp' );
    
    psth_N_weighted = N .* repmat(c{1}', [1 size(N,2)]);
    fwd_N_weighted = N .* repmat(c{2}, [size(N,1) 1]);
    
    weight_sum_fwd = sum(N,1);
    weight_sum_psth = sum(N,2);
    
    A2_psth_all_wavg_yaw_x = squeeze(sum(psth_N_weighted, 1)) ./ weight_sum_fwd;
    fwd_all_wavg_psth_x = sum(fwd_N_weighted, 2) ./ weight_sum_psth;

    A2_FR_fwd(fl,:,:) = [c{1}; fwd_all_wavg_psth_x'];
end

savestr = '/A1_A2_psth_vs_yaw_and_fwd_scatter_each_fly';
savestr = [ savestr '_' SAVE_SUFFIX_STR ];

saveas(f, [analysis_path savestr '.fig']);
saveas(f, [analysis_path savestr '.png']);


f = figure;

for fl = 1:3
    subplot(2,1,1);
    hold on;
    plot( squeeze(A1_FR_yaw(fl,1,:)), squeeze(A1_FR_yaw(fl,2,:)), 'r' );
    plot( squeeze(A2_FR_yaw(fl,1,:)), squeeze(A2_FR_yaw(fl,2,:)), 'g' );
%     plot( squeeze(mean(squeeze(A1_FR_yaw(:,1,:)))), squeeze(mean(squeeze(A1_FR_yaw(:,2,:)))), 'r', 'LineWidth', 3 );
%     plot( squeeze(mean(squeeze(A2_FR_yaw(:,1,:)))), squeeze(mean(squeeze(A2_FR_yaw(:,2,:)))), 'g', 'LineWidth', 3  );

    
    xlim([0 100]);
    ylim([-500 100]);
    ylabel('Yaw (deg/s)');
    xlabel('Firing rate (spikes/s)');
    legend('A1', 'A2');
    
    subplot(2,1,2);
    hold on;
    plot( squeeze(A1_FR_fwd(fl,1,:)), squeeze(A1_FR_fwd(fl,2,:)), 'r' );
    plot( squeeze(A2_FR_fwd(fl,1,:)), squeeze(A2_FR_fwd(fl,2,:)), 'g' );
%     plot( squeeze(mean(squeeze(A1_FR_fwd(:,1,:)))), squeeze(mean(squeeze(A1_FR_fwd(:,2,:)))), 'r', 'LineWidth', 3 );
%     plot( squeeze(mean(squeeze(A2_FR_fwd(:,1,:)))), squeeze(mean(squeeze(A2_FR_fwd(:,2,:)))), 'g', 'LineWidth', 3  );
%    legend('A1', 'A2');
    xlim([0 100]);
    ylim([-10 15]);
    ylabel('Fwd (mm/s)');
    xlabel('Firing rate (spikes/s)');
end

subplot(2,1,1);
plot( squeeze(A1_FR_yaw(1,1,:)), squeeze(mean(squeeze(A1_FR_yaw(:,2,:)))), 'r', 'LineWidth', 3 );
plot( squeeze(A2_FR_yaw(1,1,:)), squeeze(mean(squeeze(A2_FR_yaw(:,2,:)))), 'g', 'LineWidth', 3 );
legend('A1', 'A2');

subplot(2,1,2);
plot( squeeze(A1_FR_fwd(1,1,:)), squeeze(mean(squeeze(A1_FR_fwd(:,2,:)))), 'r', 'LineWidth', 3 );
plot( squeeze(A2_FR_fwd(1,1,:)), squeeze(mean(squeeze(A2_FR_fwd(:,2,:)))), 'g', 'LineWidth', 3 );


saveas(f, [analysis_path savestr '_fit.fig']);
saveas(f, [analysis_path savestr '_fit.png']);

%% Include running straight only 
YAW_THRESHOLD = 200;
SAVE_SUFFIX_STR = ['include_running_straight_fwd_back_only_yaw_cutoff_' num2str(YAW_THRESHOLD)];

% Copy all to a temporary variable
A1_psth_down_tmp    = A1_psth_all;
fwd_all_down_A1_tmp = A1_fwd_all;
yaw_all_down_A1_tmp = A1_yaw_all;

A2_psth_down_tmp    = A2_psth_all;
fwd_all_down_A2_tmp = A2_fwd_all;
yaw_all_down_A2_tmp = A2_yaw_all;


A1_psth_all = [];
A1_fwd_all = [];
A1_yaw_all = [];

for i = 1:length(fwd_all_down_A1_tmp)    
    cur_yaw = yaw_all_down_A1_tmp(i);
    
    if( ((cur_yaw > -1.0*YAW_THRESHOLD) && (cur_yaw < YAW_THRESHOLD)) )
        A1_psth_all(end+1) = A1_psth_down_tmp( i );
        A1_fwd_all(end+1) = fwd_all_down_A1_tmp( i );
        A1_yaw_all(end+1) = yaw_all_down_A1_tmp( i );
    end    
end

A2_psth_all = [];
A2_fwd_all = [];
A2_yaw_all = [];

for i = 1:length(fwd_all_down_A2_tmp)    
    cur_yaw = yaw_all_down_A2_tmp(i);
    
    if( ((cur_yaw > -1.0*YAW_THRESHOLD) && (cur_yaw < YAW_THRESHOLD)) )
        A2_psth_all(end+1) = A2_psth_down_tmp( i );
        A2_fwd_all(end+1) = fwd_all_down_A2_tmp( i );
        A2_yaw_all(end+1) = yaw_all_down_A2_tmp( i );
    end    
end



%% Plot 2D histograms - log color scale
%f = figure;
f = figure('units','normalized','outerposition',[0 0 1 1]);

CAXIS_LIMIT_YAW = 30;
CAXIS_LIMIT_FWD = CAXIS_LIMIT_YAW;

if(strcmp(SAVE_SUFFIX_STR,'exclude_standing') == 1)
    YAW_LIM = 1000;
else
    YAW_LIM = YAW_THRESHOLD;    
end
c=[1 10 100 300];

subplot(2,2,1);

[N,cen] = hist3( [A1_psth_all', A1_yaw_all'], 'NBins', [50 100] );

imagesc(cen{1}([1 end]), cen{2}([1 end]), log(N'));

colormap(gray); 
caxis(log([c(1) c(length(c))]));
h = colorbar('FontSize',11,'YTick',log(c),'YTickLabel',c);
ylabel(h, 'log(Counts)');

xlabel('A1 PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
zlabel('count');
%view(2)
%axis tight;
%ylim([-1.0*YAW_LIM YAW_LIM]);
title(['Yaw (shift: ' num2str(SHIFT_FACTOR*BIN_SIZE) ' s) Avg PSTH: ' num2str(mean(A1_psth_all)) ' Avg yaw: ' num2str(mean(A1_yaw_all))]);
%set(gca(), 'YTickLabel', flip(get(gca(), 'YTickLabel')));
set(gca(), 'ydir', 'reverse');

subplot(2,2,2);

[N,cen] = hist3([A1_psth_all', A1_fwd_all'], 'NBins', [50 100], 'CDataMode','auto','FaceColor','interp' );
imagesc(cen{1}([1 end]), cen{2}([1 end]), log(N'));

colormap(gray); 
caxis(log([c(1) c(length(c))]));
h = colorbar('FontSize',11,'YTick',log(c),'YTickLabel',c);
ylabel(h, 'log(Counts)');

xlabel('A1 PSTH (spikes/s)');
ylabel('Fwd (mm/s)');
zlabel('count');
view(2)
axis tight;
ylim([-40 80]);
zlabel('count');
title(['Fwd (fwd shift factor: ' num2str(SHIFT_FACTOR*BIN_SIZE) ' s) Avg PSTH: ' num2str(mean(A1_psth_all)) ' Avg fwd: ' num2str(mean(A1_fwd_all))]);
%set(gca(), 'YTickLabel', flip(get(gca(), 'YTickLabel')));

subplot(2,2,3);

[N,cen] = hist3([A2_psth_all', A2_yaw_all'], 'NBins', [50 100], 'CDataMode','auto','FaceColor','interp' );
imagesc(cen{1}([1 end]), cen{2}([1 end]), log(N'));

colormap(gray); 
caxis(log([c(1) c(length(c))]));
h = colorbar('FontSize',11,'YTick',log(c),'YTickLabel',c);
ylabel(h, 'log(Counts)');

xlabel('A2 PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
zlabel('count');
view(2)
axis tight;
ylim([-1.0*YAW_LIM YAW_LIM]);
title(['Avg PSTH: ' num2str(mean(A2_psth_all)) ' Avg yaw: ' num2str(mean(A2_yaw_all))]);
%set(gca(), 'YTickLabel', flip(get(gca(), 'YTickLabel')));

subplot(2,2,4);

[N,cen] = hist3([A2_psth_all', A2_fwd_all'], 'NBins', [50 100], 'CDataMode','auto','FaceColor','interp' );

imagesc(cen{1}([1 end]), cen{2}([1 end]), log(N'));
colormap(gray); 

caxis(log([c(1) c(length(c))]));
h = colorbar('FontSize',11,'YTick',log(c),'YTickLabel',c);
ylabel(h, 'log(Counts)');

xlabel('A2 PSTH (spikes/s)');
ylabel('Fwd (mm/s)');
zlabel('count');
view(2)
axis tight;
ylim([-40 80]);
zlabel('count');
title(['Avg PSTH: ' num2str(mean(A2_psth_all)) ' Avg fwd: ' num2str(mean(A2_fwd_all))]);
%set(gca(), 'YTickLabel', flip(get(gca(), 'YTickLabel')));

savestr = '/A1_A2_psth_vs_yaw_and_fwd_hist2_log_caxis_all_flies';
savestr = [ savestr '_' SAVE_SUFFIX_STR ];

saveas(f, [analysis_path savestr '.fig']);
saveas(f, [analysis_path savestr '.png']);

%% Plot 2D histograms 
%f = figure;
f = figure('units','normalized','outerposition',[0 0 1 1]);

CAXIS_LIMIT_YAW = 30;
CAXIS_LIMIT_FWD = CAXIS_LIMIT_YAW;

if(strcmp(SAVE_SUFFIX_STR,'exclude_standing') == 1)
    YAW_LIM = 1000;
else
    YAW_LIM = YAW_THRESHOLD;    
end

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
title(['Yaw (shift: ' num2str(SHIFT_FACTOR*BIN_SIZE) ' s) Avg PSTH: ' num2str(mean(A1_psth_all)) ' Avg yaw: ' num2str(mean(A1_yaw_all))]);

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
title(['Fwd (fwd shift factor: ' num2str(SHIFT_FACTOR*BIN_SIZE) ' s) Avg PSTH: ' num2str(mean(A1_psth_all)) ' Avg fwd: ' num2str(mean(A1_fwd_all))]);

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
title(['Avg PSTH: ' num2str(mean(A2_psth_all)) ' Avg yaw: ' num2str(mean(A2_yaw_all))]);

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
title(['Avg PSTH: ' num2str(mean(A2_psth_all)) ' Avg fwd: ' num2str(mean(A2_fwd_all))]);

savestr = '/A1_A2_psth_vs_yaw_and_fwd_hist2_all_flies';
savestr = [ savestr '_' SAVE_SUFFIX_STR ];

saveas(f, [analysis_path savestr '.fig']);
saveas(f, [analysis_path savestr '.png']);

%% Plot 2D scatter plots 
%f = figure;
f = figure('units','normalized','outerposition',[0 0 1 1]);

CAXIS_LIMIT_YAW = 30;
CAXIS_LIMIT_FWD = CAXIS_LIMIT_YAW;

if(strcmp(SAVE_SUFFIX_STR,'exclude_standing') == 1)
    YAW_LIM = 1500;
else
    YAW_LIM = YAW_THRESHOLD;    
end

subplot(2,2,1);
scatter( A1_psth_all', A1_yaw_all', 1);
xlabel('A1 PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
grid on;
axis tight;
ylim([-1.0*YAW_LIM YAW_LIM]);
title(['Yaw (shift: ' num2str(SHIFT_FACTOR*BIN_SIZE) ' s) Avg PSTH: ' num2str(mean(A1_psth_all)) ' Avg yaw: ' num2str(mean(A1_yaw_all))]);

subplot(2,2,2);
scatter( A1_psth_all', A1_fwd_all', 1 );
xlabel('A1 PSTH (spikes/s)');
ylabel('Fwd (mm/s)');
grid on;
axis tight;
ylim([-40 80]);
title(['Fwd (fwd shift factor: ' num2str(SHIFT_FACTOR*BIN_SIZE) ' s) Avg PSTH: ' num2str(mean(A1_psth_all)) ' Avg fwd: ' num2str(mean(A1_fwd_all))]);

subplot(2,2,3);
scatter(A2_psth_all', A2_yaw_all', 1 );
xlabel('A2 PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
grid on;
axis tight;
ylim([-1.0*YAW_LIM YAW_LIM]);
title(['Avg PSTH: ' num2str(mean(A2_psth_all)) ' Avg yaw: ' num2str(mean(A2_yaw_all))]);

subplot(2,2,4);
scatter(A2_psth_all', A2_fwd_all', 1 );
xlabel('A2 PSTH (spikes/s)');
ylabel('Fwd (mm/s)');
grid on;
axis tight;
ylim([-40 80]);
title(['Avg PSTH: ' num2str(mean(A2_psth_all)) ' Avg fwd: ' num2str(mean(A2_fwd_all))]);

savestr = '/A1_A2_psth_vs_yaw_and_fwd_scatter_all_flies';
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
hh = title(SAVE_SUFFIX_STR);
set(hh, 'Interpreter', 'none');

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

%% Plot the relationship between yaw and fwd vel

CAXIS_LIMIT_FWD1 = 100;


if 0
subplot(1,2,1);
plot( A1_yaw_all, A1_fwd_all, 'o', 'MarkerSize', 1 );
xlabel('Yaw (deg/s)');
ylabel('Fwd (mm/s)');
title('A1 (excludes standing)');
grid on;

subplot(2,2,2);
hist3([A1_yaw_all', A1_fwd_all'], 'NBins', [50 50], 'CDataMode','auto','FaceColor','interp' );
xlabel('Yaw (deg/s)');
ylabel('Fwd (mm/s)');
zlabel('count');
colormap gray;
colorbar
caxis([0 CAXIS_LIMIT_FWD1]);
view(2)
axis tight;
%ylim([-40 80]);
zlabel('count');
title(['Avg fwd: ' num2str(mean(A1_fwd_all)) ' Avg yaw: ' num2str(mean(A1_yaw_all))]);
end

if 0
f = figure;
subplot(1,1,1);
plot( A2_yaw_all, A2_fwd_all, 'o', 'MarkerSize', 1 );
xlabel('Yaw (deg/s)');
ylabel('Fwd (mm/s)');
title('A2 (excludes standing)');
grid on;
saveas(f, [analysis_path 'A1_A2_yaw_vs_fwd_scatter.fig']);
saveas(f, [analysis_path 'A1_A2_yaw_vs_fwd_scatter.png']);
end

f = figure;
c=[1 10 100];

[N,cen] = hist3( [A2_yaw_all', A2_fwd_all'], 'NBins', [50 50] );

imagesc(cen{1}([1 end]), cen{2}([1 end]), log(N'));

USE_JET = 0;
if( USE_JET == 1 )
    colormap( jet ); 
else
    colormap( gray ); 
end

caxis(log([c(1) c(length(c))]));
h = colorbar('FontSize',11,'YTick',log(c),'YTickLabel',c);
ylabel(h, 'Counts');

xlabel('Yaw (deg/s)');
ylabel('Fwd (mm/s)');
zlabel('count');
set(gca,'YDir','normal');

if( USE_JET == 1 )
    saveas(f,[analysis_path '/A1_A2_yaw_vs_fwd_scatter_logcolor_v2_jet.fig']);
    saveas(f,[analysis_path '/A1_A2_yaw_vs_fwd_scatter_logcolor_v2_jet.png']);
else
    saveas(f,[analysis_path '/A1_A2_yaw_vs_fwd_scatter_logcolor_v2_gray.fig']);
    saveas(f,[analysis_path '/A1_A2_yaw_vs_fwd_scatter_logcolor_v2_gray.png']);
end

if 0
subplot(1,1,1);
hist3([A2_yaw_all', A2_fwd_all'], 'NBins', [50 50], 'CDataMode','auto','FaceColor','interp' );
xlabel('Yaw (deg/s)');
ylabel('Fwd (mm/s)');
zlabel('count');
colormap gray;
colorbar
caxis([0 CAXIS_LIMIT_FWD1]);
view(2)
axis tight;
%ylim([-40 80]);
zlabel('count');
title(['Avg fwd: ' num2str(mean(A2_fwd_all)) ' Avg yaw: ' num2str(mean(A2_yaw_all))]);
end








%% Include running backward only 
YAW_THRESHOLD = 200;
BACK_THRESHOLD = -0.1;
SAVE_SUFFIX_STR = ['include_running_back_only_yaw_cutoff_' num2str(250) '_back_cutoff_' num2str(BACK_THRESHOLD)];

% Copy all to a temporary variable
A1_psth_down_tmp    = A1_psth_all;
fwd_all_down_A1_tmp = A1_fwd_all;
yaw_all_down_A1_tmp = A1_yaw_all;

A2_psth_down_tmp    = A2_psth_all;
fwd_all_down_A2_tmp = A2_fwd_all;
yaw_all_down_A2_tmp = A2_yaw_all;


A1_psth_all = [];
A1_fwd_all = [];
A1_yaw_all = [];

for i = 1:length(fwd_all_down_A1_tmp)    
    cur_yaw = yaw_all_down_A1_tmp(i);
    cur_fwd = fwd_all_down_A1_tmp(i);
    
    if( ((cur_yaw > -1.0*YAW_THRESHOLD) && (cur_yaw < YAW_THRESHOLD)) && (cur_fwd < BACK_THRESHOLD ))
        A1_psth_all(end+1) = A1_psth_down_tmp( i );
        A1_fwd_all(end+1) = fwd_all_down_A1_tmp( i );
        A1_yaw_all(end+1) = yaw_all_down_A1_tmp( i );
    end    
end

A2_psth_all = [];
A2_fwd_all = [];
A2_yaw_all = [];

for i = 1:length(fwd_all_down_A2_tmp)    
    cur_yaw = yaw_all_down_A2_tmp(i);
    cur_fwd = fwd_all_down_A2_tmp(i);
    
    if( ((cur_yaw > -1.0*YAW_THRESHOLD) && (cur_yaw < YAW_THRESHOLD)) && (cur_fwd < BACK_THRESHOLD ))
        A2_psth_all(end+1) = A2_psth_down_tmp( i );
        A2_fwd_all(end+1) = fwd_all_down_A2_tmp( i );
        A2_yaw_all(end+1) = yaw_all_down_A2_tmp( i );
    end    
end


%% Include running forward only 
YAW_THRESHOLD = 200;
FWD_THRESHOLD = 0.1;
SAVE_SUFFIX_STR = ['include_running_fwd_only_yaw_cutoff_' num2str(250) '_fwd_cutoff_' num2str(FWD_THRESHOLD)];

% Copy all to a temporary variable
A1_psth_down_tmp    = A1_psth_all;
fwd_all_down_A1_tmp = A1_fwd_all;
yaw_all_down_A1_tmp = A1_yaw_all;

A2_psth_down_tmp    = A2_psth_all;
fwd_all_down_A2_tmp = A2_fwd_all;
yaw_all_down_A2_tmp = A2_yaw_all;


A1_psth_all = [];
A1_fwd_all = [];
A1_yaw_all = [];

for i = 1:length(fwd_all_down_A1_tmp)    
    cur_yaw = yaw_all_down_A1_tmp(i);
    cur_fwd = fwd_all_down_A1_tmp(i);
    
    if( ((cur_yaw > -1.0*YAW_THRESHOLD) && (cur_yaw < YAW_THRESHOLD)) && (cur_fwd > FWD_THRESHOLD ))
        A1_psth_all(end+1) = A1_psth_down_tmp( i );
        A1_fwd_all(end+1) = fwd_all_down_A1_tmp( i );
        A1_yaw_all(end+1) = yaw_all_down_A1_tmp( i );
    end    
end

A2_psth_all = [];
A2_fwd_all = [];
A2_yaw_all = [];

for i = 1:length(fwd_all_down_A2_tmp)    
    cur_yaw = yaw_all_down_A2_tmp(i);
    cur_fwd = fwd_all_down_A2_tmp(i);
    
    if( ((cur_yaw > -1.0*YAW_THRESHOLD) && (cur_yaw < YAW_THRESHOLD)) && (cur_fwd > FWD_THRESHOLD ))
        A2_psth_all(end+1) = A2_psth_down_tmp( i );
        A2_fwd_all(end+1) = fwd_all_down_A2_tmp( i );
        A2_yaw_all(end+1) = yaw_all_down_A2_tmp( i );
    end    
end












