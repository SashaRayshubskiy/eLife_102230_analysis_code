%% Load all data from all the animals to be included in this analysis

clear all;
close all;

%working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/descending_neurons/';

working_dir = '/data/drive1/sasha/';

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

dirname = '180410_gfp_3G_ss730_dual_08';

analysis_path = [working_dir dirname '/analysis/'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

directories_to_analyze = { { '180410_gfp_3G_ss730_dual_08', [0], 11.5*300 } };

[t_all, t_vel_all, yaw_all, fwd_all, ephys_all_A, ephys_all_B] = load_LAL_DN_data( working_dir, directories_to_analyze, ephys_SR, ball_SR );

idx = 1;


%% Plot yaw, fwd, ephysA,B in a time window

XMIN = 648;
XMAX = 660;

figure;

ax(1) = subplot(2,1,1);

yyaxis left;
hold on;
plot(t_vel_all{1}, yaw_all{1});

yyaxis right;
hold on;
plot(t_vel_all{1}, fwd_all{1});

xlim([XMIN XMAX]);

ax(2) = subplot(2,1,2);
hold on;
plot(t_all{1}, ephys_all_A{1}, 'b');
plot(t_all{1}, ephys_all_B{1}, 'g');
xlim([XMIN XMAX]);


%% Extract Vm (filter out the spikes)
FILT_FACTOR = 0.4;

VmFilt_A2_L = medfilt1( ephys_all_A{1}, FILT_FACTOR * ephys_SR, 'truncate' );
VmFilt_A2_R = medfilt1( ephys_all_B{1}, FILT_FACTOR * ephys_SR, 'truncate' );

START = 1;
END   = 10000*10;

% figure;
% hold on;
% plot( t_all{1}(START:END), VmFilt_A2_L(START:END), 'b' );
% plot( t_all{1}(START:END), VmFilt_A2_R(START:END), 'g' );

figure;
hold on;
plot( t_all{1}(START:END), ephys_all_A{1}(START:END), 'b' );
plot( t_all{1}(START:END), VmFilt_A2_L(START:END), 'g' );



%% Hamming filter
START = 1;
END   = 10000*10;

HAM_WINDOW_SIZE = 4000;

f1 = fir1(HAM_WINDOW_SIZE, 0.0002, 'low');
VmFilt_A2_L_hamm = filter(f1, 1, VmFilt_A2_L);
VmFilt_A2_R_hamm = filter(f1, 1, VmFilt_A2_R);

% figure;
% hold on;
% plot( t_all{1}(START:END), VmFilt_A2_L(START:END), 'b' );
% plot( t_all{1}(START:END-HAM_WINDOW_SIZE/2), VmFilt_A2_L_hamm(START+HAM_WINDOW_SIZE/2:END), 'g' );

figure;
hold on;
plot( t_all{1}(START:END-HAM_WINDOW_SIZE/2), VmFilt_A2_L_hamm(START+HAM_WINDOW_SIZE/2:END), 'b' );
plot( t_all{1}(START:END-HAM_WINDOW_SIZE/2), VmFilt_A2_R_hamm(START+HAM_WINDOW_SIZE/2:END), 'g' );


%% Downsample the data

BIN_SIZE = 0.050; % s
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

t_down = squeeze(mean(reshape(t_all{1}, [DT_EPHYS, length(t_all{1})/DT_EPHYS]), 1));
A2_Vm_L_down = squeeze(mean(reshape(VmFilt_A2_L, [ DT_EPHYS, length(VmFilt_A2_L)/DT_EPHYS ] ),1));
A2_Vm_R_down = squeeze(mean(reshape(VmFilt_A2_R, [ DT_EPHYS, length(VmFilt_A2_R)/DT_EPHYS ] ),1));

yaw_all_down = squeeze(mean(reshape(yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));
fwd_all_down = squeeze(mean(reshape(fwd_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));

%% Generate Left A2 vs. Right A2 scatter plot

f = figure;

plot(A2_Vm_R_down, A2_Vm_L_down, 'o', 'MarkerSize', 3);

saveas(f,[analysis_path '/right_vs_left_A2_scatter_plot.fig']);
saveas(f,[analysis_path '/right_vs_left_A2_scatter_plot.png']);

%% Correlation in window plot

WINDOW_SIZE             = 0.25; % seconds
WINDOW_SIZE_IN_SAMPLES  = WINDOW_SIZE * ephys_SR; 

WINDOW_STEP_SIZE = 0.05; % seconds;
WINDOW_STEP_SIZE_IN_SAMPLES = WINDOW_STEP_SIZE * ephys_SR;

step_count = length(VmFilt_A2_L) / WINDOW_STEP_SIZE_IN_SAMPLES;
cur_step = 1;

A2_L_R_corr = zeros( 1, step_count);

STEPS_IN_WINDOW = WINDOW_SIZE_IN_SAMPLES /  WINDOW_STEP_SIZE_IN_SAMPLES;

for s = 1:step_count-STEPS_IN_WINDOW-1
    
    cur_win = [ cur_step : (cur_step+WINDOW_SIZE_IN_SAMPLES) ];
    
    cur_r = corr( VmFilt_A2_L(cur_win), VmFilt_A2_R(cur_win) );
    
    A2_L_R_corr(s) = cur_r;
    
    cur_step = cur_step + WINDOW_STEP_SIZE_IN_SAMPLES;
end


%% Show hamming window time courses along with correlation
START_T = 648;
END_T = 656;

START_BALL = START_T * ball_SR;
END_BALL = END_T * ball_SR;

START_EPYS = START_T * ephys_SR;
END_EPHYS = END_T * ephys_SR;

f = figure; 

ax(1) = subplot(2,1,1);
yyaxis left;
hold on;
%plot(t_vel_all{1}(START_BALL:END_BALL), yaw_all{1}(START_BALL:END_BALL));
plot(t_vel_all{1}, yaw_all{1});
ylabel('Yaw (deg/s)');

yyaxis right;
hold on;
% plot(t_vel_all{1}(START_BALL:END_BALL), fwd_all{1}(START_BALL:END_BALL));
plot(t_vel_all{1}, fwd_all{1});
ylabel('Fwd (mm/s)');
xlim([START_T, END_T]);

ax(2) = subplot(2,1,2);
yyaxis left;
hold on;
% plot(t_all{1}(START_EPYS:END_EPHYS), VmFilt_A2_L(START_EPYS:END_EPHYS), '-b');
% plot(t_all{1}(START_EPYS:END_EPHYS), VmFilt_A2_R(START_EPYS:END_EPHYS), '-g');
plot(t_all{1}, VmFilt_A2_L, '-b');
plot(t_all{1}, VmFilt_A2_R, '-g');

yyaxis right;
hold on;

t_down = squeeze(mean(reshape(t_all{1}, [DT_EPHYS, length(t_all{1})/DT_EPHYS]), 1));

plot( t_down, A2_L_R_corr, 'r');
xlim([START_T, END_T]);

linkaxes(ax, 'x');

%saveas(f,[analysis_path '/right_vs_left_A2_hamming_filtered_tc.fig']);
%saveas(f,[analysis_path '/right_vs_left_A2_hamming_filtered_tc.png']);


%% Show A2 L/R correlation histogram 

f = figure;
hist(A2_L_R_corr, 50);
xlabel('Pearsons corr');
ylabel('count');
saveas(f,[analysis_path '/A2_right_left_corr_hist.fig']);
saveas(f,[analysis_path '/A2_right_left_corr_hist.png']);


%% Calculate PSTH for both channels

SPIKE_THRESHOLD_LAL_DN = 0.25;
psth_dt_samples = ephys_SR/ball_SR;
tic; A2_L_psth = calculate_psth_A2( t_all{1}, t_vel_all{1}, ephys_all_A{1}, ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples ); toc;
tic; A2_R_psth = calculate_psth_A2( t_all{1}, t_vel_all{1}, ephys_all_B{1}, ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples ); toc;

%% Plot yaw, fwd, and PSTH for left and right A2 neuron

XMIN = 648;
XMAX = 660;

figure;

ax(1) = subplot(3,1,1);

yyaxis left;
hold on;
plot(t_vel_all{1}, yaw_all{1});

yyaxis right;
hold on;
plot(t_vel_all{1}, fwd_all{1});

xlim([XMIN XMAX]);

ax(2) = subplot(3,1,2);
hold on;
plot(t_all{1}, ephys_all_A{1}, 'b');
plot(t_all{1}, ephys_all_B{1}, 'g');
xlim([XMIN XMAX]);

ax(3) = subplot(3,1,3);
hold on;
plot(t_vel_all{1}, A2_L_psth, 'b');
plot(t_vel_all{1}, A2_R_psth, 'g');
xlim([XMIN XMAX]);

linkaxes(ax, 'x');

%% Show relationship between left and right A2 neurons PSTH difference and yaw

BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;

A2_L_psth_base = A2_L_psth - mean(A2_L_psth);
A2_R_psth_base = A2_R_psth - mean(A2_R_psth);

A2_L_psth_down = squeeze(mean(reshape( A2_L_psth_base, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ),1));
A2_R_psth_down = squeeze(mean(reshape( A2_R_psth_base, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ),1));
yaw_all_down = squeeze(mean(reshape(yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));

cnt = 0;
A2_L_psth_choosen = zeros(1, length( A2_L_psth_down ));
A2_R_psth_choosen = zeros(1, length( A2_L_psth_down ));
yaw_choosen  = zeros(1, length( yaw_all_down ));

PSTH_BASELINE_CUTOFF = 1.0;
SHIFT_FACTOR = 3;
cnt = 0;
for i = 1:length( A2_L_psth_down )-(SHIFT_FACTOR+1)
    
    cur_A2_L_psth = A2_L_psth_down( i );
    cur_A2_R_psth = A2_R_psth_down( i );
    
    if( (( cur_A2_L_psth <= -1.0*PSTH_BASELINE_CUTOFF ) || ( cur_A2_L_psth >= PSTH_BASELINE_CUTOFF )) && ...
        (( cur_A2_R_psth <= -1.0*PSTH_BASELINE_CUTOFF ) || ( cur_A2_R_psth >= PSTH_BASELINE_CUTOFF )) )        
        A2_L_psth_choosen( cnt+1 ) = cur_A2_L_psth;
        A2_R_psth_choosen( cnt+1 ) = cur_A2_R_psth;
        yaw_choosen( cnt+1 ) = yaw_all_down(i+SHIFT_FACTOR);
    end
        
    cnt = cnt + 1;
end

f = figure;

subplot(2,2,1);
plot( A2_L_psth_choosen(1:cnt), yaw_choosen(1:cnt), 'o', 'MarkerSize', 3 );
xlabel('A2 L PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
title('Left A2 vs. yaw');
grid on;

subplot(2,2,2);
plot( A2_R_psth_choosen(1:cnt), yaw_choosen(1:cnt), 'o', 'MarkerSize', 3 );
xlabel('A2 R PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
title('Right A2 vs. yaw');
grid on;

subplot(2,2,3);
A2_L_R_psth_choosen_diff = A2_L_psth_choosen(1:cnt) - A2_R_psth_choosen(1:cnt);

plot( A2_L_R_psth_choosen_diff, yaw_choosen(1:cnt), 'o', 'MarkerSize', 3 );
xlabel('A2 L-R PSTH diff (spikes/s)');
ylabel('Yaw (deg/s)');
title('diff(Left,Right) A2 vs. yaw');
grid on;

subplot(2,2,4);
%plot( A2_L_psth_choosen(1:cnt) ./ A2_R_psth_choosen(1:cnt), yaw_choosen(1:cnt), 'o', 'MarkerSize', 3 );
semilogx( A2_L_psth_choosen(1:cnt) ./ A2_R_psth_choosen(1:cnt), yaw_choosen(1:cnt), 'o', 'MarkerSize', 3 );
xlabel('A2 L/R PSTH ratio');
ylabel('Yaw (deg/s)');
title('div(Left,Right) A2 vs. yaw');
grid on;

saveas(f,[analysis_path '/A2_right_left_PSTH_diff_div_vs_yaw_scatter.fig']);
saveas(f,[analysis_path '/A2_right_left_PSTH_diff_div_vs_yaw_scatter.png']);

f = figure;
scatter3( A2_L_psth_choosen(1:cnt), yaw_choosen(1:cnt), A2_L_R_psth_choosen_diff );
xlabel('A2 L psth (spike/s)');
ylabel('yaw vel (deg/s)');
zlabel('A2 L-R psth diff(spike/s)');

saveas(f,[analysis_path '/A2_left_PSTH_diff_vs_yaw_scatter3.fig']);
saveas(f,[analysis_path '/A2_left_PSTH_diff_vs_yaw_scatter3.png']);


%% Show relationship between left and right A2 neurons PSTH difference and fwd

BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;

A2_L_R_diff_psth = A2_L_psth - A2_R_psth;

A2_L_R_diff_psth_down = squeeze(mean(reshape( A2_L_R_diff_psth, [ DT_YAW, length(A2_L_R_diff_psth)/DT_YAW ] ),1));
A2_L_psth_down = squeeze(mean(reshape( A2_L_psth, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ),1));
A2_R_psth_down = squeeze(mean(reshape( A2_R_psth, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ),1));
fwd_all_down = squeeze(mean(reshape(fwd_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));

f = figure;

subplot(1,3,1);
plot(A2_L_psth_down(1:end-3), fwd_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('A2 L PSTH (spikes/s)');
ylabel('Fwd (deg/s)');
title('Left A2 vs. fwd');
grid on;

subplot(1,3,2);
plot(A2_R_psth_down(1:end-3), fwd_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('A2 R PSTH (spikes/s)');
ylabel('Fwd (mm/s)');
title('Right A2 vs. fwd');
grid on;

subplot(1,3,3);
plot(A2_L_R_diff_psth_down(1:end-3), fwd_all_down(1:end-3), 'o', 'MarkerSize', 3);
xlabel('A2 L-R PSTH diff (spikes/s)');
ylabel('Fwd (mm/s)');
title('Left-Right A2 vs. fwd');
grid on;

saveas(f,[analysis_path '/A2_right_left_PSTH_diff_vs_fwd_scatter.fig']);
saveas(f,[analysis_path '/A2_right_left_PSTH_diff_vs_fwd_scatter.png']);


%% Plot A2 L,R PSTH with a fwd vel overlay

% Assuming fwd range [-40 60]

RBC_SIZE = 100;

BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;

A2_L_R_diff_psth = A2_L_psth - A2_R_psth;

A2_L_R_diff_psth_down = squeeze(mean(reshape( A2_L_R_diff_psth, [ DT_YAW, length(A2_L_R_diff_psth)/DT_YAW ] ),1));
A2_L_psth_down = squeeze(mean(reshape( A2_L_psth, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ),1));
A2_R_psth_down = squeeze(mean(reshape( A2_R_psth, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ),1));
fwd_all_down = squeeze(mean(reshape(fwd_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));

rbc = jet(RBC_SIZE);
SHIFT_FACTOR = 2;
hold on;

A2_fwd_colors = zeros( length(A2_L_psth_down), 3 );

for i = 1:length(A2_L_psth_down)-SHIFT_FACTOR
    
    cur_fwd = fwd_all_down( i+SHIFT_FACTOR );
    cur_rbc_index = ceil((cur_fwd + 40) * ( RBC_SIZE / 100 ));
    
    if( cur_rbc_index < 1 )
        cur_rbc_index = 1;
    elseif( cur_rbc_index > RBC_SIZE )
        cur_rbc_index = RBC_SIZE;
    end
    
    cur_clr = rbc( cur_rbc_index, : );
    
    A2_fwd_colors( i, : ) = cur_clr;
    
    %plot(A2_L_psth_down(i), A2_R_psth_down(i), 'o', 'MarkerSize', 3, 'color', cur_clr );
end

f = figure;
rbc = jet(RBC_SIZE);
scatter(A2_L_psth_down, A2_R_psth_down, 4, A2_fwd_colors);
colormap('jet')
h = colorbar;
ylabel(h, 'Fwd vel (mm/s)');
caxis([-40 60]);
grid on;
xlabel('A2 Left PSTH (spikes/s)');
ylabel('A2 Right PSTH (spikes/s)');
title('Fwd velocity overlay');

saveas(f,[analysis_path '/A2_right_left_PSTH_scatter_with_fwd_overlay.fig']);
saveas(f,[analysis_path '/A2_right_left_PSTH_scatter_with_fwd_overlay.png']);


%% Plot A2 L,R PSTH with a yaw vel overlay

% Assuming yaw range [-1000 1000]

RBC_SIZE = 2000;

BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;

A2_L_psth_base = A2_L_psth - mean(A2_L_psth);
A2_R_psth_base = A2_R_psth - mean(A2_R_psth);

A2_L_psth_down = squeeze(mean(reshape( A2_L_psth_base, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ),1));
A2_R_psth_down = squeeze(mean(reshape( A2_R_psth_base, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ),1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));

rbc = jet(RBC_SIZE);
SHIFT_FACTOR = 3;
hold on;

A2_yaw_colors = zeros( length(A2_L_psth_down), 3 );

A2_L_choosen = zeros( 1, length(A2_L_psth_down) );
A2_R_choosen = zeros( 1, length(A2_R_psth_down) );
yaw_choosen = zeros( 1, length(A2_R_psth_down) );
A2_yaw_colors_choosen = zeros( length(A2_L_psth_down), 3 );

YAW_CUTOFF = 0;
cnt = 0;
for i = 1:length(A2_L_psth_down)-SHIFT_FACTOR
    
    cur_yaw = yaw_all_down( i+SHIFT_FACTOR );
    cur_rbc_index = ceil(cur_yaw + 1000);
    
    if( cur_rbc_index < 1 )
        cur_rbc_index = 1;
    elseif( cur_rbc_index > RBC_SIZE )
        cur_rbc_index = RBC_SIZE;
    end
    
    cur_clr = rbc( cur_rbc_index, : );
    
    if( (cur_yaw < -1 * YAW_CUTOFF) || (cur_yaw > YAW_CUTOFF) )
    % if( (cur_yaw > -1 * YAW_CUTOFF) && (cur_yaw < YAW_CUTOFF) )
        A2_L_choosen( cnt+1 ) = A2_L_psth_down(i);
        A2_R_choosen( cnt+1 ) = A2_R_psth_down(i);
        yaw_choosen( cnt+1 ) = cur_yaw;
        A2_yaw_colors_choosen( cnt+1, : ) = cur_clr;
        cnt = cnt + 1;
    end    
end

f = figure;
rbc = jet(RBC_SIZE);
%scatter(A2_L_psth_down, A2_R_psth_down, 4, A2_yaw_colors);
scatter( A2_L_choosen, A2_R_choosen, 4, A2_yaw_colors_choosen);
%scatter( A2_L_choosen-A2_R_choosen, yaw_choosen, 3 );
colormap('jet')
h = colorbar;
ylabel(h, 'Yaw vel (deg/s)');
caxis([-1000 1000]);
xlabel('A2 Left PSTH (spikes/s)');
ylabel('A2 Right PSTH (spikes/s)');
grid on;
axis tight;
title(['Yaw velocity overlay: yaw cutoff: ' num2str(YAW_CUTOFF) ]);

saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay.fig']);
saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay.png']);

% Plot a 2D histrogram for A2 L PSTH vs. A2 R PSTH
f = figure;
hist3([A2_L_choosen', A2_R_choosen']);
xlabel('A2 L PSTH (spikes/s)');
ylabel('A2 R PSTH (spikes/s)');
zlabel('count');
saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay_2D_hist.fig']);
saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay_2D_hist.fig']);


%% Show relationship between left and right A2 neurons Vm difference and yaw

BIN_SIZE = 0.05;
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

A2_L_R_diff = VmFilt_A2_L - VmFilt_A2_R;

A2_L_R_diff_down = squeeze(mean(reshape( A2_L_R_diff, [ DT_EPHYS, length(A2_L_R_diff)/DT_EPHYS ] ),1));
A2_L_down = squeeze(mean(reshape( VmFilt_A2_L, [ DT_EPHYS, length(A2_L_R_diff)/DT_EPHYS ] ),1));
A2_R_down = squeeze(mean(reshape( VmFilt_A2_R, [ DT_EPHYS, length(A2_L_R_diff)/DT_EPHYS ] ),1));
yaw_all_down = squeeze(mean(reshape(yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));

f = figure;

subplot(1,3,1);
plot(A2_L_down(1:end-3), yaw_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('A2 L Vm (mV)');
ylabel('Yaw (deg/s)');
title('Left A2 vs. yaw');

subplot(1,3,2);
plot(A2_R_down(1:end-3), yaw_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('A2 R Vm (mV)');
ylabel('Yaw (deg/s)');
title('Right A2 vs. yaw');

subplot(1,3,3);
plot(A2_L_R_diff_down(1:end-3), yaw_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('A2 L-R Vm (mV)');
ylabel('Yaw (deg/s)');
title('Left-Right A2 vs. yaw');

saveas(f,[analysis_path '/A2_right_left_Vm_diff_vs_yaw_scatter.fig']);
saveas(f,[analysis_path '/A2_right_left_Vm_diff_vs_yaw_scatter.png']);


%% Show relationship between left and right A2 neurons Vm difference and fwd

BIN_SIZE = 0.05;
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

A2_L_R_diff = VmFilt_A2_L - VmFilt_A2_R;

A2_L_R_diff_down = squeeze(mean(reshape( A2_L_R_diff, [ DT_EPHYS, length(A2_L_R_diff)/DT_EPHYS ] ),1));
A2_L_down = squeeze(mean(reshape( VmFilt_A2_L, [ DT_EPHYS, length(A2_L_R_diff)/DT_EPHYS ] ),1));
A2_R_down = squeeze(mean(reshape( VmFilt_A2_R, [ DT_EPHYS, length(A2_L_R_diff)/DT_EPHYS ] ),1));
fwd_all_down = squeeze(mean(reshape(fwd_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));

f = figure;

subplot(1,3,1);
plot(A2_L_down(1:end-3), fwd_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('A2 L Vm (mV)');
ylabel('Fwd (mm/s)');
title('Left A2 vs. fwd');
grid on;

subplot(1,3,2);
plot(A2_R_down(1:end-3), fwd_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('A2 R Vm (mV)');
ylabel('Fwd (mm/s)');
title('Right A2 vs. fwd');
grid on;

subplot(1,3,3);
plot(A2_L_R_diff_down(1:end-3), fwd_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('A2 L-R Vm (mV)');
ylabel('Fwd (mm/s)');
title('Left-Right A2 vs. fwd');
grid on;

saveas(f,[analysis_path '/A2_right_left_Vm_diff_vs_fwd_scatter.fig']);
saveas(f,[analysis_path '/A2_right_left_Vm_diff_vs_fwd_scatter.png']);



%% Show relationship between left and right A2 neuron correlation and yaw
BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;

yaw_all_down = squeeze(mean(reshape(yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));

f = figure;

plot(A2_L_R_corr(1:end-3), yaw_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('Pearsons corr (L/R A2 mV)');
ylabel('Yaw (deg/s)');

saveas(f,[analysis_path '/A2_right_left_corr_vs_yaw_scatter.fig']);
saveas(f,[analysis_path '/A2_right_left_corr_vs_yaw_scatter.png']);


%% Show relationship between left and right A2 neuron correlation and fwd
BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;

fwd_all_down = squeeze(mean(reshape(fwd_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));

f = figure;

plot(A2_L_R_corr(1:end-3), fwd_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('Pearsons corr (L/R A2 mV)');
ylabel('Fwd (mm/s)');

saveas(f,[analysis_path '/A2_right_left_corr_vs_fwd_scatter.fig']);
saveas(f,[analysis_path '/A2_right_left_corr_vs_fwd_scatter.png']);

%% Multivariate regression. Is left, right or both A2 better at predicting yaw, fwd vel?

BIN_SIZE = 0.05;
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

YAW_EPHYS_SHIFT = 3;

A2_L_psth_down = squeeze(mean(reshape( A2_L_psth, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ), 1));
A2_R_psth_down = squeeze(mean(reshape( A2_R_psth, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ), 1));
fwd_all_down = squeeze(mean(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]), 1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]), 1));

% Using mvregress()
%Y = horzcat(fwd_all_down', yaw_all_down' );
%myones = ones( length( fwd_all_down ), 1 );
% X = horzcat(myones, A2_L_psth_down', A2_R_psth_down' );
% [ beta, Sigma, E, CovB, logL ] = mvregress(X,Y);


%%%%%%%%
% Using fit()
%%%%%%%%
Y = yaw_all_down(YAW_EPHYS_SHIFT:end)';

% A2 left and right
X = horzcat(A2_L_psth_down(1:end-YAW_EPHYS_SHIFT+1)', A2_R_psth_down(1:end-YAW_EPHYS_SHIFT+1)' );
[fobj, gof] = fit(X,Y, 'poly11');
gof


% A2 left 
X = A2_L_psth_down(1:end-YAW_EPHYS_SHIFT+1)';
[fobj, gof] = fit(X,Y, 'poly1');
gof


% A2 right
X = A2_R_psth_down(1:end-YAW_EPHYS_SHIFT+1)';
[fobj, gof] = fit(X,Y, 'poly1');
gof


%% Remove data points where yaw is 'small' to make it easier to fit.

BIN_SIZE = 0.05;
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

SHIFT_FACTOR = 3;

A2_L_psth_down = squeeze(mean(reshape( A2_L_psth, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ), 1));
A2_R_psth_down = squeeze(mean(reshape( A2_R_psth, [ DT_YAW, length(A2_L_psth)/DT_YAW ] ), 1));
fwd_all_down = squeeze(mean(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]), 1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]), 1));

A2_L_choosen = zeros( 1, length(A2_L_psth_down) );
A2_R_choosen = zeros( 1, length(A2_R_psth_down) );
A2_yaw_choosen = zeros( 1, length(A2_L_psth_down) );
A2_fwd_choosen = zeros( 1, length(A2_L_psth_down) );

if 0
YAW_CUTOFF = 250;
cnt = 0;
for i = 1:length(A2_L_psth_down)-SHIFT_FACTOR
    
    cur_yaw = yaw_all_down( i+SHIFT_FACTOR );
    
    if( (cur_yaw < -1 * YAW_CUTOFF) || (cur_yaw > YAW_CUTOFF) )
        A2_L_choosen( cnt+1 ) = A2_L_psth_down(i);
        A2_R_choosen( cnt+1 ) = A2_R_psth_down(i);
        A2_yaw_choosen( cnt+1 ) = cur_yaw;
        cnt = cnt + 1;
    end    
end
end

YAW_CUTOFF = 250;
BANDPASS = 100;
SHORTPASS = 101;
LONGPASS  = 102;

cutoff_type = BANDPASS;

cnt = 0;
for i = 1:length(A2_L_psth_down)-SHIFT_FACTOR
    
    cur_yaw = yaw_all_down( i+SHIFT_FACTOR );
    
    if( cutoff_type == BANDPASS )
        cut_cond = (cur_yaw < -1.0*YAW_CUTOFF) || (cur_yaw > YAW_CUTOFF);
    elseif( cutoff_type == LONGPASS )
        cut_cond = (cur_yaw > YAW_CUTOFF);
    elseif( cutoff_type == SHORTPASS )
        cut_cond = (cur_yaw < 1.0*YAW_CUTOFF);
    end
    
    if( cut_cond == 1 )
        A2_L_choosen( cnt+1 ) = A2_L_psth_down(i);
        A2_R_choosen( cnt+1 ) = A2_R_psth_down(i);
        A2_yaw_choosen( cnt+1 ) = cur_yaw;
        cnt = cnt + 1;
    end
end

% figure;
% plot( A2_L_choosen, A2_R_choosen, 'o', 'MarkerSize', 3 )

%%%%%%%%
% Using fit() for yaw
%%%%%%%%
Y = A2_yaw_choosen(1:cnt)';

% A2 left and right
X = horzcat(A2_L_choosen(1:cnt)', A2_R_choosen(1:cnt)' );
[fobj, gof_LR_yaw] = fit(X,Y, 'poly11');

% A2 left 
X = A2_L_choosen(1:cnt)';
[fobj, gof_L_yaw] = fit(X,Y, 'poly1');

% A2 right
X = A2_R_choosen(1:cnt)';
[fobj, gof_R_yaw] = fit(X,Y, 'poly1');

A2_L_choosen = zeros( 1, length(A2_L_psth_down) );
A2_R_choosen = zeros( 1, length(A2_R_psth_down) );
A2_fwd_choosen = zeros( 1, length(A2_L_psth_down) );

FWD_CUTOFF = -50;
cnt = 0;

for i = 1:length(A2_L_psth_down)-SHIFT_FACTOR
    
    cur_fwd = fwd_all_down( i+SHIFT_FACTOR );
    
    if(cur_fwd > FWD_CUTOFF)
        A2_L_choosen( cnt+1 ) = A2_L_psth_down(i);
        A2_R_choosen( cnt+1 ) = A2_R_psth_down(i);
        A2_fwd_choosen( cnt+1 ) = cur_fwd;
        cnt = cnt + 1;
    end    
end

%%%%%%%%
% Using fit() for fwd
%%%%%%%%
Y = A2_fwd_choosen(1:cnt)';

% A2 left and right
X = horzcat(A2_L_choosen(1:cnt)', A2_R_choosen(1:cnt)' );
[fobj, gof_LR_fwd] = fit(X,Y, 'poly11');

% A2 left 
X = A2_L_choosen(1:cnt)';
[fobj, gof_L_fwd] = fit(X,Y, 'poly1');

% A2 right
X = A2_R_choosen(1:cnt)';
[fobj, gof_R_fwd] = fit(X,Y, 'poly1');

f = figure;
gofs_yaw = [ gof_LR_yaw.rsquare gof_L_yaw.rsquare gof_R_yaw.rsquare ];
gofs_fwd = [ gof_LR_fwd.rsquare gof_L_fwd.rsquare gof_R_fwd.rsquare ];
bb = bar([gofs_yaw; gofs_fwd]);
ylabel('R^2 val');
    
if( cutoff_type == BANDPASS )
    hh = title(['A2 psth linear regression of ( yaw cutoff: { ' num2str( -1.0*YAW_CUTOFF ) ' , ' num2str(YAW_CUTOFF) ' }, fwd curoff: ' num2str(FWD_CUTOFF) ' )']);
elseif( cutoff_type == LONGPASS )
    hh = title(['A2 psth linear regression of ( yaw cutoff: > { ' num2str( YAW_CUTOFF ) ' } , fwd curoff: ' num2str(FWD_CUTOFF) ' )']);
elseif( cutoff_type == SHORTPASS )
    hh = title(['A2 psth linear regression of ( yaw cutoff: < { ' num2str(-1.0*YAW_CUTOFF) ' }, fwd curoff: ' num2str(FWD_CUTOFF) ' )']);
end
    
set(hh, 'Interpreter', 'none');
set(gca(), 'XTickLabel', {'Yaw', 'Fwd'})

xx = get(gca(), 'Children');
hLegend = legend(xx, {'Right', 'Left', 'Left+Right'});
    
saveas( f, [analysis_path '/A2_psth_right_vs_left_vs_both_rsquared_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.fig'] );
saveas( f, [analysis_path '/A2_psth_right_vs_left_vs_both_rsquared_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.png'] );

f = figure;
gofs_yaw = [ gof_LR_yaw.rmse gof_L_yaw.rmse gof_R_yaw.rmse ];
gofs_fwd = [ gof_LR_fwd.rmse gof_L_fwd.rmse gof_R_fwd.rmse ];
bb = bar([gofs_yaw; gofs_fwd]);
ylabel('RMSE val');
hh = title(['A2 psth linear regression of (yaw cutoff: ' num2str(YAW_CUTOFF) ', fwd curoff: ' num2str(FWD_CUTOFF) ' )']);
set(hh, 'Interpreter', 'none');
set(gca(), 'XTickLabel', {'Yaw', 'Fwd'})

xx = get(gca(), 'Children');
hLegend = legend(xx, {'Right', 'Left', 'Left+Right'});
    
saveas( f, [analysis_path '/A2_psth_right_vs_left_vs_both_rmse_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.fig'] );
saveas( f, [analysis_path '/A2_psth_right_vs_left_vs_both_rmse_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.png'] );

%% Relate the duration of firing rate burst and yaw velocity
% Strategy: find FR peaks

MIN_PEAK = 40;

[~, locs_L, width_L, ~] = findpeaks( A2_L_psth, 'MinPeakHeight', MIN_PEAK, 'Annotate', 'extents');
[~, locs_R, width_R, ~] = findpeaks( A2_R_psth, 'MinPeakHeight', MIN_PEAK, 'Annotate', 'extents');

% figure; 
% 
% yyaxis left;
% hold on;
% plot(t_vel_all{1}, A2_L_psth)
% plot(t_vel_all{1}(locs), A2_L_psth(locs), '+', 'MarkerSize', 3, 'color', 'g');
% 
% yyaxis right;
% 
% hold on;
% plot(t_vel_all{1}(locs), width, '+', 'MarkerSize', 3, 'color', 'r');

yaw_near_peak_L = zeros(1, length(locs_L));

WIN = 5; % each yaw is 0.01 s 

for i = 1:length(locs_L)
    
    cur_loc = locs_L(i);
    
    yaw_near_peak_L(i) = mean(yaw_all{1}(cur_loc-2:cur_loc+WIN));
end

yaw_near_peak_R = zeros(1, length(locs_R));

WIN = 5; % each yaw is 0.01 s 

for i = 1:length(locs_R)
    
    cur_loc = locs_R(i);
    
    yaw_near_peak_R(i) = mean(yaw_all{1}(cur_loc-2:cur_loc+WIN));
end

f = figure;
hold on;
CIRCLE_SIZE = 20;
scatter(width_L,yaw_near_peak_L, CIRCLE_SIZE, 'b', 'filled');
scatter(width_R,yaw_near_peak_R, CIRCLE_SIZE, 'r', 'filled');

xlabel('Width of FR burst');
ylabel('Yaw (deg/s)');

legend({'A2 Left', 'A2 Right'});

saveas( f, [analysis_path '/A2_psth_width_right_and_left_vs_yaw_scatter.fig'] );
saveas( f, [analysis_path '/A2_psth_width_right_and_left_vs_yaw_scatter.png'] );




















