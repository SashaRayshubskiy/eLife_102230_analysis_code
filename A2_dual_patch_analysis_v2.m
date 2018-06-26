%% Load all data from all the animals to be included in this analysis

clear all;
close all;

%working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/descending_neurons/';

working_dir = '/data/drive1/sasha/';

%%% WARNING: MAKE SURE EPHYS SAMPLE RATE IS CORRECT 10000 Hz
%%% WARNING: MAKE SURE EPHYS SAMPLE RATE IS CORRECT 10000 Hz
%%% WARNING: MAKE SURE EPHYS SAMPLE RATE IS CORRECT 10000 Hz

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

dirname = '180410_gfp_3G_ss730_dual_08';

analysis_path = [working_dir dirname '/analysis/'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

directories_to_analyze = { { '180410_gfp_3G_ss730_dual_08', [0], 11.5*300 } };

[t_all, t_vel_all, yaw_all, fwd_all, ephys_all_A, ephys_all_B] = load_LAL_DN_data( working_dir, directories_to_analyze, ephys_SR, ball_SR, 1);

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

linkaxes(ax, 'x')
%% Extract Vm (filter out the spikes)
FILT_FACTOR = 0.04;

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

%% Correlation analysis: Using diff as a proxy for correlation 
% Examine the difference of Vm

Vm_diff = VmFilt_A2_L - VmFilt_A2_R;

XLIM = 6;
NBINS = 100;

f = figure;
subplot(3,1,1)
hist(Vm_diff,NBINS);
xlim([-1*XLIM XLIM]);
title('A2 L/R Vm diff');

subplot(3,1,2)
hist(VmFilt_A2_L,NBINS);
xlim([-1*XLIM XLIM]);
title('A2 L Vm');

subplot(3,1,3)
hist(VmFilt_A2_R,NBINS);
xlim([-1*XLIM XLIM]);
title('A2 R Vm');
xlabel('Vm (mV)');

saveas(f, [analysis_path '/Vm_histograms.fig']);
saveas(f, [analysis_path '/Vm_histograms.png']);

%% Gauss unmix models
Vm_diff = VmFilt_A2_L - VmFilt_A2_R;

k = 2;
GMModel_diff = fitgmdist(Vm_diff, k, 'Options', statset('Display', 'final'));

GMModel_L = fitgmdist(VmFilt_A2_L, k, 'Options', statset('Display', 'final'));

GMModel_R = fitgmdist(VmFilt_A2_R, k, 'Options', statset('Display', 'final'));

%% Show Vm histograms for both neurons, overlayed with mixed gaussian fit

f = figure;

X = [-10:0.01:10];
my_norm = [];

subplot(3,1,1);
% Normalize the density to match the total area of the histogram
h = histogram(Vm_diff,500);
binwidth = h.BinWidth;
area = (length(Vm_diff) * binwidth) / 1.45;

X = [-10:0.01:10];
my_norm = [];
for i = 1:k
    my_norm(i,:) = area * normpdf(X, GMModel_diff.mu(i), GMModel_diff.Sigma(i));
end

xlim([-10 10]);
hold on;
for i = 1:k
    plot(X, squeeze(my_norm(i,:))); 
end
title('L-R Vm diff');

subplot(3,1,2)
% Normalize the density to match the total area of the histogram
h = histogram(VmFilt_A2_L,500);
binwidth = h.BinWidth;
area = (length(Vm_diff) * binwidth) / 1.45;

X = [-10:0.01:10];
my_norm = [];
for i = 1:k
    my_norm(i,:) = area * normpdf(X, GMModel_L.mu(i), GMModel_L.Sigma(i));
end

xlim([-10 10]);
hold on;
for i = 1:k
    plot(X, squeeze(my_norm(i,:))); 
end
title('L Vm');

subplot( 3, 1, 3 )
% Normalize the density to match the total area of the histogram
h = histogram( VmFilt_A2_R, 500 );
binwidth = h.BinWidth;
area = (length(Vm_diff) * binwidth) / 1.45;

X = [-10:0.01:10];
my_norm = [];
for i = 1:k
    my_norm(i,:) = area * normpdf(X, GMModel_R.mu(i), GMModel_R.Sigma(i));
end

xlim([-10 10]);
hold on;
for i = 1:k
    plot(X, squeeze(my_norm(i,:))); 
end
title('R Vm');

%% Correlation analysis: compare correlation and large Vm

BIN_SIZE = 0.050; % s
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

t_down = squeeze(mean(reshape(t_all{1}, [DT_EPHYS, length(t_all{1})/DT_EPHYS]), 1));
A2_Vm_L_down = squeeze(mean(reshape(VmFilt_A2_L, [ DT_EPHYS, length(VmFilt_A2_L)/DT_EPHYS ] ),1));
A2_Vm_R_down = squeeze(mean(reshape(VmFilt_A2_R, [ DT_EPHYS, length(VmFilt_A2_R)/DT_EPHYS ] ),1));

yaw_all_down = squeeze(mean(reshape(yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));
fwd_all_down = squeeze(mean(reshape(fwd_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));

% 5 Categories
% 1. Large negative diff
% 2. Large positive diff
% 3. Small diff and large postitive Vm for both L and R
% 4. Small diff and large negative Vm for both L and R
% 5. Small diff and small Vm for both L and R

CAT_COUNT = 5;
yaw_cat = cell(1,CAT_COUNT);
fwd_cat = cell(1,CAT_COUNT);
L_cell = cell(1,CAT_COUNT);
R_cell = cell(1,CAT_COUNT);
for i = 1:CAT_COUNT
    yaw_cat{i} = [];
    fwd_cat{i} = [];
    L_cell{i} = [];
    R_cell{i} = [];
end

VM_LARGE_THRESHOLD = 2.0;
VM_SMALL_THRESHOLD = 0.25;

VM_DIFF_LARGE_THRESHOLD = 2.0;
VM_DIFF_SMALL_THRESHOLD = 0.25;
LOOK_BACK = 5;
LOOK_AHEAD = 5;

SHIFT_FACTOR = 3;

L_cell_stat = [];
R_cell_stat = [];
yaw_stat = [];
fwd_stat = [];

L_cell_mv = [];
R_cell_mv = [];
yaw_mv = [];
fwd_mv = [];
STATIONARY_THRESHOLD = 0.25;
MOVING_THRESHOLD = 0.5;

for i = LOOK_BACK:length(t_down)-LOOK_AHEAD

    cur_A2_Vm_L = A2_Vm_L_down( i );
    cur_A2_Vm_R = A2_Vm_R_down( i );
    cur_Vm_diff = cur_A2_Vm_L - cur_A2_Vm_R;
    
    cur_cat = -1;
    
    if(cur_Vm_diff > VM_DIFF_LARGE_THRESHOLD)
        cur_cat = 1;
    elseif(cur_Vm_diff < -1.0*VM_DIFF_LARGE_THRESHOLD)
        cur_cat = 2;
    elseif( (cur_Vm_diff > -1.0*VM_DIFF_SMALL_THRESHOLD) && (cur_Vm_diff < VM_DIFF_SMALL_THRESHOLD))
        
        if( ( cur_A2_Vm_L > VM_LARGE_THRESHOLD ) && ( cur_A2_Vm_R > VM_LARGE_THRESHOLD ))
            cur_cat = 3;
        elseif( ( cur_A2_Vm_L < -1.0*VM_LARGE_THRESHOLD ) && ( cur_A2_Vm_R < -1.0*VM_LARGE_THRESHOLD ))
            cur_cat = 4;
        elseif( ( cur_A2_Vm_L > -1.0*VM_SMALL_THRESHOLD ) && ( cur_A2_Vm_L < VM_SMALL_THRESHOLD ) && ...
                ( cur_A2_Vm_R > -1.0*VM_SMALL_THRESHOLD ) && ( cur_A2_Vm_R < VM_SMALL_THRESHOLD ))            
            cur_cat = 5;
        end
    end        
    
    if( cur_cat == -1 )
        continue;
    end
        
    L_cell{cur_cat}(end+1,:)  = A2_Vm_L_down(i-LOOK_BACK:i+LOOK_AHEAD);
    R_cell{cur_cat}(end+1,:)  = A2_Vm_R_down(i-LOOK_BACK:i+LOOK_AHEAD);    

    i_shift = i + SHIFT_FACTOR;
    yaw_cat{cur_cat}(end+1,:) = yaw_all_down( i_shift-LOOK_BACK : i_shift+LOOK_AHEAD);
    fwd_cat{cur_cat}(end+1,:) = fwd_all_down( i_shift-LOOK_BACK : i_shift+LOOK_AHEAD);
    
    % Divide category 3 trials by stationary to moving vs. moving
    if( cur_cat == 3 )
        cur_fwd = fwd_all_down( i_shift-LOOK_BACK : i_shift+LOOK_AHEAD);
        
        if( (mean(cur_fwd(1:LOOK_BACK)) < STATIONARY_THRESHOLD ) && (mean(cur_fwd(LOOK_BACK:LOOK_AHEAD)) > MOVING_THRESHOLD ))
            L_cell_stat(end+1,:) = A2_Vm_L_down(i-LOOK_BACK:i+LOOK_AHEAD);
            R_cell_stat(end+1,:)  = A2_Vm_R_down(i-LOOK_BACK:i+LOOK_AHEAD);
            yaw_stat(end+1,:)  = yaw_all_down( i_shift-LOOK_BACK : i_shift+LOOK_AHEAD);
            fwd_stat(end+1,:)  = fwd_all_down( i_shift-LOOK_BACK : i_shift+LOOK_AHEAD);            
        elseif( (mean(cur_fwd(1:LOOK_BACK)) > MOVING_THRESHOLD ) && (mean(cur_fwd(LOOK_BACK:LOOK_AHEAD)) > MOVING_THRESHOLD ))
            L_cell_mv(end+1,:) = A2_Vm_L_down(i-LOOK_BACK:i+LOOK_AHEAD);
            R_cell_mv(end+1,:)  = A2_Vm_R_down(i-LOOK_BACK:i+LOOK_AHEAD);
            yaw_mv(end+1,:)  = yaw_all_down( i_shift-LOOK_BACK : i_shift+LOOK_AHEAD);
            fwd_mv(end+1,:)  = fwd_all_down( i_shift-LOOK_BACK : i_shift+LOOK_AHEAD);
        end
    end
end

f = figure('units','normalized','outerposition',[0 0 1 1])  
cur_subplot_base = 1;
t_event = [ -1.0*LOOK_BACK*BIN_SIZE : BIN_SIZE : LOOK_AHEAD*BIN_SIZE ];

category_text = { 'Large + diff', 'Large - diff', ... 
                  'Small diff and large + Vm, L & R', ...
                  'Small diff and large - Vm, L & R', ...
                  'Small diff and small Vm, L & R' };

for c = 1:CAT_COUNT
        
    if(c == 1)
        cur_subplots = [1 6 11 16];
    elseif(c == 2)
        cur_subplots = [2 7 12 17];
    elseif(c == 3)
        cur_subplots = [3 8 13 18];
    elseif(c == 4)
        cur_subplots = [4 9 14 19];
    elseif(c == 5)
        cur_subplots = [5 10 15 20];
    end
    
    subplot(4,5,cur_subplots(1));
    plot(t_event, mean( fwd_cat{c} ) );
    ylabel('Fwd vel (mm/s)');
    title([category_text{c} ]);
    ylim([0 10])
    legend(['Events: ' num2str(size(fwd_cat{c},1))]);

    subplot(4,5,cur_subplots(2));
    plot(t_event, mean( yaw_cat{c} ) );
    ylabel('Yaw vel (deg/s)');
    ylim([-150 150])
    
    subplot(4,5,cur_subplots(3));
    plot(t_event, mean( L_cell{c} ) );
    ylabel('Left cell Vm (mV)');
    ylim([-3.0 3.0]);
    
    subplot(4,5,cur_subplots(4));
    plot(t_event, mean( R_cell{c} ) );    
    ylabel('Right cell Vm (mV)');
    ylim([-3.0 3.0]);
    xlabel('Time (s)');
    
    cur_subplot_base = cur_subplot_base + 4;
end

saveas(f, [analysis_path '/left_right_by_category.fig']);
saveas(f, [analysis_path '/left_right_by_category.png']);

f = figure;
sub_idx = [1 3 5 7 2 4 6 8];
subplot(4,2,sub_idx(1));
plot(t_event, mean(fwd_stat));
ylabel('Fwd vel (mm/s)');
legend(['Events: ' num2str(size(fwd_stat,1))]);
title('Stationary to moving');

subplot(4,2,sub_idx(2));
plot(t_event, mean(yaw_stat));
ylabel('Yaw vel (deg/s)');

subplot(4,2,sub_idx(3));
plot(t_event, mean(L_cell_stat));
ylabel('Left cell (mV)');

subplot(4,2,sub_idx(4));
plot(t_event, mean(R_cell_stat));
ylabel('Right cell (mV)');
xlabel('Time (s)');

subplot(4,2,sub_idx(5));
plot(t_event, mean(fwd_mv));
ylabel('Fwd vel (mm/s)');
legend(['Events: ' num2str(size(fwd_mv,1))]);
title('Moving');

subplot(4,2,sub_idx(6));
plot(t_event, mean(yaw_mv));
ylabel('Yaw vel (deg/s)');

subplot(4,2,sub_idx(7));
plot(t_event, mean(L_cell_mv));
ylabel('Left cell (mV)');

subplot(4,2,sub_idx(8));
plot(t_event, mean(R_cell_mv));
ylabel('Right cell (mV)');
xlabel('Time (s)');

saveas(f, [analysis_path '/moving_vs_stationary_category.fig']);
saveas(f, [analysis_path '/moving_vs_stationary_category.png']);

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
psth_dt_samples = ephys_SR / ball_SR;
tic; A2_L_psth = calculate_psth_A2( t_all{1}, t_vel_all{1}, ephys_all_A{1}, ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples ); toc;
tic; A2_R_psth = calculate_psth_A2( t_all{1}, t_vel_all{1}, ephys_all_B{1}, ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples ); toc;

%% Test for binning artifact 

figure;
plot(A2_L_psth, 'o', 'MarkerSize', 3);

figure;
plot(A2_L_psth-A2_R_psth, 'o', 'MarkerSize', 3);


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

%A2_L_psth_base = A2_L_psth - mean(A2_L_psth);
%A2_R_psth_base = A2_R_psth - mean(A2_R_psth);

A2_L_psth_base = A2_L_psth;
A2_R_psth_base = A2_R_psth;

A2_L_psth_down = squeeze(mean(reshape( A2_L_psth_base, [ DT_YAW, length(A2_L_psth) / DT_YAW ] ), 1));
A2_R_psth_down = squeeze(mean(reshape( A2_R_psth_base, [ DT_YAW, length(A2_L_psth) / DT_YAW ] ), 1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( yaw_all{1} ) / DT_YAW ]), 1));
% fwd_all_down = squeeze(mean(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} ) / DT_YAW ]), 1));
fwd_all_down_std = squeeze(std(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} ) / DT_YAW ]), 0, 1));

% figure;
% hist(fwd_all_down_std, 1000);

cnt = 0;
A2_L_psth_choosen = zeros(1, length( A2_L_psth_down ));
A2_R_psth_choosen = zeros(1, length( A2_L_psth_down ));
yaw_choosen  = zeros(1, length( yaw_all_down ));

PSTH_BASELINE_CUTOFF = 1.0;
SHIFT_FACTOR = 3;
mcnt = 0;

WINDOW_SIZE = 5;
FWD_STD_THRESHOLD = 0.01;

for i = 1:length( A2_L_psth_down )-(SHIFT_FACTOR+1)-WINDOW_SIZE
    
    cur_A2_L_psth = A2_L_psth_down( i );
    cur_A2_R_psth = A2_R_psth_down( i );
    
    cur_fwd_std = fwd_all_down_std( i );
    
    if( abs(cur_fwd_std) < FWD_STD_THRESHOLD )
       continue; 
    end
    
%    if( (( cur_A2_L_psth <= -1.0*PSTH_BASELINE_CUTOFF ) || ( cur_A2_L_psth >= PSTH_BASELINE_CUTOFF )) && ...
%        (( cur_A2_R_psth <= -1.0*PSTH_BASELINE_CUTOFF ) || ( cur_A2_R_psth >= PSTH_BASELINE_CUTOFF )) )
        A2_L_psth_choosen( mcnt+1 ) = cur_A2_L_psth;
        A2_R_psth_choosen( mcnt+1 ) = cur_A2_R_psth;
        yaw_choosen( mcnt+1 ) = yaw_all_down(i+SHIFT_FACTOR);
        mcnt = mcnt + 1;
%    end        
end

MARKER_SIZE = 1;

f = figure;

subplot(2,2,1);
plot( A2_L_psth_choosen(1:mcnt), yaw_choosen(1:mcnt), 'o', 'MarkerSize', MARKER_SIZE );
xlabel('A2 L PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
title('Left A2 vs. yaw');
grid on;

subplot(2,2,2);
plot( A2_R_psth_choosen(1:mcnt), yaw_choosen(1:mcnt), 'o', 'MarkerSize', MARKER_SIZE );
xlabel('A2 R PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
title('Right A2 vs. yaw');
grid on;

subplot(2,2,3);
A2_L_R_psth_choosen_diff = A2_L_psth_choosen(1:mcnt) - A2_R_psth_choosen(1:mcnt);

plot( A2_L_R_psth_choosen_diff, yaw_choosen(1:mcnt), 'o', 'MarkerSize', MARKER_SIZE );
xlabel('A2 L-R PSTH diff (spikes/s)');
ylabel('Yaw (deg/s)');
title('diff(Left,Right) A2 vs. yaw');
grid on;

subplot(2,2,4);
%plot( A2_L_psth_choosen(1:cnt) ./ A2_R_psth_choosen(1:cnt), yaw_choosen(1:cnt), 'o', 'MarkerSize', 3 );
semilogx( A2_L_psth_choosen(1:mcnt) ./ A2_R_psth_choosen(1:mcnt), yaw_choosen(1:mcnt), 'o', 'MarkerSize', MARKER_SIZE );
xlabel('A2 L/R PSTH ratio');
ylabel('Yaw (deg/s)');
title('div(Left,Right) A2 vs. yaw');
grid on;

saveas(f,[analysis_path '/A2_right_left_PSTH_diff_div_vs_yaw_scatter.fig']);
saveas(f,[analysis_path '/A2_right_left_PSTH_diff_div_vs_yaw_scatter.png']);

%%
f = figure;
scatter3( A2_L_psth_choosen(1:mcnt), yaw_choosen(1:mcnt), A2_L_R_psth_choosen_diff );
xlabel('A2 L psth (spike/s)');
ylabel('yaw vel (deg/s)');
zlabel('A2 L-R psth diff(spike/s)');

saveas(f,[analysis_path '/A2_left_PSTH_diff_vs_yaw_scatter3.fig']);
saveas(f,[analysis_path '/A2_left_PSTH_diff_vs_yaw_scatter3.png']);


%% WARNING: This depends on previous block. 
% Plot A2 L/R diff vs. yaw with A2 Left PSTH overlay

f = figure;

A2_small_diff = zeros(1, length(A2_L_R_psth_choosen_diff) );
A2_small_diff_yaw = zeros(1, length(A2_L_R_psth_choosen_diff) );
A2_small_diff_L_psth_colors = zeros( length(A2_L_R_psth_choosen_diff), 3 );

DIFF_CUTOFF = 100;
RBC_SIZE = 40;
rbc = jet( RBC_SIZE );

min_PSTH = min( min(A2_L_psth_choosen), min(A2_R_psth_choosen) );
max_PSTH = max( max(A2_L_psth_choosen), max(A2_R_psth_choosen) );

cnt = 0;
for i = 1:length(A2_L_R_psth_choosen_diff)
    
    cur_diff = A2_L_R_psth_choosen_diff( i );
    
    %if( A2_L_psth_choosen( i ) > 3 )
    if( abs(cur_diff) < DIFF_CUTOFF )        
        A2_small_diff(cnt+1)        = cur_diff;
        A2_small_diff_yaw(cnt+1)    = yaw_choosen( i );                        
        
        PSTH_L_R_avg = (A2_L_psth_choosen( i ) + A2_R_psth_choosen( i ))/2.0;

        % map [-3.56 60] range to [1 40];
        %
        %cur_color = rbc( rbc_index, : );
        
        if( PSTH_L_R_avg < 0.1 )
            cur_color = rgb('Magenta');
            %rbc_index = ceil(interp1([min_PSTH, 0], [1,RBC_SIZE], PSTH_L_R_avg)); 
            %cur_color = rbc( rbc_index, : );
        else
            cur_color = rgb('DimGray');
        end

%         if( PSTH_L_R_avg < -1.3 )
%             cur_color = rgb('Magenta');
%         else
%             cur_color = rgb('DimGray');
%         end
        
        A2_small_diff_L_psth_colors(cnt+1,:) = cur_color;
               
        cnt = cnt + 1;
    end
end

scatter( A2_small_diff, A2_small_diff_yaw, 1, A2_small_diff_L_psth_colors);

if 0
colormap('jet')
h = colorbar;
caxis([0 RBC_SIZE]);
ticks = get(h, 'Ticks');
cmap_labels = cell(1, length(ticks));
for i = 1:length( ticks )
    x = interp1([0,RBC_SIZE], [min_PSTH, 0], ticks(i));
    
    cmap_labels{i} = num2str(x,2);
end
set(h, 'TickLabels', cmap_labels);
ylabel(h, 'A2 avg(L,R) PSTH below baseline (spikes/s)');
end

grid on;
xlabel('A2 L-R diff PSTH (spikes/s)');
ylabel('Yaw (deg/s)');

Y = A2_small_diff_yaw;
X = A2_small_diff;
[fobj, gof] = fit(X',Y', 'poly1');

title(['A2 left, right diff vs yaw with L,R PSTH = zero overlay (R^2 = ' num2str(gof.rsquare) ')']);


saveas(f, [analysis_path '/A2_left_right_diff_vs_yaw_with_left_PSTH_overlay_v2.fig']);
saveas(f, [analysis_path '/A2_left_right_diff_vs_yaw_with_left_PSTH_overlay_v2.png']);

%% Calculate avg yaw from a 2D histogram of the above plot

PSTH_BINS = 50;
YAW_BINS = 100;

[N, c] = hist3( [A2_small_diff', A2_small_diff_yaw'], 'NBins', [PSTH_BINS YAW_BINS], 'CDataMode', 'auto', 'FaceColor', 'interp' );

psth_N_weighted = N .* repmat(c{1}', [1 size(N,2)]);
yaw_N_weighted = N .* repmat(c{2}, [size(N,1) 1]);

weight_sum_yaw = sum(N,1);
weight_sum_psth = sum(N,2);

A1_psth_all_wavg_yaw_x = squeeze(sum(psth_N_weighted, 1)) ./ weight_sum_yaw;
yaw_all_wavg_psth_x = sum(yaw_N_weighted, 2) ./ weight_sum_psth;

f = figure;
% subplot(2,1,1);
% plot(c{2}, A1_psth_all_wavg_yaw_x);
% ylim([0 100]);
% xlim([-1000 1000]);
% xlabel('Yaw (deg/s)');
% ylabel('A2 L-R PSTH diff (spikes/s)');
% 
subplot(1,1,1);
plot(c{1}, yaw_all_wavg_psth_x');
%xlim([0 100]);
ylim([-1000 1000]);
grid on;
ylabel('Yaw (deg/s)');
xlabel('A2 L-R PSTH diff (spikes/s)');
title(['Avg yaw for each A2 L-R diff bin (n bins= ' num2str(PSTH_BINS) ')']);

saveas(f, [analysis_path '/A2_left_right_diff_vs_yaw_weighted_avg.fig']);
saveas(f, [analysis_path '/A2_left_right_diff_vs_yaw_weighted_avg.png']);

%% Compare various yaw distributions
% a. when both cells are not firing
% b. when both cells are firing a small amount (but not zero)
% c. all of yaw

f = figure('units','normalized','outerposition',[0 0 1 1]);

NBINS = 1000;
XMIN = -1000;
XMAX = 1000;
subplot(3,1,1);
hist(yaw_choosen(1:mcnt), NBINS);
title(['All yaw, mean: ' num2str(mean(yaw_choosen(1:mcnt))) '  std: ' num2str(std(yaw_choosen(1:mcnt))) ]);
xlim([XMIN XMAX]);
ylim([0 800]);

yaw_when_both_A2_PSTH_is_zero = zeros(1,length(A2_L_R_psth_choosen_diff));
yaw_when_both_A2_PSTH_is_zero_cnt = 0;

yaw_when_both_A2_PSTH_is_small = zeros(1,length(A2_L_R_psth_choosen_diff));
yaw_when_both_A2_PSTH_is_small_cnt = 0;

ZERO_CUTOFF = 0.1;
SMALL_CUTOFF = 10;

for i = 1:length(A2_L_R_psth_choosen_diff)
    
    PSTH_L_R_avg = (A2_L_psth_choosen( i ) + A2_R_psth_choosen( i ))/2.0;

    cur_diff = A2_L_R_psth_choosen_diff( i );
    
    if( (A2_L_psth_choosen( i ) < ZERO_CUTOFF) && (A2_R_psth_choosen( i ) < ZERO_CUTOFF) )

        yaw_when_both_A2_PSTH_is_zero( yaw_when_both_A2_PSTH_is_zero_cnt + 1 ) = yaw_choosen( i );
        yaw_when_both_A2_PSTH_is_zero_cnt = yaw_when_both_A2_PSTH_is_zero_cnt + 1;
    end
    
%     if( ((A2_L_psth_choosen( i ) > ZERO_CUTOFF) && (A2_L_psth_choosen( i ) < SMALL_CUTOFF)) && ...
%             ((A2_R_psth_choosen( i ) > ZERO_CUTOFF) && (A2_R_psth_choosen( i ) < SMALL_CUTOFF)) )
        
    if( ( abs(cur_diff) > 1 ) &&  ( abs(cur_diff) < SMALL_CUTOFF ) )
        yaw_when_both_A2_PSTH_is_small( yaw_when_both_A2_PSTH_is_small_cnt + 1 ) = yaw_choosen( i );
        yaw_when_both_A2_PSTH_is_small_cnt = yaw_when_both_A2_PSTH_is_small_cnt + 1;
    end
end

subplot(3,1,2);
hist(yaw_when_both_A2_PSTH_is_zero(1:yaw_when_both_A2_PSTH_is_zero_cnt), NBINS);
title(['Zero PSTH in both yaw, mean: ' num2str(mean(yaw_when_both_A2_PSTH_is_zero)) '  std: ' num2str(std(yaw_when_both_A2_PSTH_is_zero)) ]);
xlim([XMIN XMAX]);
ylim([0 550]);

subplot(3,1,3);
hist(yaw_when_both_A2_PSTH_is_small(1:yaw_when_both_A2_PSTH_is_small_cnt), NBINS);
title(['1 < L-R < ' num2str(SMALL_CUTOFF) ' spike/s, yaw mean: ' num2str(mean(yaw_when_both_A2_PSTH_is_small)) '  std: ' num2str(std(yaw_when_both_A2_PSTH_is_small)) ]);
xlim([XMIN XMAX]);
ylim([0 80]);

saveas(f, [analysis_path '/yaw_distribution_comparison.fig']);
saveas(f, [analysis_path '/yaw_distribution_comparison.png']);

[h,p] = ttest2(yaw_choosen(1:mcnt), yaw_when_both_A2_PSTH_is_small(1:yaw_when_both_A2_PSTH_is_small_cnt))
[h,p] = ttest2(yaw_choosen(1:mcnt), yaw_when_both_A2_PSTH_is_zero(1:yaw_when_both_A2_PSTH_is_zero_cnt))
[h,p] = ttest2(yaw_when_both_A2_PSTH_is_small(1:yaw_when_both_A2_PSTH_is_small_cnt), yaw_when_both_A2_PSTH_is_zero(1:yaw_when_both_A2_PSTH_is_zero_cnt))

%% Plot avg between both neurons and yaw

PSTH_L_R_avg = (A2_L_psth_choosen + A2_R_psth_choosen)/2.0;

f = figure;

plot( PSTH_L_R_avg, yaw_choosen, 'o', 'MarkerSize', 3 );
xlabel('A2 left/right avg (spikes/s)'); 
ylabel('Yaw (deg/s)');

saveas(f, [analysis_path 'A2_left_right_avg_vs_yaw_scatter.fig']);
saveas(f, [analysis_path 'A2_left_right_avg_vs_yaw_scatter.png']);

%% Plot avg between both neurons and yaw

PSTH_L_R_avg = (A2_L_psth_choosen + A2_R_psth_choosen)/2.0;

% MIN_YAW = min(yaw_choosen);
% MAX_YAW = max(yaw_choosen);
MIN_YAW = -500;
MAX_YAW =  500;

RBC_SIZE = 100;
rbc = jet( RBC_SIZE );

yaw_colors = zeros(length(A2_L_R_psth_choosen_diff),3);
for i=1:length(yaw_colors)

    cur_yaw = yaw_choosen(i);
    
    if(cur_yaw <= MIN_YAW)
        rbc_index = 1;
    elseif(cur_yaw >= MAX_YAW)
        rbc_index = RBC_SIZE;
    else
        rbc_index = ceil(interp1([MIN_YAW, MAX_YAW], [1,RBC_SIZE], cur_yaw)); 
    end
    
    yaw_colors(i,:) = rbc(rbc_index,:);
end

f = figure;

scatter( A2_L_R_psth_choosen_diff, PSTH_L_R_avg(1:length(A2_L_R_psth_choosen_diff)), 3, yaw_colors );
ylabel('A2 left/right avg (spikes/s)'); 
xlabel('A2 left/right diff (spikes/s)');
colorbar;
colormap jet;
caxis([MIN_YAW MAX_YAW]);
grid on;

saveas(f, [analysis_path 'A2_left_right_avg_vs_diff_with_yaw_overlay_scatter.fig']);
saveas(f, [analysis_path 'A2_left_right_avg_vs_diff_with_yaw_overlay_scatter.png']);


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

%A2_L_psth_base = A2_L_psth - mean(A2_L_psth);
%A2_R_psth_base = A2_R_psth - mean(A2_R_psth);
A2_L_psth_base = A2_L_psth;
A2_R_psth_base = A2_R_psth;

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

YAW_CUTOFF = 150;
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
    
    % Uncommend to remove points from plot based on YAW_CUTOFF
    if( abs(cur_yaw) < YAW_CUTOFF ) 
        A2_yaw_colors_choosen( cnt+1, : ) = rgb('Magenta');
    else
        A2_yaw_colors_choosen( cnt+1, : ) = rgb('DimGray');
    end
    
    A2_L_choosen( cnt+1 ) = A2_L_psth_down(i);
    A2_R_choosen( cnt+1 ) = A2_R_psth_down(i);
    yaw_choosen( cnt+1 ) = cur_yaw;
    cnt = cnt + 1;

% Uncommend to remove points from plot based on YAW_CUTOFF
%     if( (cur_yaw < -1 * YAW_CUTOFF) || (cur_yaw > YAW_CUTOFF) )
%     % if( (cur_yaw > -1 * YAW_CUTOFF) && (cur_yaw < YAW_CUTOFF) )
%         A2_L_choosen( cnt+1 ) = A2_L_psth_down(i);
%         A2_R_choosen( cnt+1 ) = A2_R_psth_down(i);
%         yaw_choosen( cnt+1 ) = cur_yaw;
%         A2_yaw_colors_choosen( cnt+1, : ) = cur_clr;
%         cnt = cnt + 1;
%     end    
end

f = figure;
rbc = jet(RBC_SIZE);
%scatter(A2_L_psth_down, A2_R_psth_down, 4, A2_yaw_colors);
scatter( A2_L_choosen, A2_R_choosen, 2, A2_yaw_colors_choosen);
%scatter( A2_L_choosen-A2_R_choosen, yaw_choosen, 3 );
%colormap('jet')
% h = colorbar;
% ylabel(h, 'Yaw vel (deg/s)');
% caxis([-1000 1000]);
xlabel('A2 Left PSTH (spikes/s)');
ylabel('A2 Right PSTH (spikes/s)');
grid on;
axis tight;
title(['abs(yaw velocity) < ' num2str(YAW_CUTOFF) ' = magenta, otherwise = gray']);

saveas(f,[analysis_path '/A2_right_left_PSTH_scatter_with_yaw_overlay_v2.fig']);
saveas(f,[analysis_path '/A2_right_left_PSTH_scatter_with_yaw_overlay_v2.png']);

%%
% Plot a 2D histrogram for A2 L PSTH vs. A2 R PSTH
f = figure;

hist3([A2_L_choosen', A2_R_choosen'], 'NBins', [50 50], 'CDataMode','auto','FaceColor','interp' );
xlabel('Left A2 PSTH (spikes/s)');
ylabel('Right A2 PSTH (spikes/s)');
zlabel('count');

USE_JET = 1;
if( USE_JET == 1 )
    colormap( jet ); 
else
    colormap( gray ); 
end

caxis([0 70]);
view(2)
xlim([0 90]);
ylim([0 90]);
axis tight;
h = colorbar('FontSize',11);
ylabel(h, 'Counts');

if( USE_JET == 1 )
    saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay_2D_hist_v2_jet.fig']);
    saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay_2D_hist_v2_jet.png']);
else
    saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay_2D_hist_v2_gray.fig']);
    saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay_2D_hist_v2_gray.png']);
end


%%
% Plot a 2D histrogram for A2 L PSTH vs. A2 R PSTH using log color scale
f = figure;

if 0
hist3([A2_L_choosen', A2_R_choosen'], 'NBins', [50 50], 'CDataMode','auto','FaceColor','interp' );
xlabel('Left A2 PSTH (spikes/s)');
ylabel('Right A2 PSTH (spikes/s)');
zlabel('count');
colormap jet;
colorbar;
caxis([0 70]);
view(2)
axis tight;
end

c=[1 10 100 300];


[N,cen] = hist3( [A2_L_choosen', A2_R_choosen'], 'NBins', [50 50] );

imagesc(cen{1}([1 end]), cen{2}([1 end]), log(N));

USE_JET = 0;
if( USE_JET == 1 )
    colormap( jet ); 
else
    colormap( gray ); 
end

caxis(log([c(1) c(length(c))]));
h = colorbar('FontSize',11,'YTick',log(c),'YTickLabel',c);
ylabel(h, 'Counts');

xlabel('Left A2 PSTH (spikes/s)');
ylabel('Right A2 PSTH (spikes/s)');
zlabel('count');
set(gca,'YDir','normal')
xlim([0 85]);
ylim([0 85]);

if(USE_JET == 1)
    saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay_2D_hist_logcolor_v2_jet.fig']);
    saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay_2D_hist_logcolor_v2_jet.png']);
else
    saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay_2D_hist_logcolor_v2_gray.fig']);
    saveas(f,[analysis_path '/A2_right_left_diff_PSTH_scatter_with_yaw_overlay_2D_hist_logcolor_v2_gray.png']);
end

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


%% Plot R^2 as a functin of shift factor
BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;


A2_L_psth_base = A2_L_psth;
A2_R_psth_base = A2_R_psth;

A2_L_psth_down = squeeze(mean(reshape( A2_L_psth_base, [ DT_YAW, length(A2_L_psth) / DT_YAW ] ), 1));
A2_R_psth_down = squeeze(mean(reshape( A2_R_psth_base, [ DT_YAW, length(A2_L_psth) / DT_YAW ] ), 1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( yaw_all{1} ) / DT_YAW ]), 1));
fwd_all_down = squeeze(mean(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} ) / DT_YAW ]), 1));
fwd_all_down_std = squeeze(std(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} ) / DT_YAW ]), 0, 1));

SHIFT_MAX = 50;

FWD_STD_THRESHOLD = 0.01;

rsqr_per_shift = zeros( 1, SHIFT_MAX );

for ss = 1:SHIFT_MAX
    
    cnt = 0;
    A2_L_psth_choosen = zeros(1, length( A2_L_psth_down ));
    A2_R_psth_choosen = zeros(1, length( A2_L_psth_down ));
    yaw_choosen  = zeros(1, length( yaw_all_down ));
    fwd_choosen  = zeros(1, length( yaw_all_down ));
    
    mcnt = 0;
    
    for i = 1:length( A2_L_psth_down )-(ss+1)
        
        cur_A2_L_psth = A2_L_psth_down( i );
        cur_A2_R_psth = A2_R_psth_down( i );
        
        cur_fwd_std = fwd_all_down_std( i );
        
        if( abs(cur_fwd_std) < FWD_STD_THRESHOLD )
            continue;
        end
        
        A2_L_psth_choosen( mcnt+1 ) = cur_A2_L_psth;
        A2_R_psth_choosen( mcnt+1 ) = cur_A2_R_psth;
        yaw_choosen( mcnt+1 ) = yaw_all_down(i+ss);
        fwd_choosen( mcnt+1 ) = fwd_all_down(i+ss);
        mcnt = mcnt + 1;
    end
    
    Y = yaw_choosen(1:mcnt);
    X = A2_L_psth_choosen(1:mcnt) - A2_R_psth_choosen(1:mcnt);
    [fobj, gof] = fit(X',Y', 'poly1');

    rsqr_per_shift(ss) = gof.rsquare;
end

f = figure;
plot( rsqr_per_shift );
xlabel(['Shift amount (' num2str(BIN_SIZE) ' sec per unit)']);
ylabel('R squared');
title('R^2 value vs. shift factor');

saveas(f, [analysis_path '/rsquare_vs_shift_factor.fig']);
saveas(f, [analysis_path '/rsquare_vs_shift_factor.png']);

%% 

GREATER_PSTH_CUTOFF = 15;
%LESSER_PSTH_CUTOFF = 10;
LEFT_PSTH_CUTOFF = 40;
RIGHT_PSTH_CUTOFF = LEFT_PSTH_CUTOFF;
OTHER_PSTH_CUTOFF = 5;

yaw_in_greater_cutoff = zeros(1, length(yaw_choosen));
fwd_in_greater_cutoff = zeros(1, length(fwd_choosen));

% yaw_in_lesser_cutoff = zeros(1, length(yaw_choosen));
% fwd_in_lesser_cutoff = zeros(1, length(fwd_choosen));

yaw_in_left_cutoff = zeros(1, length(yaw_choosen));
fwd_in_left_cutoff = zeros(1, length(fwd_choosen));

yaw_in_right_cutoff = zeros(1, length(yaw_choosen));
fwd_in_right_cutoff = zeros(1, length(fwd_choosen));

cur_g_idx = 1;
cur_l_idx = 1;
cur_left_idx = 1;
cur_right_idx = 1;
for i = 1:length(A2_L_psth_choosen)
    
    cur_L_psth = A2_L_psth_choosen(i);
    cur_R_psth = A2_R_psth_choosen(i);
    
    if( (cur_L_psth > GREATER_PSTH_CUTOFF) && (cur_R_psth > GREATER_PSTH_CUTOFF) )
        yaw_in_greater_cutoff(cur_g_idx) = yaw_choosen(i);
        fwd_in_greater_cutoff(cur_g_idx) = fwd_choosen(i);  
        cur_g_idx = cur_g_idx + 1;
    end
    
    if( (cur_L_psth > LEFT_PSTH_CUTOFF) && (cur_R_psth < OTHER_PSTH_CUTOFF) )
        yaw_in_left_cutoff(cur_left_idx) = yaw_choosen(i);
        fwd_in_left_cutoff(cur_left_idx) = fwd_choosen(i);  
        cur_left_idx = cur_left_idx + 1;
    end
    
    if( (cur_R_psth > RIGHT_PSTH_CUTOFF) && (cur_L_psth < OTHER_PSTH_CUTOFF))
        yaw_in_right_cutoff(cur_right_idx) = yaw_choosen(i);
        fwd_in_right_cutoff(cur_right_idx) = fwd_choosen(i);  
        cur_right_idx = cur_right_idx + 1;
    end
end

NBINS = 60;
f = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,3,1);
hist( yaw_in_greater_cutoff(1:cur_g_idx), NBINS );
xlabel('Yaw vel (deg/s)');
ylabel('Count');
title(['A2 left and right PSTH > ' num2str(GREATER_PSTH_CUTOFF)]);
legend(['Mean: ' num2str(mean(yaw_in_greater_cutoff(1:cur_g_idx))) ' std: ' num2str(std(yaw_in_greater_cutoff(1:cur_g_idx))) ]);


subplot(2,3,2);
hist( yaw_in_left_cutoff(1:cur_left_idx), NBINS );
xlabel('Yaw vel (deg/s)');
ylabel('Count');
title(['A2 left PSTH > ' num2str( LEFT_PSTH_CUTOFF ) ' right PSTH < ' num2str(OTHER_PSTH_CUTOFF)]);
legend(['Mean: ' num2str(mean(yaw_in_left_cutoff(1:cur_left_idx))) ' std: ' num2str(std(yaw_in_left_cutoff(1:cur_left_idx))) ]);


subplot(2,3,3);
hist( yaw_in_right_cutoff(1:cur_right_idx), NBINS );
xlabel('Yaw vel (deg/s)');
ylabel('Count');
title(['A2 right PSTH > ' num2str( RIGHT_PSTH_CUTOFF ) ' left PSTH < ' num2str(OTHER_PSTH_CUTOFF)]);
legend(['Mean: ' num2str(mean(yaw_in_right_cutoff(1:cur_right_idx))) ' std: ' num2str(std(yaw_in_right_cutoff(1:cur_right_idx))) ]);


subplot(2,3,4);
hist( fwd_in_greater_cutoff(1:cur_g_idx), NBINS );
xlabel('Fwd vel (mm/s)');
ylabel('Count');
legend(['Mean: ' num2str(mean(fwd_in_greater_cutoff(1:cur_g_idx))) ' std: ' num2str(std(fwd_in_greater_cutoff(1:cur_g_idx))) ]);


subplot(2,3,5);
hist( fwd_in_left_cutoff(1:cur_left_idx), NBINS );
xlabel('Fwd vel (mm/s)');
ylabel('Count');
legend(['Mean: ' num2str(mean(fwd_in_left_cutoff(1:cur_left_idx))) ' std: ' num2str(std(fwd_in_left_cutoff(1:cur_left_idx))) ]);


subplot(2,3,6);
hist( fwd_in_right_cutoff(1:cur_right_idx), NBINS );
xlabel('Fwd vel (mm/s)');
ylabel('Count');
legend(['Mean: ' num2str(mean(fwd_in_right_cutoff(1:cur_right_idx))) ' std: ' num2str(std(fwd_in_right_cutoff(1:cur_right_idx))) ]);


saveas(f, [analysis_path '/yaw_fwd_hist_in_A2_left_right_mutual_PSTH_gcutoff_' num2str(GREATER_PSTH_CUTOFF) '_lcutoff_' num2str(LESSER_PSTH_CUTOFF) '.fig']);
saveas(f, [analysis_path '/yaw_fwd_hist_in_A2_left_right_mutual_PSTH_gcutoff_' num2str(GREATER_PSTH_CUTOFF) '_lcutoff_' num2str(LESSER_PSTH_CUTOFF) '.png']);














