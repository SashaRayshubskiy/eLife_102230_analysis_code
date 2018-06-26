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

% A1 dual patching
% dirname = '180603_gfp_3G_ss731_dual_07';
%dirname = '180614_gfp_3G_ss731_dual_09';
% cell_1_label = 'A1 Left Vm';
% cell_2_label = 'A1 Right Vm';

% A2 dual patching 
dirname = '180410_gfp_3G_ss730_dual_08';
%dirname = '180430_gfp_3G_ss730_dual_12';   
cell_1_label = 'A2 Left Vm';
cell_2_label = 'A2 Right Vm';

% A1/A2 dual patching on the left
%dirname = '170816_2xGFP_ss731_75C10_01';
%dirname = '180503_2xGFP_ss731_75C10_03';
%dirname = '180503_2xGFP_ss731_75C10_02';
% cell_1_label = 'A2 Left Vm';
% cell_2_label = 'A1 Left Vm';
% cell_1_label = 'A1 Left Vm';
% cell_2_label = 'A2 Left Vm';

analysis_path = [working_dir dirname '/analysis/'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

%directories_to_analyze = { { '170816_2xGFP_ss731_75C10_01', [0], 11.5*150 } };
directories_to_analyze = { { '180410_gfp_3G_ss730_dual_08', [0], 11.5*300 } };
%directories_to_analyze = { { '180430_gfp_3G_ss730_dual_12', [0], 11.5*70 } };
%directories_to_analyze = { { '180603_gfp_3G_ss731_dual_07', [0], 11.5*140 } };
%directories_to_analyze = { { '180614_gfp_3G_ss731_dual_09', [0], 11.5*125 } }; % Good up to 270
%directories_to_analyze = { { '180503_2xGFP_ss731_75C10_02', [0], 11.5*125 } };
%directories_to_analyze = { { '180503_2xGFP_ss731_75C10_03', [0], 11.5*215 } };

[t_all, t_vel_all, yaw_all, fwd_all, ephys_all_A, ephys_all_B] = load_LAL_DN_data( working_dir, directories_to_analyze, ephys_SR, ball_SR, 1);

idx = 1;


%%

SAMPLING_RATE = 10000;
TRIAL_CNT_FACTOR = 4;
dt = SAMPLING_RATE*11.5*TRIAL_CNT_FACTOR;

voltage_A_r = reshape(ephys_all_A{1}, [dt (length(ephys_all_A{1}) / dt)]);
baseline_avg = squeeze(mean(voltage_A_r));
voltageA_meansub = voltage_A_r - repmat(baseline_avg, [dt, 1]);
voltageA = reshape(voltageA_meansub, [1 length(ephys_all_A{1})]);

voltage_B_r = reshape(ephys_all_B{1}, [dt (length(ephys_all_B{1}) / dt)]);
baseline_avg = squeeze(mean(voltage_B_r));
voltageB_meansub = voltage_B_r - repmat(baseline_avg, [dt, 1]);
voltageB = reshape(voltageB_meansub, [1 length(ephys_all_B{1})]);
        
FILT_FACTOR = 0.04;
VmFilt_A2_L = medfilt1( voltageA, FILT_FACTOR * ephys_SR, 'truncate' );
VmFilt_A2_R = medfilt1( voltageB, FILT_FACTOR * ephys_SR, 'truncate' );

START = 1;
END   = 10000*10;

figure;
hold on;
plot( t_all{1}(START:END), voltageA(START:END), 'b' );
plot( t_all{1}(START:END), VmFilt_A2_L(START:END), 'g' );

FILT_FACTOR = 0.1;

fwd_filt = medfilt1( fwd_all{1}, FILT_FACTOR * ball_SR, 'truncate' );
accel  = diff(fwd_filt);
accel(end+1) = 0;

dT = 1.0 / ball_SR;
accel = accel ./ dT;

accel_filt = medfilt1( accel, 0.05 * ball_SR, 'truncate' );

%% Calculate delta Vm

WINDOW_TRIAL_CNT = 10;
BASELINE_WINDOW = ephys_SR * 11.5 * WINDOW_TRIAL_CNT;

ephys_all_A_dVm = zeros(1,length(ephys_all_A{1}));
ephys_all_B_dVm = zeros(1,length(ephys_all_B{1}));

init_base_A = mean(ephys_all_B{1}(1:BASELINE_WINDOW));
init_base_B = mean(ephys_all_B{1}(1:BASELINE_WINDOW));

for i=1:length(ephys_all_A{1})
    
    if( i < (length(ephys_all_A{1})-BASELINE_WINDOW) )        
        ephys_all_A_dVm(i) = ephys_all_A{1}(i) - mean(ephys_all_A{1}(i:i+BASELINE_WINDOW));
        ephys_all_B_dVm(i) = ephys_all_A{1}(i) - mean(ephys_all_B{1}(i:i+BASELINE_WINDOW));    
    else
        ephys_all_A_dVm(i) = ephys_all_A{1}(i) - mean(ephys_all_A{1}(i-BASELINE_WINDOW:i));
        ephys_all_B_dVm(i) = ephys_all_A{1}(i) - mean(ephys_all_B{1}(i-BASELINE_WINDOW:i));            
    end
end

FILT_FACTOR = 0.04;

VmFilt_A2_L = medfilt1( ephys_all_A_dVm, FILT_FACTOR * ephys_SR, 'truncate' );
VmFilt_A2_R = medfilt1( ephys_all_B_dVm, FILT_FACTOR * ephys_SR, 'truncate' );

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

FILT_FACTOR = 0.1;

fwd_filt = medfilt1( fwd_all{1}, FILT_FACTOR * ball_SR, 'truncate' );
accel  = diff(fwd_filt);
accel(end+1) = 0;

dT = 1.0 / ball_SR;
accel = accel ./ dT;

accel_filt = medfilt1( accel, 0.05 * ball_SR, 'truncate' );

START = 1;
END   = 100*10;

figure;
yyaxis left;
hold on;
plot( t_vel_all{1}(START:END), fwd_all{1}(START:END), 'b' );
plot( t_vel_all{1}(START:END), fwd_filt(START:END), 'g' );

yyaxis right;
hold on;
plot( t_vel_all{1}(START:END), accel(START:END), 'r' );
plot( t_vel_all{1}(START:END), accel_filt(START:END), 'k' );


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

FILT_FACTOR = 0.1;

fwd_filt = medfilt1( fwd_all{1}, FILT_FACTOR * ball_SR, 'truncate' );
accel  = diff(fwd_filt);
accel(end+1) = 0;

dT = 1.0 / ball_SR;
accel = accel ./ dT;

accel_filt = medfilt1( accel, 0.05 * ball_SR, 'truncate' );

START = 1;
END   = 100*10;

figure;
yyaxis left;
hold on;
plot( t_vel_all{1}(START:END), fwd_all{1}(START:END), 'b' );
plot( t_vel_all{1}(START:END), fwd_filt(START:END), 'g' );

yyaxis right;
hold on;
plot( t_vel_all{1}(START:END), accel(START:END), 'r' );
plot( t_vel_all{1}(START:END), accel_filt(START:END), 'k' );


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

%% Plot L, R Vm overlayed with yaw

BIN_SIZE = 0.050; % s
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

t_down = squeeze(mean(reshape(t_all{1}, [DT_EPHYS, length(t_all{1})/DT_EPHYS]), 1));
A2_Vm_L_down = squeeze(mean(reshape(VmFilt_A2_L, [ DT_EPHYS, length(VmFilt_A2_L)/DT_EPHYS ] ),1));
A2_Vm_R_down = squeeze(mean(reshape(VmFilt_A2_R, [ DT_EPHYS, length(VmFilt_A2_R)/DT_EPHYS ] ),1));

yaw_all_down = squeeze(mean(reshape(yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));
fwd_all_down = squeeze(mean(reshape(fwd_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));

yaw_colors = zeros( length(A2_Vm_L_down), 3 );
fwd_colors = zeros( length(A2_Vm_L_down), 3 );

SHIFT_FACTOR = 3;
RBC_SIZE = 300;
rbc = jet(RBC_SIZE);
YAW_MAX = 500;
YAW_RANGE = [-1.0*YAW_MAX YAW_MAX];

FWD_MIN = -5;
FWD_MAX = 20;
FWD_RANGE = [FWD_MIN FWD_MAX];

for i = 1:length( A2_Vm_L_down ) - SHIFT_FACTOR
    cur_yaw = yaw_all_down( i + SHIFT_FACTOR );
    cur_fwd = fwd_all_down( i + SHIFT_FACTOR );
        
    if(cur_yaw < YAW_RANGE(1))
        cur_rbc_yaw_index = 1;
    elseif( cur_yaw > YAW_RANGE(2))
        cur_rbc_yaw_index = RBC_SIZE;    
    else
        cur_rbc_yaw_index = ceil(interp1( YAW_RANGE, [1,RBC_SIZE], cur_yaw ));
    end
    
    if(cur_fwd < FWD_RANGE(1))
        cur_rbc_fwd_index = 1;
    elseif( cur_fwd > FWD_RANGE(2))
        cur_rbc_fwd_index = RBC_SIZE;    
    else
        cur_rbc_fwd_index = ceil(interp1( FWD_RANGE, [1,RBC_SIZE], cur_fwd ));
    end
    
    yaw_colors(i,:) = rbc( cur_rbc_yaw_index, : );    
    fwd_colors(i,:) = rbc( cur_rbc_fwd_index, : );    
end

%% Show scatter plots
MARKER_SIZE = 2;
f = figure;

A2_Vm_RANGE = [-10 10];

subplot(1,2,1);
scatter( A2_Vm_L_down, A2_Vm_R_down, MARKER_SIZE, yaw_colors );
colormap('jet');
axis image;
xlim( A2_Vm_RANGE );
ylim( A2_Vm_RANGE );
xlabel('A2 Left Vm');
ylabel('A2 Right Vm');
caxis( YAW_RANGE );
colorbar;
title('Yaw overlay');
grid on;
refline(1,0);

subplot(1,2,2);
hold on;
scatter( A2_Vm_L_down, A2_Vm_R_down, MARKER_SIZE, fwd_colors );
colormap('jet');
axis image;
xlim( A2_Vm_RANGE );
ylim( A2_Vm_RANGE );
xlabel('A2 Left Vm');
ylabel('A2 Right Vm');
caxis( FWD_RANGE );
colorbar;
title('Fwd overlay');
grid on;
refline(1,0);

saveas(f, [analysis_path '/A2_left_vs_right_with_yaw_and_fwd_overlay.fig']);
saveas(f, [analysis_path '/A2_left_vs_right_with_yaw_and_fwd_overlay.png']);

%% Convert fwd vel to fwd accel
FILT_FACTOR = 0.1;

fwd_filt = medfilt1( fwd_all{1}, FILT_FACTOR * ball_SR, 'truncate' );
accel  = diff(fwd_filt);
accel(end+1) = 0;

dT = 1.0 / ball_SR;
accel = accel ./ dT;

accel_filt = medfilt1( accel, 0.05 * ball_SR, 'truncate' );

START = 1;
END   = 100*10;

figure;
yyaxis left;
hold on;
plot( t_vel_all{1}(START:END), fwd_all{1}(START:END), 'b' );
plot( t_vel_all{1}(START:END), fwd_filt(START:END), 'g' );

yyaxis right;
hold on;
plot( t_vel_all{1}(START:END), accel(START:END), 'r' );
plot( t_vel_all{1}(START:END), accel_filt(START:END), 'k' );



%% Bin the scatter plots
SHIFT_START = 0;
SHIFT_FINISH = 40;
EXPLORATION_COUNT = SHIFT_FINISH-SHIFT_START+1;

% LA_START = 5;
% LA_FINISH = 60;
% EXPLORATION_COUNT = LA_FINISH - LA_START + 1;

BIN_SIZE = 0.020; % s
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

t_down = squeeze(mean(reshape(t_all{1}, [DT_EPHYS, length(t_all{1})/DT_EPHYS]), 1));
A2_Vm_L_down = squeeze(mean(reshape(VmFilt_A2_L, [ DT_EPHYS, length(VmFilt_A2_L)/DT_EPHYS ] ),1));
A2_Vm_R_down = squeeze(mean(reshape(VmFilt_A2_R, [ DT_EPHYS, length(VmFilt_A2_R)/DT_EPHYS ] ),1));

yaw_all_down = squeeze(mean(reshape(yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));
% fwd_all_down = squeeze(mean(reshape(fwd_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));
fwd_all_down = squeeze(mean(reshape(fwd_filt, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));
accel_down = squeeze(mean(reshape(accel_filt, [ DT_YAW, length( accel_filt )/DT_YAW ]),1));

COMPUTE_MEAN_ACCEL = 0;
EXAMINE_LOOKAHEAD = 0;
EXAMINE_SHIFT_FACTOR = 1;
COMPUTE_FWD_ACCEL = 0;
LOOK_BACK = 5;
%LOOK_AHEAD = 35;

if( COMPUTE_FWD_ACCEL == 1 )
    % Single yaw point and accelration for fwd
    binned_type = 'yaw_SP_fwd_accel';
else
    binned_type = 'single_point';
end

VM_BIN_SIZE = 0.25;
VM_MAX = 10;
VM_RANGE = [ -1.0*VM_MAX VM_MAX ];

VM_BIN_EDGES = [ VM_RANGE(1) : VM_BIN_SIZE : VM_RANGE(2) ];
BIN_COUNT = length(VM_BIN_EDGES);

yaw_binned_final = zeros( BIN_COUNT, BIN_COUNT, EXPLORATION_COUNT );
fwd_binned_final = zeros( BIN_COUNT, BIN_COUNT, EXPLORATION_COUNT );

SHIFT_FACTOR = 8;    
sf = SHIFT_FACTOR;
LA_OPTIMAL = 23;

BLACKOUT = 999999;
LA = LA_OPTIMAL;

for sf = 8:8
%for sf = SHIFT_START:SHIFT_FINISH
%for sf = SHIFT_FACTOR:SHIFT_FACTOR
%for LA = LA_START:LA_FINISH
%for LA = LA_OPTIMAL:LA_OPTIMAL
    SHIFT_FACTOR = sf;

    t_accel = [ -1.0*LOOK_BACK*BIN_SIZE : BIN_SIZE : LA*BIN_SIZE ];
    yaw_binned_tmp = cell( BIN_COUNT, BIN_COUNT );
    fwd_binned_tmp = cell( BIN_COUNT, BIN_COUNT );
    
    for i = BIN_COUNT
        for j = BIN_COUNT
            yaw_binned_tmp{i,j} = [];
            fwd_binned_tmp{i,j} = [];
        end
    end
    
    %for i = 1:length( A2_Vm_L_down ) - SHIFT_FACTOR
    for i = (LOOK_BACK+1):( length( A2_Vm_L_down ) - SHIFT_FACTOR - LA-1)
        
        cur_yaw   = yaw_all_down( i + SHIFT_FACTOR );
        cur_fwd   = fwd_all_down( i + SHIFT_FACTOR );
        cur_accel = accel_down( i + SHIFT_FACTOR );
                
        
        cur_Vm_L = A2_Vm_L_down( i );
        cur_Vm_R = A2_Vm_R_down( i );
        
        cur_Vm_L_index = map_range( VM_RANGE, [1,BIN_COUNT], cur_Vm_L );
        cur_Vm_R_index = map_range( VM_RANGE, [1,BIN_COUNT], cur_Vm_R );
        
        % disp(['Current indecies: [ ' num2str(cur_Vm_L_index) ',' num2str(cur_Vm_R_index) ']']);
        
        yaw_binned_tmp{ cur_Vm_L_index, cur_Vm_R_index }( end+1 ) = cur_yaw;
        
        if( COMPUTE_FWD_ACCEL == 1 )
            if( COMPUTE_MEAN_ACCEL == 1 )                
                fwd_binned_tmp{ cur_Vm_L_index, cur_Vm_R_index }( end+1, : ) = fwd_all_down(i+SHIFT_FACTOR-LOOK_BACK:i+SHIFT_FACTOR+LA);
            else
                fwd_binned_tmp{ cur_Vm_L_index, cur_Vm_R_index }( end+1, : ) = cur_accel;
            end
        else
            fwd_binned_tmp{ cur_Vm_L_index, cur_Vm_R_index }( end+1 ) = cur_fwd;
        end
    end
    
    if( EXAMINE_LOOKAHEAD == 1)
        
        cur_la_index = LA-(LA_START-1);
        
        for i = 1:BIN_COUNT
            for j = 1:BIN_COUNT
                yaw_binned_final( i, j, cur_la_index ) = mean( yaw_binned_tmp{i,j} );
                                
                if( COMPUTE_FWD_ACCEL == 1)
                    if(  length(fwd_binned_tmp{i,j}) == 0 )
                        fwd_binned_final( i, j, cur_la_index ) = BLACKOUT;
                    else
                        fwd_binned_final( i, j, cur_la_index ) = get_accel( t_accel, fwd_binned_tmp{i,j} );
                    end
                else
                    fwd_binned_final( i, j, cur_la_index ) = mean( fwd_binned_tmp{i,j} );
                end
            end
        end
    elseif( EXAMINE_SHIFT_FACTOR == 1 )
        for i = 1:BIN_COUNT
            for j = 1:BIN_COUNT
                
                cur_bin_cnt = size(yaw_binned_tmp{i,j},2);
                

                BIN_CNT_THRESHOLD = 0;
                if( cur_bin_cnt < BIN_CNT_THRESHOLD  )
                    fwd_binned_final( i, j, sf+1 ) = BLACKOUT;
                    yaw_binned_final( i, j, sf+1 ) = BLACKOUT;
                else
                    fwd_binned_final( i, j, sf+1 ) = mean( fwd_binned_tmp{i,j} );                    
                    yaw_binned_final( i, j, sf+1 ) = mean( yaw_binned_tmp{i,j} );                    
                end
                
%                 if( COMPUTE_FWD_ACCEL == 1)
%                     if(  length(fwd_binned_tmp{i,j}) == 0 )
%                         fwd_binned_final( i, j, sf+1 ) = BLACKOUT;
%                     else
%                         if( COMPUTE_MEAN_ACCEL == 1 )         
%                             fwd_binned_final( i, j, sf+1 ) = get_accel( t_accel, fwd_binned_tmp{i,j} );
%                         else
%                             fwd_binned_final( i, j, sf+1 ) = mean( fwd_binned_tmp{i,j} ); % this is the instantaneous accel value
%                         end
%                     end
%                 else
%                     fwd_binned_final( i, j, sf+1 ) = mean( fwd_binned_tmp{i,j} );
%                 end
            end
        end
    end

%yaw_binned_final( find( isnan(yaw_binned_final) ) ) = BLACKOUT;
%fwd_binned_final( find( isnan(fwd_binned_final) ) ) = BLACKOUT;

yaw_binned_final( find( yaw_binned_final == BLACKOUT ) ) = NaN;
fwd_binned_final( find( fwd_binned_final == BLACKOUT ) ) = NaN;

% Display above
    
    SHIFT_TIME = SHIFT_FACTOR * BIN_SIZE;
    GAUSS_SIGMA = 0.25;
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1);
    hold on;
    %imagesc( VM_BIN_EDGES, VM_BIN_EDGES, imgaussfilt(squeeze( yaw_binned_final(:,:,sf+1)),GAUSS_SIGMA ));
    %imagesc( VM_BIN_EDGES, VM_BIN_EDGES, squeeze( yaw_binned_final(:,:,sf+1)));
    h = pcolor( VM_BIN_EDGES, VM_BIN_EDGES, squeeze( yaw_binned_final(:,:,sf+1)) );
    set(h, 'EdgeColor', 'none');
    xlabel( cell_2_label );
    ylabel( cell_1_label );
    axis image;
    colormap(jet);
    caxis([-300 300]);
    colorbar;
    % grid on;
    h = refline(1,0);
    h.LineWidth = 1.5;
    refline(-1,0);
    refline(0,0);
    plot([0 0 ],ylim);
    title(['Shift amount: ' num2str(SHIFT_TIME) 's -- look ahead: ' num2str(LA*BIN_SIZE) 's']);
    
    subplot(1,2,2);
    hold on;
    %imagesc( VM_BIN_EDGES, VM_BIN_EDGES, imgaussfilt(squeeze( fwd_binned_final(:,:,sf+1)), GAUSS_SIGMA) );
    %imagesc( VM_BIN_EDGES, VM_BIN_EDGES, squeeze( fwd_binned_final(:,:,cur_la_index)) );
    h = pcolor( VM_BIN_EDGES, VM_BIN_EDGES, squeeze( fwd_binned_final(:,:, sf+1)) );
    set(h, 'EdgeColor', 'none');
    
    xlabel( cell_2_label );
    ylabel( cell_1_label );
    axis image;
    colormap('jet');
    caxis([-15 15]);
    colorbar;
    % grid on;
    refline(1,0);
    refline(-1,0);
    refline(0,0);
    plot([0 0 ],ylim);
    
     saveas(f, [analysis_path '/A1_left_vs_right_Vm_with_yaw_and_fwd_overlay_' binned_type '_bin_size_' num2str(VM_BIN_SIZE) '_shift_time_' num2str(SHIFT_TIME) '_look_ahead_' num2str(LA) '.fig']);
     saveas(f, [analysis_path '/A1_left_vs_right_Vm_with_yaw_and_fwd_overlay_' binned_type '_bin_size_' num2str(VM_BIN_SIZE) '_shift_time_' num2str(SHIFT_TIME) '_look_ahead_' num2str(LA) '.png']);
%     saveas(f, [analysis_path '/A1_left_vs_right_Vm_with_yaw_and_fwd_overlay_' binned_type '_bin_size_' num2str(VM_BIN_SIZE) '_shift_time_' num2str(SHIFT_TIME) '_prediction_fwd.fig']);


%close(f);
end

%%

fwd_tmp = squeeze( fwd_binned_final(:,:,cur_la_index));

fwd_tmp( find( fwd_tmp < 0.5 ) ) = 0.0;

GAUSS_SIGMA = 1.5;

figure;
hold on;
imagesc( VM_BIN_EDGES, VM_BIN_EDGES, imgaussfilt(fwd_tmp, GAUSS_SIGMA) );
%imagesc( VM_BIN_EDGES, VM_BIN_EDGES, fwd_tmp );
xlabel('A2 Right Vm');
ylabel('A2 Left Vm');
axis image;
colormap('jet');
caxis([0 1]);
colorbar;
grid on;
refline(1,0);
refline(-1,0);
refline(0,0);
plot([0 0 ],ylim);

%% Play a movie of yaw and fwd relationships to L and R neuron as a function of shift factor

v = VideoWriter( [analysis_path '/A2_left_vs_right_Vm_with_yaw_and_accel_overlay_' binned_type '_bin_size_' num2str(VM_BIN_SIZE) '_shift_time.avi'] );
%v = VideoWriter( [analysis_path '/A2_left_vs_right_Vm_with_yaw_and_fwd_overlay_' binned_type '_bin_size_' num2str(VM_BIN_SIZE) '_shift_time.avi'] );
%v = VideoWriter( [analysis_path '/A2_left_vs_right_Vm_with_yaw_and_fwd_overlay_' binned_type '_bin_size_' num2str(VM_BIN_SIZE) '_look_ahead.avi'] );
v.FrameRate = 1;
open(v);
GAUSS_SIGMA = 1.75;

f = figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:size(yaw_binned_final,3)
    
    SHIFT_TIME = (i-1) * BIN_SIZE;
    
    subplot(1,2,1);
    hold on;
    %imagesc( VM_BIN_EDGES, VM_BIN_EDGES, imgaussfilt(squeeze( yaw_binned_final(:,:,sf+1)),GAUSS_SIGMA ));
    %imagesc( VM_BIN_EDGES, VM_BIN_EDGES, squeeze( yaw_binned_final(:,:,sf+1)));
    h = pcolor( VM_BIN_EDGES, VM_BIN_EDGES, squeeze( yaw_binned_final(:,:,i)) );
    set(h, 'EdgeColor', 'none');
    xlabel(cell_2_label);
    ylabel(cell_1_label);
    axis image;
    colormap(jet);
    caxis([-300 300]);
    colorbar;
    % grid on;
    h = refline(1,0);
    h.LineWidth = 1.5;
    refline(-1,0);
    refline(0,0);
    plot([0 0 ],ylim);
    title(['Shift amount: ' num2str(SHIFT_TIME) 's -- look ahead: ' num2str(LA*BIN_SIZE) 's'])
    
    subplot(1,2,2);
    hold on;
    %imagesc( VM_BIN_EDGES, VM_BIN_EDGES, imgaussfilt(squeeze( fwd_binned_final(:,:,sf+1)), GAUSS_SIGMA) );
    %imagesc( VM_BIN_EDGES, VM_BIN_EDGES, squeeze( fwd_binned_final(:,:,cur_la_index)) );
    h = pcolor( VM_BIN_EDGES, VM_BIN_EDGES, squeeze( fwd_binned_final(:,:, i)) );
    set(h, 'EdgeColor', 'none');
    
    xlabel(cell_2_label);
    ylabel(cell_1_label);
    axis image;
    colormap('jet');
    caxis([-15 15]);
    colorbar;
    % grid on;
    refline(1,0);
    refline(-1,0);
    refline(0,0);
    plot([0 0 ],ylim);
  
%     subplot(1,2,1);
%     colormap('jet');
%     hold on;
%     imagesc(VM_BIN_EDGES, VM_BIN_EDGES, squeeze( yaw_binned_final(:,:,i)) );
%     %imagesc( VM_BIN_EDGES, VM_BIN_EDGES, imgaussfilt(squeeze( yaw_binned_final(:,:,i)),GAUSS_SIGMA) );
%     xlabel('A2 Right Vm');
%     ylabel('A2 Left Vm');
%     axis image;
%     caxis([-300 300]);
%     colorbar;
%     grid on;
%     refline(1,0);
%     refline(-1,0);
%     refline(0,0);
%     plot([0 0 ],ylim);
%     title(['Shift amount: ' num2str(SHIFT_TIME) 's'])
%     
%     subplot(1,2,2);
%     hold on;
%     imagesc( VM_BIN_EDGES, VM_BIN_EDGES, squeeze( fwd_binned_final(:,:,i)) );
%     %imagesc( VM_BIN_EDGES, VM_BIN_EDGES, imgaussfilt(squeeze( fwd_binned_final(:,:,i)),GAUSS_SIGMA) );
%     xlabel('A2 Right Vm');
%     ylabel('A2 Left Vm');
%     axis image;
%     colormap('jet');
%     caxis([-8 8]);
%     colorbar;
%     grid on;
%     refline(1,0);
%     refline(-1,0);
%     refline(0,0);
%     plot([0 0 ],ylim);    
%     title(['Look ahead: ' num2str((LA_START+i-1)* BIN_SIZE) 's']);    
    
    frame = getframe(f);   
    writeVideo(v,frame);
    
    pause(0.1);
end

close(v);


%% Plot 2D histograms of yaw vs. fwd vel

f = figure;
c=[1 10 100];

% h = refline(1,0);
% h.LineWidth = 1.5;
% refline(-1,0);
% refline(0,0);
% plot([0 0 ],ylim);

[N,cen] = hist3( [yaw_all{1}, fwd_all{1}], 'NBins', [50 50] );

%hold on;
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
grid on;


if( USE_JET == 1 )
    saveas(f,[analysis_path '/A1_A2_yaw_vs_fwd_scatter_logcolor_v2_jet.fig']);
    saveas(f,[analysis_path '/A1_A2_yaw_vs_fwd_scatter_logcolor_v2_jet.png']);
else
    saveas(f,[analysis_path '/A1_A2_yaw_vs_fwd_scatter_logcolor_v2_gray.fig']);
    saveas(f,[analysis_path '/A1_A2_yaw_vs_fwd_scatter_logcolor_v2_gray.png']);
end














