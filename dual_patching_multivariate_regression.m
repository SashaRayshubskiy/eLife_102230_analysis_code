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

%% Save data for python

fwd_vel = fwd_all{1};
yaw_vel = yaw_all{1};

save([analysis_path dirname '_data_for_python.mat'], 'fwd_vel', 'yaw_vel', 'ephys_all_A', 'ephys_all_B', '-v7');

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

%% Plot auto correlation of yaw and fwd vel

figure;

subplot(2,1,1)
autocorr(fwd_all{1});
%plot(xcorr(fwd_all{1}));
title('Fwd vel');

subplot(2,1,2)
autocorr(yaw_all{1});
%plot(xcorr(yaw_all{1}));
title('Yaw vel');


%% Multi-variate regression 
% Set data

cell_1_data = VmFilt_A2_L;
cell_2_data = VmFilt_A2_R;


BIN_SIZE = 0.020; % s
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

t_down = squeeze(mean(reshape(t_all{1}, [DT_EPHYS, length(t_all{1})/DT_EPHYS]), 1));

cell_1_down = squeeze(mean(reshape(cell_1_data, [ DT_EPHYS, length(cell_2_data)/DT_EPHYS ] ),1));
cell_2_down = squeeze(mean(reshape(cell_2_data, [ DT_EPHYS, length(cell_2_data)/DT_EPHYS ] ),1));

yaw_all_down = squeeze(mean(reshape(yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));
fwd_all_down = squeeze(mean(reshape(fwd_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));
fwd_all_std = squeeze(std(reshape(fwd_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));

% Exclude trials when the fly is standing
cell_1_down_temp = cell_1_down;
cell_2_down_temp = cell_2_down;
yaw_all_down_temp = yaw_all_down;
fwd_all_down_temp = fwd_all_down;

cell_1_down = [];
cell_2_down = [];
yaw_all_down = [];
fwd_all_down = [];

FWD_STD_THRESHOLD = 0.01;

i=1;
while( i < length(fwd_all_down_temp) )
    
    cur_fwd_std = fwd_all_std(i);
    
    if( cur_fwd_std > FWD_STD_THRESHOLD )
        cell_1_down(end+1) = cell_1_down_temp(i);
        cell_2_down(end+1) = cell_2_down_temp(i);
        yaw_all_down(end+1) = yaw_all_down_temp(i);
        fwd_all_down(end+1) = fwd_all_down_temp(i);
    end
    
    i = i + 1;
end

%% Plot sanity check variables

figure;

subplot(2,2,1);
scatter(cell_1_down, yaw_all_down);
grid on;

subplot(2,2,2);
scatter(cell_1_down, fwd_all_down);
grid on;

subplot(2,2,3);
scatter(cell_2_down, yaw_all_down);
grid on;

subplot(2,2,4);
scatter(cell_2_down, fwd_all_down);
grid on;

%% Multi-variate regression using mvregress

 % X = n x p
 % Y = n x d
 % beta = p x d

%X = [cell_1_down', cell_2_down'];
Y = [fwd_all_down', yaw_all_down' ];
X = [ cell_1_down', cell_2_down'];
%Y = [ fwd_all_down' ];
[ beta, Sigma, E, CovB, logL ] = mvregress( X, Y );

figure;
plot(E)

% calculate_r_squared_adjusted( beta, Y, X );
% 1 = fwd, 2 = yaw
SStotal = [];
SSresid = [];
R2adjusted = [];

n = size(Y,1);
p = size(X,2);

for j = 1:size(Y,2)
    SSresid(j) = sum(squeeze(E(:,j)).^2);
    %SStotal(j) = (length(Y(:,j))-1) * var(Y(:,j));
    
    SStotal(j) = sum((Y(:,j) - mean(squeeze(Y(:,j)))).^2);
    
    %R2adjusted(j) = 1 - ((SSresid(j) / SStotal(j))*((n-1)/(n-p-1)));
    R2adjusted(j) = 1 - (SSresid(j) / SStotal(j));
end

figure;
bar(R2adjusted);







































accel_filt = medfilt1( accel, 0.05 * ball_SR, 'truncate' );