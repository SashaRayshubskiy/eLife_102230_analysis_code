%% Load all data from all the animals to be included in this analysis

clear all;
close all;

%working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/descending_neurons/';

working_dir = '/data/drive1/sasha/';

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

dirname = '170816_2xGFP_ss731_75C10_01';

analysis_path = [working_dir '/summary_analysis/A1_A2_dual_patch_analysis/'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

directories_to_analyze = { { '170816_2xGFP_ss731_75C10_01', [0,1], 2*11.5*150 } };
% directories_to_analyze = { { '180503_2xGFP_ss731_75C10_02', [0], 11.5*126 } };
% directories_to_analyze = { { '180503_2xGFP_ss731_75C10_03', [0], 11.5*216 } };

[t_all, t_vel_all, yaw_all, fwd_all, ephys_all_A, ephys_all_B] = load_LAL_DN_data( working_dir, directories_to_analyze, ephys_SR, ball_SR, 1 );

idx = 1;

if( strcmp(directories_to_analyze{1}{1}, '180503_2xGFP_ss731_75C10_03') == 1 )
    % For this dataset A1 and A2 are switched
    tmpA        = ephys_all_A;
    ephys_all_A = ephys_all_B;
    ephys_all_B = tmpA;
end

%% Show histgram of fwd vel
figure;
hist(fwd_all{1}, 1000);

%% Convert physiology to PSTH and process Vm
SPIKE_THRESHOLD_LAL_DN = 0.25;
psth_dt_samples = ephys_SR/ball_SR;
tic; A2_psth = calculate_psth_A2( t_all{1}, t_vel_all{1}, ephys_all_A{1}, ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples ); toc;

SPIKE_THRESHOLD_LAL_DN = 2.0;
tic; A1_psth = calculate_psth_A1( t_all{1}, t_vel_all{1}, ephys_all_B{1}, ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples ); toc;

BIN_SIZE = 0.050; % ms
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

% Downsample 
yaw_t_down = squeeze(mean(reshape(t_vel_all{1}, [DT_YAW, length(t_vel_all{1})/DT_YAW]),1));
A2_psth_down = squeeze(mean(reshape(A2_psth, [DT_YAW, length(A2_psth)/DT_YAW]),1));
A1_psth_down = squeeze(mean(reshape(A1_psth, [DT_YAW, length(A1_psth)/DT_YAW]),1));

% Process Vm
% Remove spikes first
FILT_FACTOR = 0.04;
VmFilt_A2_L_tmp = medfilt1( ephys_all_A{1}, FILT_FACTOR * ephys_SR, 'truncate' );
VmFilt_A2_R_tmp = medfilt1( ephys_all_B{1}, FILT_FACTOR * ephys_SR, 'truncate' );

% computer dVm
TRIAL_CNT_FACTOR = 6; 
dt = ephys_SR*11.5*TRIAL_CNT_FACTOR;

voltage_A_r = reshape(VmFilt_A2_L_tmp, [dt (length(ephys_all_A{1}) / dt)]);
baseline_avg = squeeze(mean(voltage_A_r));
voltageA_meansub = voltage_A_r - repmat(baseline_avg, [dt, 1]);
A2_Vm_b = reshape(voltageA_meansub, [1 length(ephys_all_A{1})]);

voltage_B_r = reshape(VmFilt_A2_R_tmp, [dt (length(ephys_all_B{1}) / dt)]);
baseline_avg = squeeze(mean(voltage_B_r));
voltageB_meansub = voltage_B_r - repmat(baseline_avg, [dt, 1]);
A1_Vm_b = reshape(voltageB_meansub, [1 length(ephys_all_B{1})]);       

% Must get Vm on the same timebase as ball data
BIN_SIZE = 0.01; % s
DT_EPHYS = ephys_SR * BIN_SIZE;

A2_Vm = squeeze(mean(reshape( A2_Vm_b, [ DT_EPHYS, length(A2_Vm_b)/DT_EPHYS ] ),1));
A1_Vm = squeeze(mean(reshape( A1_Vm_b, [ DT_EPHYS, length(A1_Vm_b)/DT_EPHYS ] ),1));


%% Isolate transitions from stationary to moving

SAMPLE_CNT = length(A2_psth);
%SAMPLE_CNT = 10000;

FROM_STATIONARY_THRESHOLD = 2.0;
FROM_MOVING_THRESHOLD = 0.001;
STATIONARY_STATE = 11;
MOVING_STATE = 12;

FWD_LOOK_AHEAD = 50;
WINDOW_SIZE = 5;

fwd_stationary_to_moving = [];
yaw_stationary_to_moving = [];
A1_FR_stationary_to_moving = [];
A2_FR_stationary_to_moving = [];
A1_Vm_stationary_to_moving = [];
A2_Vm_stationary_to_moving = [];

LOOK_BACK = 70;
LOOK_AHEAD = 50;
   
cur_fwd_std = std(fwd_all{1}(1:1+FWD_LOOK_AHEAD));
cur_fwd_mean = mean(fwd_all{1}(1:1+FWD_LOOK_AHEAD));
if(cur_fwd_std < FROM_MOVING_THRESHOLD)
    cur_state = STATIONARY_STATE;
else
    cur_state = MOVING_STATE;
end

stationary_points = [];
moving_points = [];

transition_points = [];
STATIONARY_STATE_CNT_THRESHOLD = 50;
cur_stationary_state_cnt = 0;
i = LOOK_BACK+1;
FWD_MOVEMENT_THRESHOLD = 0.1;

f = figure;
while( i < (SAMPLE_CNT-LOOK_AHEAD) )
    
    cur_fwd = fwd_all{1}(i:i+FWD_LOOK_AHEAD);
    cur_fwd_std = std(cur_fwd);

    if(cur_state == STATIONARY_STATE)
        
        stationary_points(end+1) = i;
        
        if( cur_fwd_std > FROM_STATIONARY_THRESHOLD ) 
            
            % Mark as an event only if the fly is standing for some time
            % And fwd velocity is positive
            if( (cur_stationary_state_cnt > STATIONARY_STATE_CNT_THRESHOLD ) && (mean(cur_fwd) > 0.05 ))
                
                % Find the exact point of change
                cur_change_interval = cur_fwd;
                
                hold on;
                plot(cur_change_interval);
                
                prev_fwd_val = cur_change_interval(1);
                center = i + 40;
%                 for j = 2:length( cur_change_interval )
%                     cur_fwd_val = cur_change_interval(j);
%                     
%                     if((cur_fwd_val - prev_fwd_val) > FWD_MOVEMENT_THRESHOLD )
%                         center = center + j -1;
%                         break;
%                     end
%                     prev_fwd_val = cur_fwd_val;
%                 end
                                
                fwd_stationary_to_moving(end+1,:)    = fwd_all{1}( center-LOOK_BACK : center+LOOK_AHEAD );
                yaw_stationary_to_moving(end+1,:)    = yaw_all{1}( center-LOOK_BACK : center+LOOK_AHEAD ) - mean(yaw_all{1}( center-LOOK_BACK : center));
                A1_FR_stationary_to_moving(end+1,:)  = A1_psth( center-LOOK_BACK : center+LOOK_AHEAD );
                A2_FR_stationary_to_moving(end+1,:)  = A2_psth( center-LOOK_BACK : center+LOOK_AHEAD );
                A1_Vm_stationary_to_moving(end+1,:)  = A1_Vm( center-LOOK_BACK : center+LOOK_AHEAD );
                A2_Vm_stationary_to_moving(end+1,:)  = A2_Vm( center-LOOK_BACK : center+LOOK_AHEAD );
                transition_points(end+1) = i;
            end
            i = i + FWD_LOOK_AHEAD;            
            cur_state = MOVING_STATE;
        else
            i = i + 1;
        end
        
        cur_stationary_state_cnt = cur_stationary_state_cnt + 1;
        
    elseif(cur_state == MOVING_STATE) 
        
        moving_points(end+1) = i;

        if( cur_fwd_std < FROM_MOVING_THRESHOLD )
            cur_state = STATIONARY_STATE;
            cur_stationary_state_cnt = 0;
        end
        i = i + 1;        
    else
        disp(['ERROR: state ' num2str(cur_state) ' not recognized']);
    end
end

if 0
f = figure;

BEGIN = 1;
END = SAMPLE_CNT;

my_ax(1) = subplot(3,1,1);
hold on;
plot(t_vel_all{1}(BEGIN:END)', fwd_all{1}(BEGIN:END)');
plot(t_vel_all{1}(stationary_points)', fwd_all{1}(stationary_points)', 'x', 'color', 'g');
plot(t_vel_all{1}(moving_points)', fwd_all{1}(moving_points)', 'x', 'color', 'r');
plot(t_vel_all{1}(transition_points)', fwd_all{1}(transition_points)', 'x', 'color', rgb('Magenta'), 'MarkerSize', 6, 'LineWidth', 3);

my_ax(2) = subplot(3,1,2);
plot(t_vel_all{1}(BEGIN:END)', yaw_all{1}(BEGIN:END)');

my_ax(3) = subplot(3,1,3);
hold on
plot( t_vel_all{1}(BEGIN:END), A2_psth(BEGIN:END)' );
plot( t_vel_all{1}(BEGIN:END), A1_psth(BEGIN:END)', 'g' );

linkaxes(my_ax, 'x');
axis tight;
end

%% Requires previous statement
f = figure('units','normalized','outerposition',[0 0 1 1])

event_cnt = size(fwd_stationary_to_moving, 1);

t_event = [ -1.0*LOOK_BACK*0.01 : 0.01 : LOOK_AHEAD*0.01 ];

left_yaw = [];
left_fwd = [];
left_A1 = [];
left_A2 = [];
left_A1_Vm = [];
left_A2_Vm = [];

right_yaw = [];
right_fwd = [];
right_A1 = [];
right_A2 = [];
right_A1_Vm = [];
right_A2_Vm = [];

for i = 1:event_cnt

    cur_yaw = mean(yaw_stationary_to_moving(i,LOOK_BACK+1:LOOK_BACK+20));
    if( cur_yaw > 0 )
        cur_clr = rgb('Green');
        left_yaw(end+1, :) = yaw_stationary_to_moving(i,:);
        left_fwd(end+1, :) = fwd_stationary_to_moving(i,:);
        left_A1(end+1, :)  = A1_FR_stationary_to_moving(i,:);
        left_A2(end+1, :)  = A2_FR_stationary_to_moving(i,:);
        left_A1_Vm(end+1, :)  = A1_Vm_stationary_to_moving(i,:);
        left_A2_Vm(end+1, :)  = A2_Vm_stationary_to_moving(i,:);
    else
        cur_clr = rgb('Red');
        right_yaw(end+1, :) = yaw_stationary_to_moving(i,:);
        right_fwd(end+1, :) = fwd_stationary_to_moving(i,:);
        right_A1(end+1, :)  = A1_FR_stationary_to_moving(i,:);
        right_A2(end+1, :)  = A2_FR_stationary_to_moving(i,:);
        right_A1_Vm(end+1, :)  = A1_Vm_stationary_to_moving(i,:);
        right_A2_Vm(end+1, :)  = A2_Vm_stationary_to_moving(i,:);
    end
    
    if 0
    subplot(4,1,1);
    hold on;
    plot( t_event, fwd_stationary_to_moving(i,:), 'color', cur_clr );
    
    subplot(4,1,2);
    hold on;
    plot( t_event, yaw_stationary_to_moving(i,:), 'color', cur_clr );
    
    subplot(4,1,3);
    hold on;
    plot( t_event, A1_FR_stationary_to_moving(i,:), 'color', cur_clr );
        
    subplot(4,1,4);
    hold on;
    plot( t_event, A2_FR_stationary_to_moving(i,:), 'color', cur_clr );    
    end
end

subplot(4,1,1);
hold on;

avg_fwd_left = squeeze(mean(left_fwd));
sem_fwd_left = get_sem( left_fwd, 1 );

avg_fwd_right = squeeze(mean(right_fwd));
sem_fwd_right = get_sem( right_fwd, 1 );

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_fwd_left+sem_fwd_left) fliplr((avg_fwd_left-sem_fwd_left))], ...
        rgb('Salmon'));
set(fh, 'EdgeColor', 'None');

pl_0 = plot( t_event, avg_fwd_left, 'color', rgb('Red'));

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_fwd_right+sem_fwd_right) fliplr((avg_fwd_right-sem_fwd_right))], ...
        rgb('PaleGreen'));
set(fh, 'EdgeColor', 'None');

pl_1 = plot( t_event, avg_fwd_right, 'color', rgb('Green'));
ylabel('Fwd vel (mm/s)');
axis tight;
legend([pl_0, pl_1], ['Left turn events (' num2str(size(left_fwd,1)) ')'], ['Right turn events (' num2str(size(right_fwd,1)) ')'], 'Location', 'northwest');


subplot(4,1,2);
hold on;
avg_yaw_left = squeeze(mean(left_yaw));
sem_yaw_left = get_sem( left_yaw, 1 );

avg_yaw_right = squeeze(mean(right_yaw));
sem_yaw_right = get_sem( right_yaw, 1 );

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_yaw_left+sem_yaw_left) fliplr((avg_yaw_left-sem_yaw_left))], ...
        rgb('Salmon'));
set(fh, 'EdgeColor', 'None');

pl_0 = plot( t_event, avg_yaw_left, 'color', rgb('Red'));

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_yaw_right+sem_yaw_right) fliplr((avg_yaw_right-sem_yaw_right))], ...
        rgb('PaleGreen'));
set(fh, 'EdgeColor', 'None');

pl_1 = plot( t_event, avg_yaw_right, 'color', rgb('Green'));
ylabel('Yaw vel (deg/s)');
axis tight;

SHOW_VM = 1;
if (SHOW_VM == 1)
subplot(4,1,3);
hold on;

avg_A1_right = squeeze(mean(right_A1_Vm));
sem_A1_right = get_sem( right_A1_Vm, 1 );

avg_A2_right = squeeze(mean(right_A2_Vm));
sem_A2_right = get_sem( right_A2_Vm, 1 );

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_A1_right+sem_A1_right) fliplr((avg_A1_right-sem_A1_right))], ...
        rgb('Gray'));
set(fh, 'EdgeColor', 'None');

pl1 = plot( t_event, avg_A1_right, 'color', rgb('Black'));

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_A2_right+sem_A2_right) fliplr((avg_A2_right-sem_A2_right))], ...
        rgb('Violet'));
set(fh, 'EdgeColor', 'None');

pl2 = plot( t_event, avg_A2_right, 'color', rgb('Magenta'));
            
ylabel('Vm (mV)');
title( 'Right turning' );
% axis tight;
xlim( [t_event(1) t_event(end)] );

subplot(4,1,4);
hold on;

avg_A1_left = squeeze(mean(left_A1_Vm));
sem_A1_left = get_sem( left_A1_Vm, 1 );

avg_A2_left = squeeze(mean(left_A2_Vm));
sem_A2_left = get_sem( left_A2_Vm, 1 );

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_A1_left+sem_A1_left) fliplr((avg_A1_left-sem_A1_left))], ...
        rgb('Gray'));
set(fh, 'EdgeColor', 'None');

pl1 = plot( t_event, avg_A1_left, 'color', rgb('Black'));

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_A2_left+sem_A2_left) fliplr((avg_A2_left-sem_A2_left))], ...
        rgb('Violet'));
set(fh, 'EdgeColor', 'None');

pl2 = plot( t_event, avg_A2_left, 'color', rgb('Magenta'));
ylabel('Vm (mV)');
title( 'Left turning' );
% axis tight;
xlabel('Time (s)');
xlim( [t_event(1) t_event(end)] );

legend( [pl1, pl2], ['A1'], ['A2'], 'Location', 'northwest' );

dirpath = directories_to_analyze{1}{1};

saveas(f, [ analysis_path '/' dirpath '_A1_A2_stationary_to_moving_analysis_Vm.fig']);
saveas(f, [ analysis_path '/' dirpath '_A1_A2_stationary_to_moving_analysis_Vm.png']);

else
subplot(4,1,3);
hold on;

avg_A1_right = squeeze(mean(right_A1));
sem_A1_right = get_sem( right_A1, 1 );

avg_A2_right = squeeze(mean(right_A2));
sem_A2_right = get_sem( right_A2, 1 );

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_A1_right+sem_A1_right) fliplr((avg_A1_right-sem_A1_right))], ...
        rgb('Gray'));
set(fh, 'EdgeColor', 'None');

pl1 = plot( t_event, avg_A1_right, 'color', rgb('Black'));

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_A2_right+sem_A2_right) fliplr((avg_A2_right-sem_A2_right))], ...
        rgb('Violet'));
set(fh, 'EdgeColor', 'None');

pl2 = plot( t_event, avg_A2_right, 'color', rgb('Magenta'));

ylabel('Firing rate (sp/s)');
title( 'Right turning' );
% axis tight;
xlim( [t_event(1) t_event(end)] );

subplot(4,1,4);
hold on;

avg_A1_left = squeeze(mean(left_A1));
sem_A1_left = get_sem( left_A1, 1 );

avg_A2_left = squeeze(mean(left_A2));
sem_A2_left = get_sem( left_A2, 1 );

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_A1_left+sem_A1_left) fliplr((avg_A1_left-sem_A1_left))], ...
        rgb('Gray'));
set(fh, 'EdgeColor', 'None');

pl1 = plot( t_event, avg_A1_left, 'color', rgb('Black'));

fh = fill( [t_event, fliplr(t_event)], ... 
        [(avg_A2_left+sem_A2_left) fliplr((avg_A2_left-sem_A2_left))], ...
        rgb('Violet'));
set(fh, 'EdgeColor', 'None');

pl2 = plot( t_event, avg_A2_left, 'color', rgb('Magenta'));
ylabel('Firing rate (sp/s)');
title( 'Left turning' );
% axis tight;
xlabel('Time (s)');
xlim( [t_event(1) t_event(end)] );

legend( [pl1, pl2], ['A1'], ['A2'], 'Location', 'northwest' );

dirpath = directories_to_analyze{1}{1};

saveas(f, [ analysis_path '/' dirpath '_A1_A2_stationary_to_moving_analysis_PSTH.fig']);
saveas(f, [ analysis_path '/' dirpath '_A1_A2_stationary_to_moving_analysis_PSTH.png']);
end



%% Plot A1 vs A2 psth

f = figure;

plot(A1_psth_down, A2_psth_down, 'o', 'MarkerSize', 3);
xlabel('A1 PSTH (spikes/s)');
ylabel('A2 PSTH (spikes/s)');
title('A1 vs A2 PSTH');

saveas(f, [analysis_path '/A1_vs_A2_psth_scatter.fig']);
saveas(f, [analysis_path '/A1_vs_A2_psth_scatter.png']);

%%
FILT_FACTOR = 0.04;

VmFilt_A2 = medfilt1( ephys_all_A{1}, FILT_FACTOR * ephys_SR, 'truncate' );
VmFilt_A1 = medfilt1( ephys_all_B{1}, FILT_FACTOR * ephys_SR, 'truncate' );

VmFilt_A2_corr = VmFilt_A2 - mean(VmFilt_A2);
VmFilt_A1_corr = VmFilt_A1 - mean(VmFilt_A1);

%VmFilt_A2_corr = VmFilt_A2;
%VmFilt_A1_corr = VmFilt_A1;

START = 1;
END = 60 * ephys_SR;

% figure;
% hold on;
% plot(t_all{1}(START:END), ephys_all_A{1}(START:END), 'b' );
% plot(t_all{1}(START:END), VmFilt_A2(START:END), 'g' );

BIN_SIZE = 0.050; % s
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

t_down = squeeze(mean(reshape(t_all{1}, [DT_EPHYS, length(t_all{1})/DT_EPHYS]), 1));
A2_Vm_down = squeeze(mean(reshape(VmFilt_A2_corr, [ DT_EPHYS, length(VmFilt_A2)/DT_EPHYS ] ),1));
A1_Vm_down = squeeze(mean(reshape(VmFilt_A1_corr, [ DT_EPHYS, length(VmFilt_A1)/DT_EPHYS ] ),1));

if 0
f = figure;

plot( A1_Vm_down, A2_Vm_down, 'o', 'MarkerSize', 3 );

xlabel('A1 Vm (delta mV)');
ylabel('A2 Vm (delta mV)');
title(['Bin size: ' num2str(BIN_SIZE) ' ms']);
saveas(f, [summary_analysis_path '/A1_A2_scatter_plot.fig']);
saveas(f, [summary_analysis_path '/A1_A2_scatter_plot.png']);
end

yaw_all_down = squeeze(mean(reshape(yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));
fwd_all_down = squeeze(mean(reshape(fwd_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));

figure;
histogram(yaw_all_down, 200);

%b = glmfit( A1_Vm_down, A2_Vm_down, 'normal' );
%plot(b);

%% Plots with a yaw overlay
yaw_range = [-300 300];
RBC_SIZE = 300;
rbc = jet(RBC_SIZE);

f = figure; 
hold on;

YAW_CUTOFF = 150;

SHIFT_FACTOR = 3;
rbc_index_array = [];
for ii = 1:(length(A2_Vm_down)-SHIFT_FACTOR)
    cur_yaw = yaw_all_down( ii+SHIFT_FACTOR );
    cur_rbc_index = ceil((cur_yaw + 301)/2.0);
    %rbc_index_array(ii) = cur_rbc_index;
    
     if(( cur_yaw > -1.0*YAW_CUTOFF ) && ( cur_yaw < YAW_CUTOFF ))
     %if(( cur_yaw < -1.0*YAW_CUTOFF ) && ( cur_yaw > YAW_CUTOFF ))
         continue;
     end
    
    rbc_index_array(end+1) = cur_rbc_index;
    
    if( cur_rbc_index < 1 )
        cur_rbc_index = 1;
    elseif( cur_rbc_index > RBC_SIZE )
        cur_rbc_index = RBC_SIZE;
    end
    
    cur_clr = rbc( cur_rbc_index, : );
    plot(A1_Vm_down(ii), A2_Vm_down(ii), 'o', 'MarkerSize', 3, 'color', cur_clr);
end

% ff = figure;
% histogram(rbc_index_array, 200);
% waitforbuttonpress;
% close(ff);

xlabel('A1 Vm (delta mV)');
ylabel('A2 Vm (delta mV)');
title(['Yaw cutoff: ' num2str(YAW_CUTOFF) '  bin size: ' num2str( BIN_SIZE ) ' s']);
saveas(f, [summary_analysis_path '/A1_A2_scatter_plot_with_yaw_clr_ ' num2str(idx) '.fig']);
saveas(f, [summary_analysis_path '/A1_A2_scatter_plot_with_yaw_clr_ ' num2str(idx) '.png']);

idx = idx+1;

%% Predicting yaw from A1/A2 ephys

%% 1. How well does A1 predict yaw?

cur_yaw = yaw_all_down( SHIFT_FACTOR:end );
cur_A1  = A1_Vm_down(1:end-SHIFT_FACTOR+1);

YAW_CUTOFF = 75;
yaw_idx = find( ( cur_yaw < -1.0*YAW_CUTOFF ) | ( cur_yaw > YAW_CUTOFF ) );

figure;
plot( cur_A1(yaw_idx), cur_yaw(yaw_idx), 'x' );

xlabel('A1 Vm (delta mV)');
ylabel('Yaw (deg/s)');

%% 2. How well does A2 predict yaw?
cur_yaw = yaw_all_down( SHIFT_FACTOR:end );
cur_A2  = A2_Vm_down(1:end-SHIFT_FACTOR+1);

figure;
plot( cur_A2, cur_yaw, 'x' );
xlabel('A2 Vm (delta mV)');
ylabel('Yaw (deg/s)');


%% 3. How well does A1/A2 predict yaw?
cur_yaw = yaw_all_down( SHIFT_FACTOR:end );
cur_A1  = A1_Vm_down(1:end-SHIFT_FACTOR+1);
cur_A2  = A2_Vm_down(1:end-SHIFT_FACTOR+1);

%% 
END = length(cur_yaw)/64;

Y = cur_yaw(1:END)';
g1 = cur_A1(1:END);
g2 = cur_A2(1:END);

[p, tbl, stats ] = anovan(Y, {g1,g2}, 'model', 'interaction');



%% GLM fit of yaw
% A1
Y = cur_yaw';
x = [ cur_A1' ];

distr = 'normal';
[b,dev_A1,stats] = glmfit( x, Y, distr );

% A2
Y = cur_yaw';
x = [ cur_A2' ];

[b,dev_A2,stats] = glmfit( x, Y, distr );

% A1+A2
Y = cur_yaw';
x = [ cur_A1', cur_A2' ];

[b,dev_A1_A2,stats] = glmfit( x, Y, distr );

format long;
disp(['dev(A1): ' num2str(dev_A1)]);
disp(['dev(A2): ' num2str(dev_A2)]);
disp(['dev(A1+A2): ' num2str(dev_A1_A2)]);

%% Plots with a yaw overlay as lines

% divide the data by bouts of left/right turns (yaw data), plot these as
% lines

[left_yaw_ids, right_yaw_ids ] = parse_left_right_turns( t_down, yaw_all_down );

yaw_range = [-300 300];
RBC_SIZE = 300;
rbc = jet(RBC_SIZE);

if 0
YAW_CUTOFF = 120;
SHIFT_FACTOR = 3;
rbc_index_array = [];
for ii = 1:(length(A2_Vm_down)-SHIFT_FACTOR)
    cur_yaw = yaw_all_down( ii+SHIFT_FACTOR );
    cur_rbc_index = ceil((cur_yaw + 301)/2.0);
    %rbc_index_array(ii) = cur_rbc_index;
    
    if(( cur_yaw > -1.0*YAW_CUTOFF ) && ( cur_yaw < YAW_CUTOFF ))
        continue;
    end
    
    rbc_index_array(end+1) = cur_rbc_index;
    
    if( cur_rbc_index < 1 )
        cur_rbc_index = 1;
    elseif( cur_rbc_index > RBC_SIZE )
        cur_rbc_index = RBC_SIZE;
    end
    
    cur_clr = rbc( cur_rbc_index, : );
    plot(A1_Vm_down(ii), A2_Vm_down(ii), 'o', 'MarkerSize', 3, 'color', cur_clr);
end
end

% ff = figure;
% histogram(rbc_index_array, 200);
% waitforbuttonpress;
% close(ff);

xlabel('A1 Vm (delta mV)');
ylabel('A2 Vm (delta mV)');
title(['Bin size: ' num2str( BIN_SIZE ) ' ms']);
saveas(f, [summary_analysis_path '/A1_A2_scatter_plot_with_yaw_as_lines.fig']);
saveas(f, [summary_analysis_path '/A1_A2_scatter_plot_with_yaw_as_lines.png']);

%% Plots with a fwd vel overlay

fwd_range = [0 30];
RBC_SIZE = 30;
rbc = jet(RBC_SIZE);

f = figure; 
hold on;

SHIFT_FACTOR = 3;
rbc_index_array = [];
for ii = 1:(length(A2_Vm_down)-SHIFT_FACTOR)
    cur_fwd = fwd_all_down( ii+SHIFT_FACTOR );
    cur_rbc_index = ceil(cur_fwd);
        
    if( cur_rbc_index < 1 )
        cur_rbc_index = 1;
    elseif( cur_rbc_index > RBC_SIZE )
        cur_rbc_index = RBC_SIZE;
    end
    
    cur_clr = rbc( cur_rbc_index, : );
    plot(A1_Vm_down(ii), A2_Vm_down(ii), 'o', 'MarkerSize', 3, 'color', cur_clr);
end

% ff = figure;
% histogram(rbc_index_array, 200);
% waitforbuttonpress;
% close(ff);

xlabel('A1 Vm (delta mV)');
ylabel('A2 Vm (delta mV)');
title(['Forward overlay, bin size: ' num2str( BIN_SIZE ) ' ms']);
colorbar;
saveas(f, [summary_analysis_path '/A1_A2_scatter_plot_with_fwd_clr.fig']);
saveas(f, [summary_analysis_path '/A1_A2_scatter_plot_with_fwd_clr.png']);

%% Generate the scatter plot
figure;

plot(A1_psth_down, A2_psth_down, 'o', 'MarkerSize', 3 );

xlabel('A1 Firing Rate (spikes/s)');
ylabel('A2 Firing Rate (spikes/s)');

%% Remove data points where yaw is 'small' to make it easier to fit.

BIN_SIZE = 0.05;
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

SHIFT_FACTOR = 3;

A1_psth_down = squeeze(mean(reshape( A1_psth, [ DT_YAW, length(A1_psth)/DT_YAW ] ), 1));
A2_psth_down = squeeze(mean(reshape( A2_psth, [ DT_YAW, length(A2_psth)/DT_YAW ] ), 1));
fwd_all_down = squeeze(mean(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]), 1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]), 1));

A1_choosen = zeros( 1, length(A2_psth_down) );
A2_choosen = zeros( 1, length(A2_psth_down) );
yaw_choosen = zeros( 1, length(A2_psth_down) );

YAW_CUTOFF = 250;
cnt = 0;
for i = 1:length(A2_psth_down)-SHIFT_FACTOR-1
    
    cur_yaw = yaw_all_down( i+SHIFT_FACTOR );
    
    if( (cur_yaw < -1 * YAW_CUTOFF) || (cur_yaw > YAW_CUTOFF) )
        A1_choosen( cnt+1 ) = A1_psth_down(i);
        A2_choosen( cnt+1 ) = A2_psth_down(i);
        yaw_choosen( cnt+1 ) = cur_yaw;
        cnt = cnt + 1;
    end    
end


%%%%%%%%
% Using fit() for yaw
%%%%%%%%
Y = yaw_choosen(1:cnt)';

% A1 + A2
X = horzcat(A1_choosen(1:cnt)', A2_choosen(1:cnt)' );
[fobj, gof_A12_yaw] = fit(X,Y, 'poly11');

% A1 
X = A1_choosen(1:cnt)';
[fobj, gof_A1_yaw] = fit(X,Y, 'poly1');

% A2
X = A2_choosen(1:cnt)';
[fobj, gof_A2_yaw] = fit(X,Y, 'poly1');

A1_choosen = zeros( 1, length(A2_psth_down) );
A2_choosen = zeros( 1, length(A2_psth_down) );
fwd_choosen = zeros( 1, length(A2_psth_down) );

FWD_CUTOFF = -2;
cnt = 0;
for i = 1:length(A2_psth_down)-SHIFT_FACTOR
    
    cur_fwd = fwd_all_down( i+SHIFT_FACTOR );
    
    if(cur_fwd > FWD_CUTOFF)
        A1_choosen( cnt+1 ) = A1_psth_down(i);
        A2_choosen( cnt+1 ) = A2_psth_down(i);
        fwd_choosen( cnt+1 ) = cur_fwd;
        cnt = cnt + 1;
    end    
end


%%%%%%%%
% Using fit() for fwd
%%%%%%%%
Y = fwd_choosen(1:cnt)';

% A1 + A2
X = horzcat(A1_choosen(1:cnt)', A2_choosen(1:cnt)' );
[fobj, gof_A12_fwd] = fit(X,Y, 'poly11');

% A1 
X = A1_choosen(1:cnt)';
[fobj, gof_A1_fwd] = fit(X,Y, 'poly1');

figure;
plot( X, Y, 'o', 'MarkerSize', 3 )
xlabel('A1 (spikes/s)');
ylabel('yaw (deg/s)');


% A2
X = A2_choosen(1:cnt)';
[fobj, gof_A2_fwd] = fit(X,Y, 'poly1');

f = figure;
gofs_yaw = [ gof_A12_yaw.rsquare gof_A1_yaw.rsquare gof_A2_yaw.rsquare ];
gofs_fwd = [ gof_A12_fwd.rsquare gof_A1_fwd.rsquare gof_A2_fwd.rsquare ];
bb = bar([gofs_yaw; gofs_fwd]);
ylabel('R^2 val');
hh = title(['A1 and A2 psth linear regression of (yaw cutoff: ' num2str(YAW_CUTOFF) ', fwd curoff: ' num2str(FWD_CUTOFF) ' )']);
set(hh, 'Interpreter', 'none');
set(gca(), 'XTickLabel', {'Yaw', 'Fwd'})

xx = get(gca(), 'Children');
hLegend = legend(xx, {'A2', 'A1', 'A1+A2'});
    
saveas( f, [analysis_path '/A2_A1_psth_right_vs_left_vs_both_rsquared_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.fig'] );
saveas( f, [analysis_path '/A2_A1_psth_right_vs_left_vs_both_rsquared_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.png'] );

f = figure;
yaw_rmse = [ gof_A12_yaw.rmse gof_A1_yaw.rmse gof_A2_yaw.rmse ];
fwd_rmse = [ gof_A12_fwd.rmse gof_A1_fwd.rmse gof_A2_fwd.rmse ];
bb = bar([yaw_rmse; fwd_rmse]);
ylabel('RMSE val');
hh = title(['A1 and A2 PSTH Root Mean Square Error (rmse) (yaw cutoff: ' num2str(YAW_CUTOFF) ', fwd curoff: ' num2str(FWD_CUTOFF) ' )']);
set(hh, 'Interpreter', 'none');
set(gca(), 'XTickLabel', {'Yaw', 'Fwd'})                

xx = get(gca(), 'Children');
hLegend = legend(xx, {'A2', 'A1', 'A1+A2'});
    
saveas( f, [analysis_path '/A2_A1_PSTH_rmse_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.fig'] );
saveas( f, [analysis_path '/A2_A1_PSTH_rmse_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.png'] );



%% Vm vs yaw, fwd A1 and A2

BIN_SIZE = 0.05;
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

SHIFT_FACTOR = 3;

A2_Vm = medfilt1(ephys_all_A{1}, 0.4 * ephys_SR, 'truncate');
A1_Vm = medfilt1(ephys_all_B{1}, 0.4 * ephys_SR, 'truncate');

A1_A2_diff = A1_Vm - A2_Vm;

A1_Vm_down = squeeze(mean(reshape( A1_Vm, [ DT_EPHYS, length(A1_A2_diff)/DT_EPHYS ] ),1));
A2_Vm_down = squeeze(mean(reshape( A2_Vm, [ DT_EPHYS, length(A1_A2_diff)/DT_EPHYS ] ),1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));
fwd_all_down = squeeze(mean(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]), 1));

A1_choosen = zeros( 1, length(A2_psth_down) );
A2_choosen = zeros( 1, length(A2_psth_down) );
yaw_choosen = zeros( 1, length(A2_psth_down) );

YAW_CUTOFF = 250;
cnt = 0;
for i = 1:length(A2_psth_down)-SHIFT_FACTOR-1
    
    cur_yaw = yaw_all_down( i+SHIFT_FACTOR );
    
    if( (cur_yaw < -1 * YAW_CUTOFF) || (cur_yaw > YAW_CUTOFF) )
        A1_choosen( cnt+1 ) = A1_Vm_down(i);
        A2_choosen( cnt+1 ) = A2_Vm_down(i);
        yaw_choosen( cnt+1 ) = cur_yaw;
        cnt = cnt + 1;
    end    
end


%%%%%%%%
% Using fit() for yaw
%%%%%%%%
Y = yaw_choosen(1:cnt)';

% A1 + A2
X = horzcat(A1_choosen(1:cnt)', A2_choosen(1:cnt)' );
[fobj, gof_A12_yaw] = fit(X,Y, 'poly11');
cc = coeffvalues(fobj);
slope_A12_yaw = cc(2);

% A1 
X = A1_choosen(1:cnt)';
[fobj, gof_A1_yaw] = fit(X,Y, 'poly1');
cc = coeffvalues(fobj);
slope_A1_yaw = cc(1);

% A2
X = A2_choosen(1:cnt)';
[fobj, gof_A2_yaw] = fit(X,Y, 'poly1');
cc = coeffvalues(fobj);
slope_A2_yaw = cc(1);

A1_choosen = zeros( 1, length(A2_psth_down) );
A2_choosen = zeros( 1, length(A2_psth_down) );
fwd_choosen = zeros( 1, length(A2_psth_down) );

FWD_CUTOFF = -2;
cnt = 0;
for i = 1:length(A2_psth_down)-SHIFT_FACTOR
    
    cur_fwd = fwd_all_down( i+SHIFT_FACTOR );
    
    if(cur_fwd > FWD_CUTOFF)
        A1_choosen( cnt+1 ) = A1_Vm_down(i);
        A2_choosen( cnt+1 ) = A2_Vm_down(i);
        fwd_choosen( cnt+1 ) = cur_fwd;
        cnt = cnt + 1;
    end    
end


%%%%%%%%
% Using fit() for fwd
%%%%%%%%
Y = fwd_choosen(1:cnt)';

% A1 + A2
X = horzcat(A1_choosen(1:cnt)', A2_choosen(1:cnt)' );
[fobj, gof_A12_fwd] = fit(X,Y, 'poly11');
cc = coeffvalues(fobj);
slope_A12_fwd = cc(2);

% A1 
X = A1_choosen(1:cnt)';
[fobj, gof_A1_fwd] = fit(X,Y, 'poly1');
cc = coeffvalues(fobj);
slope_A1_fwd = cc(1);

% A2
X = A2_choosen(1:cnt)';
[fobj, gof_A2_fwd] = fit(X,Y, 'poly1');
cc = coeffvalues(fobj);
slope_A2_fwd = cc(1);

f = figure;
gofs_yaw = [ gof_A12_yaw.rsquare gof_A1_yaw.rsquare gof_A2_yaw.rsquare ];
gofs_fwd = [ gof_A12_fwd.rsquare gof_A1_fwd.rsquare gof_A2_fwd.rsquare ];
bb = bar([gofs_yaw; gofs_fwd]);
ylabel('R^2 val');
hh = title(['A1 and A2 Vm linear regression of (yaw cutoff: ' num2str(YAW_CUTOFF) ', fwd curoff: ' num2str(FWD_CUTOFF) ' )']);
set(hh, 'Interpreter', 'none');
set(gca(), 'XTickLabel', {'Yaw', 'Fwd'})                

xx = get(gca(), 'Children');
hLegend = legend(xx, {'A2', 'A1', 'A1+A2'});
    
saveas( f, [analysis_path '/A2_A1_Vm_rsquared_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.fig'] );
saveas( f, [analysis_path '/A2_A1_Vm_rsquared_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.png'] );

f = figure;
yaw_rmse = [ gof_A12_yaw.rmse gof_A1_yaw.rmse gof_A2_yaw.rmse ];
fwd_rmse = [ gof_A12_fwd.rmse gof_A1_fwd.rmse gof_A2_fwd.rmse ];
bb = bar([yaw_rmse; fwd_rmse]);
ylabel('RMSE val');
hh = title(['A1 and A2 Vm Root Mean Square Error (rmse) (yaw cutoff: ' num2str(YAW_CUTOFF) ', fwd curoff: ' num2str(FWD_CUTOFF) ' )']);
set(hh, 'Interpreter', 'none');
set(gca(), 'XTickLabel', {'Yaw', 'Fwd'})                

xx = get(gca(), 'Children');
hLegend = legend(xx, {'A2', 'A1', 'A1+A2'});
    
saveas( f, [analysis_path '/A2_A1_Vm_rmse_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.fig'] );
saveas( f, [analysis_path '/A2_A1_Vm_rmse_fit()_regress_yaw_cutoff_' num2str(YAW_CUTOFF) '_fwd_cutoff_' num2str(FWD_CUTOFF) '_bars.png'] );




%% Plot A1 and A2 PSTH with a yaw vel overlay

% Assuming yaw range [-1000 1000]

RBC_SIZE = 2000;

BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;

A1_psth_down = squeeze(mean(reshape( A1_psth, [ DT_YAW, length(A2_psth)/DT_YAW ] ),1));
A2_psth_down = squeeze(mean(reshape( A2_psth, [ DT_YAW, length(A2_psth)/DT_YAW ] ),1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));

rbc = jet(RBC_SIZE);
SHIFT_FACTOR = 3;
hold on;

A2_yaw_colors = zeros( length(A2_psth_down), 3 );

A1_choosen = zeros( 1, length(A2_psth_down) );
A2_choosen = zeros( 1, length(A2_psth_down) );
A2_yaw_colors_choosen = zeros( length(A2_psth_down), 3 );

YAW_CUTOFF = 250;
cnt = 0;
for i = 1:length(A2_psth_down)-SHIFT_FACTOR
    
    cur_yaw = yaw_all_down( i+SHIFT_FACTOR );
    cur_rbc_index = ceil(cur_yaw + 1000);
    
    if( cur_rbc_index < 1 )
        cur_rbc_index = 1;
    elseif( cur_rbc_index > RBC_SIZE )
        cur_rbc_index = RBC_SIZE;
    end
    
    cur_clr = rbc( cur_rbc_index, : );
    
    if( (cur_yaw < -1 * YAW_CUTOFF) || (cur_yaw > YAW_CUTOFF) )
        A1_choosen( cnt+1 ) = A1_psth_down(i);
        A2_choosen( cnt+1 ) = A2_psth_down(i);
        A2_yaw_colors_choosen( cnt+1, : ) = cur_clr;
        cnt = cnt + 1;
    end    
end

f = figure;
rbc = jet(RBC_SIZE);
% scatter(A2_L_psth_down, A2_R_psth_down, 4, A2_yaw_colors);
scatter( A1_choosen, A2_choosen, 4, A2_yaw_colors_choosen);
colormap('jet')
h = colorbar;
ylabel(h, 'Yaw vel (deg/s)');
caxis([-1000 1000]);
grid on;
xlabel('A1 PSTH (spikes/s)');
ylabel('A2 PSTH (spikes/s)');
title(['Yaw velocity overlay: yaw cutoff: ' num2str(YAW_CUTOFF) ]);

saveas(f,[analysis_path '/A1_A2_PSTH_scatter_with_yaw_overlay_cutoff_' num2str(YAW_CUTOFF) '.fig']);
saveas(f,[analysis_path '/A1_A2_PSTH_scatter_with_yaw_overlay_cutoff_' num2str(YAW_CUTOFF) '.png']);


%% Show relationship between left and right A2 neurons Vm difference and yaw

BIN_SIZE = 0.05;
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

A2_Vm = medfilt1(ephys_all_A{1}, 0.4 * ephys_SR, 'truncate');
A1_Vm = medfilt1(ephys_all_B{1}, 0.4 * ephys_SR, 'truncate');

A1_A2_diff = A1_Vm - A2_Vm;

A1_A2_diff_down = squeeze(mean(reshape( A1_A2_diff, [ DT_EPHYS, length(A1_A2_diff)/DT_EPHYS ] ),1));
A1_Vm_down = squeeze(mean(reshape( A1_Vm, [ DT_EPHYS, length(A1_A2_diff)/DT_EPHYS ] ),1));
A2_Vm_down = squeeze(mean(reshape( A2_Vm, [ DT_EPHYS, length(A1_A2_diff)/DT_EPHYS ] ),1));
yaw_all_down = squeeze(mean(reshape(yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));

f = figure;

subplot(1,3,1);
plot( A1_Vm_down(1:end-3), yaw_all_down(4:end), 'o', 'MarkerSize', 3 );
xlabel('A1 Vm (mV)');
ylabel('Yaw (deg/s)');
title('A1 vs. yaw');

subplot(1,3,2);
plot( A2_Vm_down(1:end-3), yaw_all_down(4:end), 'o', 'MarkerSize', 3 );
xlabel('A2 Vm (mV)');
ylabel('Yaw (deg/s)');
title('A2 vs. yaw');

subplot(1,3,3);
plot( A1_A2_diff_down(1:end-3), yaw_all_down(4:end), 'o', 'MarkerSize', 3 );
xlabel('A1-A2 Vm (mV)');
ylabel('Yaw (deg/s)');
title('A1-A2 vs. yaw');

saveas(f,[analysis_path '/A1_A2_Vm_diff_vs_yaw_scatter.fig']);
saveas(f,[analysis_path '/A1_A2_Vm_diff_vs_yaw_scatter.png']);


%% Show relationship between left and right A2 neurons PSTH difference and yaw

BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;

A1_psth_down = squeeze(mean(reshape( A1_psth, [ DT_YAW, length(A1_psth)/DT_YAW ] ), 1));
A2_psth_down = squeeze(mean(reshape( A2_psth, [ DT_YAW, length(A2_psth)/DT_YAW ] ), 1));
fwd_all_down = squeeze(mean(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]), 1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]), 1));


A1_A2_diff_psth = A1_psth_down - A2_psth_down;

f = figure;

subplot(1,3,1);
plot(A1_psth_down(1:end-3), yaw_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('A1 PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
title('A1 PSTH vs. yaw');
grid on;

subplot(1,3,2);
plot(A2_psth_down(1:end-3), yaw_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('A2 PSTH (spikes/s)');
ylabel('Yaw (deg/s)');
title('A2 PSTH vs. yaw');
grid on;

subplot(1,3,3);
plot( A1_A2_diff_psth(1:end-3), yaw_all_down(4:end), 'o', 'MarkerSize', 3);
xlabel('A1-A2 PSTH diff (spikes/s)');
ylabel('Yaw (deg/s)');
title('A1-A2 vs. yaw');
grid on;

saveas( f, [analysis_path '/A1_A2_PSTH_diff_vs_yaw_scatter.fig'] );
saveas( f, [analysis_path '/A1_A2_PSTH_diff_vs_yaw_scatter.png'] );


%% Show relationship between left and right A2 neurons PSTH difference and fwd

BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;
FWD_MIN = -50;
FWD_MAX = 100;

A1_psth_down = squeeze(mean(reshape( A1_psth, [ DT_YAW, length(A1_psth)/DT_YAW ] ), 1));
A2_psth_down = squeeze(mean(reshape( A2_psth, [ DT_YAW, length(A2_psth)/DT_YAW ] ), 1));
fwd_all_down = squeeze(mean(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]), 1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]), 1));

A1_A2_diff_psth = A1_psth_down - A2_psth_down;

f = figure;

subplot(1,3,1);
plot(A1_psth_down(1:end-3), fwd_all_down(4:end), 'o', 'MarkerSize', 3);
ylim([FWD_MIN FWD_MAX]);
xlabel('A1 PSTH (spikes/s)');
ylabel('Fwd (mm/s)');
title('A1 PSTH vs. fwd');
grid on;            


subplot(1,3,2);
plot(A2_psth_down(1:end-3), fwd_all_down(4:end), 'o', 'MarkerSize', 3);
ylim([FWD_MIN FWD_MAX]);
xlabel('A2 PSTH (spikes/s)');
ylabel('Fwd (mm/s)');
title('A2 PSTH vs. fwd');
grid on;

subplot(1,3,3);
plot( A1_A2_diff_psth(1:end-3), fwd_all_down(4:end), 'o', 'MarkerSize', 3);
ylim([FWD_MIN FWD_MAX]);
xlabel('A1-A2 PSTH diff (spikes/s)');
ylabel('Fwd (mm/s)');
title('A1-A2 vs. fwd');
grid on;

saveas( f, [analysis_path '/A1_A2_PSTH_diff_vs_fwd_scatter.fig'] );
saveas( f, [analysis_path '/A1_A2_PSTH_diff_vs_fwd_scatter.png'] );

%% 

BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;

%A2_L_psth_base = A2_L_psth - mean(A2_L_psth);
%A2_R_psth_base = A2_R_psth - mean(A2_R_psth);

A2_psth_base = A2_psth;
A1_psth_base = A1_psth;

A2_psth_down = squeeze(mean(reshape( A2_psth_base, [ DT_YAW, length(A2_psth_base) / DT_YAW ] ), 1));
A1_psth_down = squeeze(mean(reshape( A1_psth_base, [ DT_YAW, length(A2_psth_base) / DT_YAW ] ), 1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( yaw_all{1} ) / DT_YAW ]), 1));
fwd_all_down = squeeze(mean(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} ) / DT_YAW ]), 1));
% fwd_all_down = squeeze(mean(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} ) / DT_YAW ]), 1));
fwd_all_down_std = squeeze(std(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} ) / DT_YAW ]), 0, 1));

% figure;
% hist(fwd_all_down_std, 1000);

cnt = 0;
A2_psth_choosen = zeros(1, length( A2_psth_down ));
A1_psth_choosen = zeros(1, length( A2_psth_down ));
yaw_choosen  = zeros(1, length( yaw_all_down ));
fwd_choosen  = zeros(1, length( yaw_all_down ));

PSTH_BASELINE_CUTOFF = 1.0;
SHIFT_FACTOR = 3;
mcnt = 0;

WINDOW_SIZE = 5;
FWD_STD_THRESHOLD = 0.01;

for i = 1:length( A2_psth_down )-(SHIFT_FACTOR+1)-WINDOW_SIZE
    
    cur_A2_psth = A2_psth_down( i );
    cur_A1_psth = A1_psth_down( i );
    
    cur_fwd_std = fwd_all_down_std( i );
    
    if( abs(cur_fwd_std) < FWD_STD_THRESHOLD )
       continue; 
    end
    
%    if( (( cur_A2_L_psth <= -1.0*PSTH_BASELINE_CUTOFF ) || ( cur_A2_L_psth >= PSTH_BASELINE_CUTOFF )) && ...
%        (( cur_A2_R_psth <= -1.0*PSTH_BASELINE_CUTOFF ) || ( cur_A2_R_psth >= PSTH_BASELINE_CUTOFF )) )
        A2_psth_choosen( mcnt+1 ) = cur_A2_psth;
        A1_psth_choosen( mcnt+1 ) = cur_A1_psth;
        yaw_choosen( mcnt+1 ) = yaw_all_down(i+SHIFT_FACTOR);
        fwd_choosen( mcnt+1 ) = fwd_all_down(i+SHIFT_FACTOR);
        mcnt = mcnt + 1;
%    end        
end


GREATER_PSTH_CUTOFF = 25;
LESSER_PSTH_CUTOFF = 5;
A1_PSTH_CUTOFF = 25;
A2_PSTH_CUTOFF = A1_PSTH_CUTOFF;
OTHER_PSTH_CUTOFF = 5;

yaw_in_greater_cutoff = zeros(1, length(yaw_choosen));
fwd_in_greater_cutoff = zeros(1, length(fwd_choosen));

yaw_in_lesser_cutoff = zeros(1, length(yaw_choosen));
fwd_in_lesser_cutoff = zeros(1, length(fwd_choosen));

yaw_in_A1_cutoff = zeros(1, length(yaw_choosen));
fwd_in_A1_cutoff = zeros(1, length(fwd_choosen));

yaw_in_A2_cutoff = zeros(1, length(yaw_choosen));
fwd_in_A2_cutoff = zeros(1, length(fwd_choosen));

cur_g_idx = 1;
cur_l_idx = 1;
cur_A1_idx = 1;
cur_A2_idx = 1;
for i = 1:length(A2_L_psth_choosen)
    
    cur_A2_psth = A2_psth_choosen(i);
    cur_A1_psth = A1_psth_choosen(i);
    
    if( (cur_A2_psth > GREATER_PSTH_CUTOFF) && (cur_A1_psth > GREATER_PSTH_CUTOFF) )
        yaw_in_greater_cutoff(cur_g_idx) = yaw_choosen(i);
        fwd_in_greater_cutoff(cur_g_idx) = fwd_choosen(i);  
        cur_g_idx = cur_g_idx + 1;
    end
    
    if( (cur_A2_psth < LESSER_PSTH_CUTOFF) && (cur_A1_psth < LESSER_PSTH_CUTOFF) )
        yaw_in_lesser_cutoff(cur_l_idx) = yaw_choosen(i);
        fwd_in_lesser_cutoff(cur_l_idx) = fwd_choosen(i);  
        cur_l_idx = cur_l_idx + 1;
    end
    
    if( (cur_A2_psth > A2_PSTH_CUTOFF) )
        yaw_in_A2_cutoff(cur_A2_idx) = yaw_choosen(i);
        fwd_in_A2_cutoff(cur_A2_idx) = fwd_choosen(i);  
        cur_A2_idx = cur_A2_idx + 1;
    end
    
    if( (cur_A1_psth > A1_PSTH_CUTOFF) )
        yaw_in_A1_cutoff(cur_A1_idx) = yaw_choosen(i);
        fwd_in_A1_cutoff(cur_A1_idx) = fwd_choosen(i);  
        cur_A1_idx = cur_A1_idx + 1;
    end
end

NBINS = 60;
f = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,4,1);
hist( yaw_in_greater_cutoff(1:cur_g_idx), NBINS );
xlabel('Yaw vel (deg/s)');
ylabel('Count');
title(['A2 and A1 PSTH > ' num2str(GREATER_PSTH_CUTOFF)]);
legend(['Mean: ' num2str(mean(yaw_in_greater_cutoff(1:cur_g_idx))) ' std: ' num2str(std(yaw_in_greater_cutoff(1:cur_g_idx))) ]);


subplot(2,4,2);
hist( yaw_in_lesser_cutoff(1:cur_l_idx), NBINS );
xlabel('Yaw vel (deg/s)');
ylabel('Count');
title(['A2 and A1 PSTH < ' num2str(LESSER_PSTH_CUTOFF)]);
legend(['Mean: ' num2str(mean(yaw_in_left_cutoff(1:cur_l_idx))) ' std: ' num2str(std(yaw_in_left_cutoff(1:cur_l_idx))) ]);

subplot(2,4,3);
hist( yaw_in_A2_cutoff(1:cur_g_idx), NBINS );
xlabel('Yaw vel (deg/s)');
ylabel('Count');
title(['A2 PSTH > ' num2str(A2_PSTH_CUTOFF)]);
legend(['Mean: ' num2str(mean(yaw_in_A2_cutoff(1:cur_A2_idx))) ' std: ' num2str(std(yaw_in_A2_cutoff(1:cur_A2_idx))) ]);

subplot(2,4,4);
hist( yaw_in_A1_cutoff(1:cur_g_idx), NBINS );
xlabel('Yaw vel (deg/s)');
ylabel('Count');
title(['A1 PSTH > ' num2str(A2_PSTH_CUTOFF)]);
legend(['Mean: ' num2str(mean(yaw_in_A1_cutoff(1:cur_A1_idx))) ' std: ' num2str(std(yaw_in_A1_cutoff(1:cur_A1_idx))) ]);



subplot(2,4,5);
hist( fwd_in_greater_cutoff(1:cur_g_idx), NBINS );
xlabel('Fwd vel (mm/s)');
ylabel('Count');
legend(['Mean: ' num2str(mean(fwd_in_greater_cutoff(1:cur_g_idx))) ' std: ' num2str(std(fwd_in_greater_cutoff(1:cur_g_idx))) ]);


subplot(2,4,6);
hist( fwd_in_lesser_cutoff(1:cur_l_idx), NBINS );
xlabel('Fwd vel (mm/s)');
ylabel('Count');
legend(['Mean: ' num2str(mean(fwd_in_lesser_cutoff(1:cur_l_idx))) ' std: ' num2str(std(fwd_in_lesser_cutoff(1:cur_l_idx))) ]);

subplot(2,4,7);
hist( fwd_in_A2_cutoff(1:cur_g_idx), NBINS );
xlabel('Fwd vel (mm/s)');
ylabel('Count');
legend(['Mean: ' num2str(mean(fwd_in_A2_cutoff(1:cur_A2_idx))) ' std: ' num2str(std(fwd_in_A2_cutoff(1:cur_A2_idx))) ]);

subplot(2,4,8);
hist( fwd_in_A1_cutoff(1:cur_g_idx), NBINS );
xlabel('Fwd vel (mm/s)');
ylabel('Count');
legend(['Mean: ' num2str(mean(fwd_in_A1_cutoff(1:cur_A1_idx))) ' std: ' num2str(std(fwd_in_A1_cutoff(1:cur_A1_idx))) ]);

saveas(f, [analysis_path '/yaw_fwd_hist_in_A2_A1_mutual_PSTH_gcutoff_' num2str(GREATER_PSTH_CUTOFF) '_lcutoff_' num2str(LESSER_PSTH_CUTOFF) '.fig']);
saveas(f, [analysis_path '/yaw_fwd_hist_in_A2_A1_mutual_PSTH_gcutoff_' num2str(GREATER_PSTH_CUTOFF) '_lcutoff_' num2str(LESSER_PSTH_CUTOFF) '.png']);


%%

f = figure;
scatter(yaw_in_greater_cutoff, fwd_in_greater_cutoff, 3);
xlabel('Yaw (deg/s)');
ylabel('Fwd (mm/s)');
title(['A2 and A1 FR > ' num2str(GREATER_PSTH_CUTOFF)]);
saveas(f, [analysis_path '/yaw_fwd_scatter_when_A2_A1_FR_gcutoff_' num2str(GREATER_PSTH_CUTOFF) '.fig']);
saveas(f, [analysis_path '/yaw_fwd_scatter_when_A2_A1_FR_gcutoff_' num2str(GREATER_PSTH_CUTOFF) '.png']);

%% Plot yaw vs. fwd with A1 FR overlay

RBC_SIZE_A1 = 30;
RBC_SIZE_A2 = 50;

BIN_SIZE = 0.05;
DT_YAW   = ball_SR * BIN_SIZE;

A1_psth_down = squeeze(mean(reshape( A1_psth, [ DT_YAW, length(A2_psth)/DT_YAW ] ),1));
A2_psth_down = squeeze(mean(reshape( A2_psth, [ DT_YAW, length(A2_psth)/DT_YAW ] ),1));
yaw_all_down = squeeze(mean(reshape( yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));
fwd_all_down = squeeze(mean(reshape( fwd_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));

rbc_A1 = jet(RBC_SIZE_A1);
rbc_A2 = jet(RBC_SIZE_A2);
SHIFT_FACTOR = 3;
hold on;

fwd_choosen = zeros( 1, length(A2_psth_down) );
yaw_choosen = zeros( 1, length(A2_psth_down) );
A2_colors_choosen = zeros( length(A2_psth_down), 3 );
A1_colors_choosen = zeros( length(A2_psth_down), 3 );

YAW_CUTOFF = 250;
cnt = 0;
for i = 1:length(A2_psth_down)-SHIFT_FACTOR
    
    cur_A1_fr = A1_psth_down( i );
    cur_rbc_index_A1 = ceil(cur_A1_fr);
    
    if( cur_rbc_index_A1 < 1 )
        cur_rbc_index_A1 = 1;
    elseif( cur_rbc_index_A1 > RBC_SIZE_A1 )
        cur_rbc_index_A1 = RBC_SIZE_A1;
    end
    
    cur_clr_A1 = rbc_A1( cur_rbc_index_A1, : );

    cur_A2_fr = A2_psth_down( i );
    cur_rbc_index_A2 = ceil(cur_A2_fr);
    
    if( cur_rbc_index_A2 < 1 )
        cur_rbc_index_A2 = 1;
    elseif( cur_rbc_index_A2 > RBC_SIZE_A2 )
        cur_rbc_index_A2 = RBC_SIZE_A2;
    end
    
    cur_clr_A2 = rbc_A2( cur_rbc_index_A2, : );
    
    A1_colors_choosen( cnt+1, : ) = cur_clr_A1;
    A2_colors_choosen( cnt+1, : ) = cur_clr_A2;
    fwd_choosen(cnt+1) = fwd_all_down( i + SHIFT_FACTOR );
    yaw_choosen(cnt+1) = yaw_all_down( i + SHIFT_FACTOR );
    cnt = cnt + 1;
end

f = figure;
subplot(1,2,1)
rbc = jet(RBC_SIZE_A1);
scatter3( yaw_choosen, fwd_choosen, A1_psth_down);
colormap('jet')
h = colorbar;
view(2)
caxis([0 RBC_SIZE_A1]);
grid on;
xlabel('Yaw (deg/s)');
ylabel('Fwd (mm/s)');
title(['A1 FR color overlay ']);

subplot(1,2,2)
rbc = jet(RBC_SIZE_A2);
scatter3( yaw_choosen, fwd_choosen, A2_psth_down);
view(2);
colormap('jet')
h = colorbar;
caxis([0 RBC_SIZE_A2]);
grid on;
xlabel('Yaw (deg/s)');
ylabel('Fwd (mm/s)');
title(['A2 FR color overlay ']);

%saveas(f,[analysis_path '/A1_A2_PSTH_scatter_with_yaw_overlay_cutoff_' num2str(YAW_CUTOFF) '.fig']);
%saveas(f,[analysis_path '/A1_A2_PSTH_scatter_with_yaw_overlay_cutoff_' num2str(YAW_CUTOFF) '.png']);












































