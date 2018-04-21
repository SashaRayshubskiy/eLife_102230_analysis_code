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

directories_to_analyze = { { '170816_2xGFP_ss731_75C10_01', [0, 1], 2*11.5*150 } };
%directories_to_analyze = { { '170816_2xGFP_ss731_75C10_01', [0, 1], 11.5*10 } };

[t_all, t_vel_all, yaw_all, fwd_all, ephys_all_A, ephys_all_B] = load_LAL_DN_data( working_dir, directories_to_analyze, ephys_SR, ball_SR, 1 );

idx = 1;

%% Show histgram of fwd vel
figure;
hist(fwd_all{1}, 1000);

%% Convert physiology to PSTH
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

% VmFilt_A2_corr = VmFilt_A2 - mean(VmFilt_A2);
% VmFilt_A1_corr = VmFilt_A1 - mean(VmFilt_A1);

VmFilt_A2_corr = VmFilt_A2;
VmFilt_A1_corr = VmFilt_A1;

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

