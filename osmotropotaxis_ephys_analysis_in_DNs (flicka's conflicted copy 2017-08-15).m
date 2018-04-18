%% Load data
clear all;

global slash;
slash = '/';

%working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/descending_neurons/';
working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/SC_to_LAL_neurons/';
%working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/osmotropotaxis_behavior_only/';

cur_dir = '170812_lexAOpCsChrimson_2x_gfp_83blexA_R30C01_01';

datapath = [ working_dir cur_dir ];

analysis_path = [datapath '/analysis/'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

analysis_tt_path = [datapath '/analysis/trial_by_trial/'];

if(~exist(analysis_tt_path, 'dir'))
    mkdir(analysis_tt_path);
end

sids = [ 0 ];
trial_type_cnt = 3;

first_stim_t = 3.0;
last_stim_t =  3.5;

TRIAL_CNT_MAX = 170;
[ bdata_raw, bdata_time, trial_metadata ] = load_behavioral_data( sids, datapath, trial_type_cnt );

TIME_OFFSET = bdata_time(1);

% Transform raw data into velocity
t_volts = {};
volts_all = {};

t_vel_all = {};
yaw_vel_all = {};
forward_vel_all = {};     

for tt = 1:trial_type_cnt
    for trial = 1:size(bdata_raw{tt},1)
        
        tid = trial_metadata{tt}(trial,2);
         if(tid > TRIAL_CNT_MAX)
             continue;
         end
        
        t = bdata_time;
        D = squeeze(bdata_raw{tt}(trial,:,:));
        
        [currentA, voltageA, currentB, voltageB] = get_dual_scaled_voltage_and_current( D );

        
        volts_all{tt}(trial, :) = voltageA;        
        current_all{tt}(trial, :) = currentA;
        %volts_all{tt}(trial, :) = filtfilt(d,D(:,2));        
        t_volts{tt}(trial,:) = t-TIME_OFFSET;

        [ t_vel, vel_forward, vel_side, vel_yaw ] = get_velocity_ephysrig(t, D, datapath, 1 );
        
        t_vel_all{tt}( trial, : ) = t_vel-TIME_OFFSET;
        
        yaw_vel_all{tt}( trial, : ) = vel_yaw;
        forward_vel_all{tt}( trial, : ) = vel_forward;        
    end
end

%% Use this for single type runs

%display_avg_data(sids, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, trial_type_cnt );

display_avg_data_delta_v(sids, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, trial_type_cnt );
%display_avg_data_no_voltage(sids, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, trial_type_cnt );

%display_mean_vel_per_trial(sids, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, trial_type_cnt );

%%
path1 = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/SC_to_LAL_neurons/170717_R53H03_lexA_CsCr_lexAOp_01/analysis/delta_volt_mean.mat';
path2 = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/SC_to_LAL_neurons/170718_UAS-GFP_ss730_control_01/analysis/delta_volt_mean.mat';

p1 = load(path1);
p2 = load(path2);

f = figure;
hold on;
t = [0:1/4000:11.5];
t = t(1:end-1);

plot(t,p1.save_avg(1,:), 'DisplayName', 'R53H03-Chrimson/ss730-GFP', 'LineWidth', 3.0);
plot(t,p2.save_avg(1,:), 'DisplayName', 'ss730-GFP');
legend('show');

first_stim_t = 3.0;
last_stim_t = 3.5;

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

saveas(f, ['/Users/sasha/Dropbox/Wilson_lab/paper_1/data/SC_to_LAL_neurons/R53H03_vs_control.fig']);
saveas(f, ['/Users/sasha/Dropbox/Wilson_lab/paper_1/data/SC_to_LAL_neurons/R53H03_vs_control.png']);

%%
display_all_data_v2(sids, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, trial_type_cnt );

%% Display vel data and voltage data, show all data
display_all_data_3_LEFT_RIGHT_BOTH(sids, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path);

%%
display_avg_data_3_LEFT_RIGHT_BOTH(sids, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path);    

%%
display_avg_data_3_LEFT_RIGHT_BOTH_LED_mod(sids, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path);    

%% generate Vm vs. yaw plot for the stimulation period
CALC_PSTH = 1;

generate_odor_evoked_Vm_vs_yaw_plot_LAL_DNs(sids, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH);

%% Generate expected turning vs. ignoring stimulus data
turn_metadata = generate_turning_metadata_v2( sids, t_vel_all, yaw_vel_all,  analysis_path );

[ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list_v3( t_vel_all, forward_vel_all, turn_metadata );

%% Sort trials by Moving vs. stationary 

[ condition_trials, condition_trials_str, condition_str ] = generate_moving_vs_not_moving_trial_list( t_vel_all, forward_vel_all );

%% Turning vs. ignoring voltage averages
asid = 0;

avg_cond_btrace_trace_filepath = [ analysis_path '/ephys_' condition_str '_asid_' num2str( asid ) '_sid_' num2str(sids) ];

%display_two_condition_trials_ephys( condition_trials, condition_trials_str, t_volts, volts_all, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath);

MAX_TRIALS = 300;
display_two_condition_trials_avg_ephys_v2( condition_trials, condition_trials_str, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, avg_cond_btrace_trace_filepath, MAX_TRIALS);

%% Turning vs. ignoring PSTH averages
asid = 0;

avg_cond_btrace_trace_filepath = [ analysis_path '/psth_' condition_str '_asid_' num2str( asid ) '_sid_' num2str(sids) ];

%display_two_condition_trials_psth( condition_trials, condition_trials_str, t_vel_all, volt_PSTH, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, 1);

display_two_condition_trials_avg_psth_v2( condition_trials, condition_trials_str, t_vel_all, volt_PSTH, t_vel_all, yaw_vel_all, forward_vel_all,  avg_cond_btrace_trace_filepath, 1);

%% Spike sort
tt = 1;
trial_cnt = 1;
cur_volt = squeeze(volts_all{ tt }( trial_cnt, : ));
t = t_volts{ tt }( trial_cnt, : );

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

psth_dt_samples = ephys_SR/ball_SR;
psth_dt = psth_dt_samples / (1.0*ephys_SR);

VmFilt = medfilt1(cur_volt, 0.08 * ephys_SR, 'truncate');        
delta_Vm = cur_volt - VmFilt;

figure;
subplot(3,1,1);
hold on;
plot(t, cur_volt, 'b');
plot(t, VmFilt, 'g');

subplot(3,1,2);

d = diff(delta_Vm);

volts_thresh = d;
volts_thresh(find(d < 0.75)) = 0;

[~, locs] = findpeaks(volts_thresh, 'MinPeakProminence',0.02, 'Annotate','extents');

hold on;
plot(t(2:end), d, 'b');
plot(t(locs), d(locs), 'xg');

subplot(3,1,3);

spikes = zeros(1, length(delta_Vm));
spikes(locs+1) = 1;

cur_psth = sum(reshape(spikes, [psth_dt_samples, size(volts_all{ tt },2)/psth_dt_samples]),1) / psth_dt;
cur_psth = hanningsmooth(cur_psth, 50);

hold on;
plot(t_vel_all{tt}(1,:), cur_psth);

%% Plot raster for each trial along with yaw velocity
ac = get_analysis_constants();

first_stim_t = 3.0;
last_stim_t =  3.5;

GAP_LINES = 3;

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

for tt=1:3
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    
    BLOCK_SIZE = 2 + GAP_LINES;
    
    image_volt_yaw_fwd = zeros( BLOCK_SIZE * size(yaw_vel_all{tt},1), size(yaw_vel_all{tt},2) );
    
    trial_cnt = 1;
    for tr = 1:BLOCK_SIZE:size(yaw_vel_all{tt},1)*BLOCK_SIZE
    
        %image_volt_yaw_fwd( tr, : ) = forward_vel_all{ tt }( trial_cnt, : );
        image_volt_yaw_fwd( tr, : ) = yaw_vel_all{ tt }( trial_cnt, : );
        
        % downsample t_volts
        dt = ephys_SR / ball_SR;
        cur_volt = squeeze(volts_all{ tt }( trial_cnt, : ));               
        
        t_volts_down = squeeze(mean(reshape(cur_volt, [dt length(cur_volt)/dt])));
        
        VmFilt = medfilt1(t_volts_down, 0.08 * ball_SR, 'truncate');
        
        delta_Vm = t_volts_down-VmFilt;
                
%         if (tt==1 && tr == 1)
%             figure;
%             plot(delta_Vm);
%             figure(f);
%         end
        
        THRESHOLD = 0.15;
        delta_Vm(find(delta_Vm < THRESHOLD)) = 0;

        ff = figure;
        hold on;
        plot(t_volts{1}(1,:), cur_volt-mean(cur_volt), 'color', 'g');
        plot(t_vel_all{1}( 1, : ), delta_Vm);
        axis tight;
        waitforbuttonpress;
        close(ff);
        figure(f);
        
        
        image_volt_yaw_fwd( tr+1, : ) = delta_Vm;
        image_volt_yaw_fwd( tr+2, : ) = delta_Vm;
        image_volt_yaw_fwd( tr+3, : ) = delta_Vm;
    
        trial_cnt = trial_cnt + 1;
    end
    
    %subplot(1,2,tt);
    imagesc('XData', t_vel_all{tt}(1,:), 'YData', [0:size(image_volt_yaw_fwd,1)], 'CData', image_volt_yaw_fwd);
    colormap polarmap;
    caxis([-0.1 0.1]);
    axis tight;
    colorbar    
    xlabel( 'Time (s)' );
    title(ac.task_str{tt});
    
    raster_file = [analysis_path 'raster_and_yaw_sid_ ' num2str(sids) ' _' ac.task_str{tt} '.fig'];
    saveas(f, raster_file);
end

%% Calculate PSTH for a voltage trace, LAL-DNs
volt_PSTH = calculate_PSTH_for_LAL_DN( t_volts, t_vel_all, volts_all );

%% Calculate PSTH for a voltage trace SMP to LAL
volt_PSTH = cell(2,1);

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

psth_dt_samples = ephys_SR/ball_SR;
psth_dt = psth_dt_samples / (1.0*ephys_SR);

for tt=1:2
    
    volt_PSTH{tt} = zeros(size(volts_all{ tt },1), size(volts_all{ tt },2)/psth_dt_samples); 
    
    for tr=1:size(volts_all{ tt },1)
        
        %%%
        % Get spikes from a voltage trace (current clamp)
        %%%
        cur_volt = squeeze(volts_all{ tt }( tr, : ));
        VmFilt = medfilt1(cur_volt, 0.08 * ephys_SR, 'truncate');        
        delta_Vm = cur_volt - VmFilt;

        d = diff([delta_Vm(1) delta_Vm],1);
        
        volts_thresh = d;
        volts_thresh(find(d < 0.15)) = 0;        
                
        if 0
            ff = figure;
            hold on;
            
            ax1 = subplot(2,1,1);            
            plot(t_volts{tt}(1,1:end), cur_volt);            
            
            ax2 = subplot(2,1,2);            
            plot(t_volts{tt}(1,1:end), volts_thresh);            
            linkaxes([ax1 ax2],'x')
            
            break;
            %waitforbuttonpress;
            %close(ff);
        end
        
        if 0
            ff = figure;
            hold on;
            plot(t_volts{tt}(1,2:end), d);            
            
            waitforbuttonpress;
            close(ff);
        end
        
        
        %[~, locs] = findpeaks(volts_thresh, 'MinPeakProminence', 0.3, 'Annotate','extents');
        [~, locs] = findpeaks(volts_thresh, 'MinPeakDistance', 0.005*ephys_SR, 'Annotate','extents');
        
        spikes = zeros(1, length(delta_Vm));
        spikes(locs+1) = 1;
        %%%%%%%%%%%%%%
        
        if 0
            ff = figure;
            hold on;
            plot(t_volts{tt}(1,1:end), volts_thresh);            
            plot(t_volts{tt}(1,1:end), volts_thresh .* spikes(1:end), 'x', 'color', 'g', 'LineWidth', 2.0);
            
            %waitforbuttonpress;
            %close(ff);
            break;
        end        
        
        %%%
        % Get spikes from a voltage trace (current clamp)
        %%%
        cur_psth = sum(reshape(spikes, [psth_dt_samples, size(volts_all{ tt },2)/psth_dt_samples]),1) / psth_dt;
        cur_psth = hanningsmooth(cur_psth, 30);
        %%%%%%%%%%%%%%
        
        if 0
            ff = figure;
            subplot(2,1,1);
            hold on;
            plot(t_volts{tt}(1,:), volts_all{ tt }(tr,:));
            plot(t_volts{tt}(1,:), volts_all{ tt }(tr,:) .* spikes, 'x', 'color', 'g');
        
            subplot(2,1,2);
            plot(t_vel_all{tt}( 1, : ), cur_psth);

%             waitforbuttonpress;
%             close(ff);            
            break;
        end        
        volt_PSTH{tt}(tr,:) = cur_psth;
    end
    
    %break;
end


%% Calculate PSTH for a current trace
volt_PSTH = cell(2,1);

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

psth_dt_samples = ephys_SR/ball_SR;
psth_dt = psth_dt_samples / (1.0*ephys_SR);

for tt=1:2
    
    current_PSTH{tt} = zeros(size(current_all{ tt },1), size(current_all{ tt },2)/psth_dt_samples); 
    
    for tr=1:size(current_all{ tt },1)
        
        %%%
        % Get spikes from a voltage trace (current clamp)
        %%%
        cur_current = squeeze(current_all{ tt }( tr, : ));
        
         if 0
            ff = figure;
            hold on;
            CurFilt = medfilt1(cur_current, 100, 'truncate'); 
            plot(t_volts{tt}(1,1:end), CurFilt);            
            
            waitforbuttonpress;
            close(ff);
         end       
        
        CurFilt = medfilt1(cur_current, 0.08 * ephys_SR, 'truncate');        
        delta_cur = cur_current - CurFilt;

        d = diff(delta_cur);
        
        volts_thresh = d;
        volts_thresh(find(d < 0.3)) = 0;

        if 1
            ff = figure;
            hold on;
            plot(t_volts{tt}(1,2:end), d);            
            
            waitforbuttonpress;
            close(ff);
        end
        
        
        [~, locs] = findpeaks(volts_thresh, 'MinPeakProminence',0.02, 'Annotate','extents');
        
        spikes = zeros(1, length(delta_Vm));
        spikes(locs+1) = 1;
        %%%%%%%%%%%%%%
        
        
        %%%
        % Get spikes from a voltage trace (current clamp)
        %%%
        cur_psth = sum(reshape(spikes, [psth_dt_samples, size(volts_all{ tt },2)/psth_dt_samples]),1) / psth_dt;
        cur_psth = hanningsmooth(cur_psth, 30);
        %%%%%%%%%%%%%%
        
        if 1
            ff = figure;
            subplot(2,1,1);
            hold on;
            plot(t_volts{tt}(1,:), volts_all{ tt }(tr,:));
            plot(t_volts{tt}(1,:), volts_all{ tt }(tr,:) .* spikes, 'x', 'color', 'g');
        
            subplot(2,1,2);
            plot(t_vel_all{tt}( 1, : ), cur_psth);

            %waitforbuttonpress;
            %close(ff);
            break;
        end
        
        volt_PSTH{tt}(tr,:) = cur_psth;
    end
    
end

%% Plot PSTH for each trial along with yaw velocity
ac = get_analysis_constants();

first_stim_t = 3.0;
last_stim_t =  3.5;

GAP_LINES = 3;

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

for tt=1:2
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    
    BLOCK_SIZE = 2 + GAP_LINES;
    
    image_psth_yaw_fwd = zeros( BLOCK_SIZE * size(yaw_vel_all{tt},1), size(yaw_vel_all{tt},2) );
    
    trial_cnt = 1;
    for tr = 1:BLOCK_SIZE:size(yaw_vel_all{tt},1)*BLOCK_SIZE
    
        %image_volt_yaw_fwd( tr, : ) = forward_vel_all{ tt }( trial_cnt, : );
        image_psth_yaw_fwd( tr, : ) = yaw_vel_all{ tt }( trial_cnt, : );
                        
        image_psth_yaw_fwd( tr+1, : ) = volt_PSTH{tt}(trial_cnt,:)./600.0;
    
        trial_cnt = trial_cnt + 1;
    end
    
    imagesc('XData', t_vel_all{tt}(1,:), 'YData', [0:size(image_psth_yaw_fwd,1)], 'CData', image_psth_yaw_fwd);
    colormap polarmap;
    caxis([-0.15 0.15]);
    axis tight;
    colorbar    
    xlabel( 'Time (s)' );
    title(ac.task_str{tt});
    
    raster_file = [analysis_path 'psth_and_yaw_sid_ ' num2str(sids) ' _' ac.task_str{tt} '.fig'];
    saveas(f, raster_file);
end

%% Plot PSTH vs. yaw and fwd vel for each trial
ac = get_analysis_constants();
settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;
Fs = ephys_SR;
Fs_vel = ball_SR;

first_stim_t = 3.0;
last_stim_t =  3.5;

for tt=1:2
    
    for tr = 1:size(yaw_vel_all{tt},1)
        
        f = figure;
        
        trial_name = [ ac.task_str{tt} '_sid_ ' num2str(sids) '_tr_' num2str(tr)];
        
        subplot(3,1,1);
        hold on;
        plot(t_vel_all{tt}(1,:),  forward_vel_all{tt}(tr, :));
        xlim([0 t_vel_all{tt}(1,end)]);
        ylim([-0.4 0.4]);
        ylabel('Fwd vel (mm/s)');
        
        l = title( trial_name );
        set(l, 'Interpreter', 'none')
        
        yy = ylim;
        y_min = yy(1); y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');
        
        subplot(3,1,2);
        hold on;
        plot(t_vel_all{tt}(1,:), yaw_vel_all{tt}( tr, : ));
        xlim([0 t_vel_all{tt}(1,end)]);
        ylabel('Yaw vel (deg/s)');
        ylim([-0.4 0.4]);
        
        yy = ylim;
        y_min = yy(1); y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');
        
        subplot(3,1,3);
        hold on;
        plot(t_vel_all{tt}(1,:), volt_PSTH{tt}(tr,:));
        xlim([0 t_vel_all{tt}(1,end)]);
        ylim([0 120]);
        
        ylabel('Firing rate (Spikes/s)');
        xlabel('Time (s)');
        
        yy = ylim;
        y_min = yy(1); y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');
        
        saveas(f, [analysis_tt_path '/psdh_vs_behavior_' trial_name '.fig']);
        saveas(f, [analysis_tt_path '/psdh_vs_behavior_' trial_name '.png']);
        close(f);
    end
end


%% Show wind stim
wind_on = 2.0;
wind_off = 5.5;
yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ wind_on wind_on wind_off wind_off ],[y_min y_max y_max y_min ], rgb('Lavender'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

%% Show odor stim
odor_on = 3.0;
odor_off = 3.5;
yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ odor_on odor_on odor_off odor_off ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

%%
avidir = '/Users/sasha/Dropbox/Wilson_lab/data/160719_ablation_gfp_vt043158_04/analysis/stack_.avi';
vidObj=VideoReader(avidir);

i = 1;
clear vid;
while hasFrame(vidObj)
    vidFrame = rgb2gray(readFrame(vidObj));
    vid(:,:,i) = vidFrame;
    i = i + 1;
end

%% Generate expected turning vs. ignoring stimulus data
[bdata_vel_time, bdata_vel] = reformat_raw_behavioral_data_ephys( bdata_time, bdata_raw, datapath );
bdata_vel_time = bdata_vel_time - TIME_OFFSET;

turn_metadata = generate_turning_metadata( sids, bdata_vel_time, bdata_vel, analysis_path );

[ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list_v2( sids, bdata_vel_time, bdata_vel, turn_metadata, analysis_path );
