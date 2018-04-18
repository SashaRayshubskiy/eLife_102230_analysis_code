function [ image_cnts ] = generate_Vm_vs_yaw_trajectory_LAL_DNs(sids, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH)

file_str = 'Vm';
WINDOW_EPHYS_AVG = 0.2;
if (CALC_PSTH == 1)
    volts_PSTH = calculate_PSTH_for_LAL_DN( t_volts, t_vel_all, volts_all );
    file_str = 'PSTH';
    WINDOW_EPHYS_AVG = 0.3;
end

OVERLAY_ODOR_EVOKED_DATA = 1;

yaw_per_trial = {};
ephys_per_trial = {};
fwd_vels = [];

settings = sensor_settings;
stim_start_t = settings.pre_stim;
stim_end_t = settings.pre_stim + settings.stim;

volt_avg_win_start_t = stim_start_t + 0.050;
volt_avg_win_end_t = volt_avg_win_start_t + WINDOW_EPHYS_AVG;

yaw_avg_win_start_t = stim_start_t + 0.050;
yaw_avg_win_end_t = stim_end_t;

SHOW_BOTH = 1;
FWD_CUTOFF = 0.12;

ac = get_analysis_constants();

f = figure;

cur_trial_cnt = 1;

cm = colormap(jet(575));

WIN_BEFORE_STIM_ONSET = 1.0;

% Image X: [0 250], Y: [-2000 2000], 
image_cnts = zeros(250, 4000);

for tt = 1:length(t_volts)
%for tt = 1:2
    
    if(tt == ac.LEFT)
        cur_color_single = ac.LEFT_CLR_SINGLE;
        cur_color_avg = ac.LEFT_CLR;
    elseif(tt == ac.RIGHT)
        cur_color_single = ac.RIGHT_CLR_SINGLE;
        cur_color_avg = ac.RIGHT_CLR;
    elseif(tt == ac.BOTH)
        cur_color_single = ac.BOTH_CLR_SINGLE;
        cur_color_avg = ac.BOTH_CLR;
    end
    
    for trial = 1:size(t_volts{tt},1)

        % Throw away trials where the fly is not moving during the
        % stimulation period
        
        cur_fwd = forward_vel_all{tt}(trial,:);
        [ ~, fwd_avg_win_t ] = find( (t_vel_all{tt} >= stim_start_t-WIN_BEFORE_STIM_ONSET) & (t_vel_all{tt} <= stim_start_t) );
        
        % Try 1. Just averate the yaw in the stim window 
        avg_fwd = mean( cur_fwd(fwd_avg_win_t) );
        
        fwd_vels(end+1) = avg_fwd;
        
        if (avg_fwd < FWD_CUTOFF)
            continue;
        end        
        
        DT = 5;
        
        % downsample yaw        
        cur_yaw = yaw_vel_all{tt}(trial,:);
        
        yaw_t_down = squeeze(mean(reshape(t_vel_all{tt}(1,:), [DT, length(cur_yaw)/DT]),1));
        cur_yaw_down = squeeze(mean(reshape(cur_yaw, [DT, length(cur_yaw)/DT]),1));
        
        [ ~, yaw_baseline_t ] = find( ( yaw_t_down <= stim_end_t ));
        [ ~, yaw_before_stim_t ] = find( yaw_t_down <= stim_start_t );
        [ ~, yaw_during_stim_t ] = find( (yaw_t_down >= stim_start_t ) & ((yaw_t_down <= stim_end_t )));
        
        if(CALC_PSTH == 1)
            cur_psth = volts_PSTH{tt}(trial,:);
            
            cur_psth_down = squeeze(mean(reshape(cur_psth, [DT, length(cur_psth)/DT]),1));
            
        else            
            % downsample volt
            dt_v = 40;
            DT_v = dt_v*DT;
            
            cur_volt = volts_all{tt}(trial,:) - mean(volts_all{tt}(trial,:));
            
            volt_t_down = squeeze(mean(reshape(t_volts{tt}(1,:), [DT_v, length(cur_volt)/DT_v]),1));
            cur_volt_down = squeeze(mean(reshape(cur_volt, [DT_v, length(cur_volt)/DT_v]),1));
            
            [ ~, volt_baseline_t ] = find( ( volt_t_down <= stim_end_t ));
            [ ~, volt_before_stim_t ] = find( volt_t_down <= stim_start_t );
            [ ~, volt_during_stim_t ] = find( (volt_t_down >= stim_start_t ) & ((volt_t_down <= stim_end_t )));        
        end     
        
        shift_factor = 3;
        if(CALC_PSTH == 1)
            hold on;
            
            cur_window_t = yaw_baseline_t;
            
            act_yaw = cur_yaw_down(cur_window_t);
            act_psth = cur_psth_down(cur_window_t);
            act_t = yaw_t_down(cur_window_t);
            
            % z = zeros(size(act_yaw));
            % colormap jet;
            % S = surface([act_yaw;act_yaw], [act_psth;act_psth], [z;z], [act_t;act_t], 'facecol', 'no', 'edgecol', 'interp', 'linew', 2);
            % set(S, 'EdgeAlpha', 0.1);
                        
            cur_yaw_plot = cur_yaw_down(shift_factor:end);
            cur_psth_plot = cur_psth_down(1:end-shift_factor+1);
            
            %plot(cur_yaw_plot, cur_psth_plot, 'o', 'MarkerSize', 3, 'Color', rgb('Maroon') );
            plot(cur_yaw_plot, cur_psth_plot, 'o', 'MarkerSize', 3, 'Color', rgb('DarkGray') );
                        
            % Digitize the FR vs. Yaw plot and return 
            for ii = 1:length(cur_yaw_plot)
                cur_y_idx = floor( cur_yaw_plot( ii ) ) + 2000; % range is [-2000 2000] => [0 4000]
                cur_x_idx = floor( cur_psth_plot( ii ) ) + 1;
                
                image_cnts(cur_x_idx, cur_y_idx) = image_cnts(cur_x_idx, cur_y_idx) + 1;
            end
        else
            hold on;
            
            plot(cur_yaw_down(shift_factor:end), cur_volt_down(1:end-shift_factor+1), 'o', 'MarkerSize', 3, 'Color', rgb('Maroon') );            
        end
        
        cur_trial_cnt = cur_trial_cnt + 1;
        %waitforbuttonpress;
    end
end

if (CALC_PSTH == 1 )
    ylabel('Firing rate (spikes/s)');
    title(['FWD vel cutoff: ' num2str(FWD_CUTOFF) '  FR vs. yaw for ( ' num2str(cur_trial_cnt) ' ) trials - full trial, window: ' num2str(DT/settings.sensorPollFreq) ]);
else
    ylabel('delta Voltage (mV)');    
    title(['FWD vel cutoff: ' num2str(FWD_CUTOFF) '  Vm vs. yaw for ( ' num2str(cur_trial_cnt) ' ) trials - full trial, window: ' num2str(DT/settings.sensorPollFreq) ]);
end

xlabel('Yaw (deg/s)');


if 0
    % Show fwd velocity histogram
    f1 = figure;
    hist(fwd_vels,100000);
    ylabel('Counts');
    xlabel('Fwd velocity');
    xlim([-0.5 0.5]);
end

if 0
if (CALC_PSTH == 0 )
    ACTION_TAKEN_LIMIT = 50;
    CHOICE_MADE_LIMIT = 1.0;
    
    yy = ylim;
    y_min = yy(1); y_max = yy(2);
    hh = fill([ -1.0*ACTION_TAKEN_LIMIT -1.0*ACTION_TAKEN_LIMIT ACTION_TAKEN_LIMIT ACTION_TAKEN_LIMIT ],[y_min y_max y_max y_min ], rgb('Wheat'));
    set(gca,'children',circshift(get(gca,'children'),-1));
    set(hh, 'EdgeColor', 'None');
    
    xx = xlim;
    x_min = xx(1); x_max = xx(2);
    hh = fill( [ x_min x_min x_max x_max ], [ -1.0*CHOICE_MADE_LIMIT CHOICE_MADE_LIMIT CHOICE_MADE_LIMIT -1.0*CHOICE_MADE_LIMIT ], rgb('Lavender'));
    set(gca,'children',circshift(get(gca,'children'),-1));
    set(hh, 'EdgeColor', 'None');
end
end

title(['FWD vel cutoff: ' num2str(FWD_CUTOFF) ' for ' num2str(WIN_BEFORE_STIM_ONSET) ' before stim']);

if(OVERLAY_ODOR_EVOKED_DATA == 1)
    
    if(CALC_PSTH == 1)        
        oed_datapath = [analysis_path '/odor_evoked_PSTH_metadata.mat'];
    else
        oed_datapath = [analysis_path '/odor_evoked_Vm_metadata.mat'];        
    end
    
    oed = load(oed_datapath);
    
    yaw_per_trial   = oed.yaw_per_trial
    ephys_per_trial = oed.ephys_per_trial;
    
    hold on;
    ac = get_analysis_constants();
    sct_hdl(1) = scatter( yaw_per_trial{1}, ephys_per_trial{1}, [], ac.LEFT_CLR, 'filled', 'LineWidth', 0.4 );        
    sct_hdl(2) = scatter( yaw_per_trial{2}, ephys_per_trial{2}, [], ac.RIGHT_CLR, 'filled', 'LineWidth', 0.4 );        
    legend([sct_hdl(1), sct_hdl(2) ], ...
        ['Left Odor ( ' num2str(size(ephys_per_trial{1},2)) ' )'], ...
        ['Right Odor ( ' num2str(size(ephys_per_trial{2},2)) ' )'] );
end

if 0
if (SHOW_BOTH == 1)
    legend([sct_hdl{1}(1), sct_hdl{2}(1), sct_hdl{3}(1)], ...
        ['Left Odor ( ' num2str(size(ephys_per_trial{1},2)) ' )'], ...
        ['Right Odor ( ' num2str(size(ephys_per_trial{2},2)) ' )'], ...
        ['Both Odor ( ' num2str(size(ephys_per_trial{3},2)) ' )']);
else
    legend([sct_hdl{1}(1), sct_hdl{2}(1) ], ...
        ['Left Odor ( ' num2str(size(ephys_per_trial{1},2)) ' )'], ...
        ['Right Odor ( ' num2str(size(ephys_per_trial{2},2)) ' )'] );
end
end
   
set(gca, 'FontSize', 16);

saveas(f, [analysis_path '/trajectory_' file_str '_vs_yaw.fig' ]);
saveas(f, [analysis_path '/trajectory_' file_str '_vs_yaw.png' ]);

end