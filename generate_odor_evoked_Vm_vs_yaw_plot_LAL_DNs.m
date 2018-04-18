function [ output_args ] = generate_odor_evoked_Vm_vs_yaw_plot_LAL_DNs(sids, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH)

file_str = 'Vm';
WINDOW_EPHYS_AVG = 0.2;
if (CALC_PSTH == 1)
    volts_PSTH = calculate_PSTH_for_LAL_DN(  t_volts, t_vel_all, volts_all );
    file_str = 'PSTH';
    WINDOW_EPHYS_AVG = 0.3;
end

yaw_per_trial = {};
ephys_per_trial = {};
fwd_vels = [];

settings = sensor_settings;
stim_start_t = settings.pre_stim;
stim_end_t = settings.pre_stim + settings.stim;

volt_avg_win_start_t = stim_start_t + 0.050;
volt_avg_win_end_t = volt_avg_win_start_t + WINDOW_EPHYS_AVG;

yaw_avg_win_start_t = stim_start_t + 0.150;
yaw_avg_win_end_t = yaw_avg_win_start_t + 0.3;

SHOW_BOTH = 0;
FWD_CUTOFF = 0.12;

ac = get_analysis_constants();

WIN_BEFORE_STIM_ONSET = 1.0;

for tt = 1:length(t_volts)
    
    cur_trial_cnt = 1;
    
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
        
        if(CALC_PSTH == 1)
            cur_psth = volts_PSTH{tt}(trial,:);
            [~,psth_avg_win_t] = find( (t_vel_all{tt} >= volt_avg_win_start_t) & (t_vel_all{tt} <= volt_avg_win_end_t) );
            volt = mean( cur_psth(psth_avg_win_t) );            
        else            
            cur_volt = volts_all{tt}(trial,:);
            
            [~,volt_baseline_t] = find( (t_volts{tt} >= (volt_avg_win_start_t-1.0)) & (t_volts{tt} <= volt_avg_win_start_t) );
            volt_baseline = mean(cur_volt(volt_baseline_t));
            
            cur_volt_baseline = cur_volt - volt_baseline;
            
            [~,volt_avg_win_t] = find( (t_volts{tt} >= volt_avg_win_start_t) & (t_volts{tt} <= volt_avg_win_end_t) );
            volt = mean( cur_volt_baseline(volt_avg_win_t) );
        end
        
        cur_yaw = yaw_vel_all{tt}(trial,:);
        [~,yaw_avg_win_t] = find( (t_vel_all{tt} >= yaw_avg_win_start_t) & (t_vel_all{tt} <= yaw_avg_win_end_t) );
        
        % Try 1. Just averate the yaw in the stim window 
        yaw = mean( cur_yaw(yaw_avg_win_t) );
               
        % Try 3. findpeaks();
        % [~, locs] = findpeaks( cur_yaw(yaw_avg_win_t), 'NPeaks', 1, 'SortStr', 'descend', 'Annotate','extents');
        % avg_yaw = cur_yaw(locs(1));
         
        yaw_per_trial{tt}(cur_trial_cnt) = yaw;
        ephys_per_trial{tt}(cur_trial_cnt) = volt;     
        cur_trial_cnt = cur_trial_cnt + 1;
    end
end

if 0
    % Show fwd velocity histogram
    f1 = figure;
    hist(fwd_vels,100000);
    ylabel('Counts');
    xlabel('Fwd velocity');
    xlim([-0.5 0.5]);
end

f = figure;

for tt = 1:length(t_volts)
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
   
    if(( tt == ac.BOTH ) && (SHOW_BOTH == 0)) 
        continue;
    end
    
    hold on;
    if(( tt == ac.BOTH ) || ( tt == ac.RIGHT ))
        sct_hdl{tt}(1) = scatter( yaw_per_trial{tt}, ephys_per_trial{tt}, [], cur_color_avg );        
    else
        sct_hdl{tt}(1) = scatter( yaw_per_trial{tt}, ephys_per_trial{tt}, [], cur_color_avg, 'filled' );
    end
    
    if(CALC_PSTH == 1)
        ylabel('Firing rate (spikes/s)');
    else
        ylabel('delta Volt (mV)');        
    end
    
    xlabel('Yaw velocity (deg/s)');    
end

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

title(['FWD vel cutoff: ' num2str(FWD_CUTOFF) '  using avg for yaw during stim']);

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
   
save([analysis_path '/odor_evoked_' file_str '_metadata.mat'], 'yaw_per_trial', 'ephys_per_trial');

timenow_str = datestr(datetime, 'yymmdd_HHMMSS');
saveas(f, [analysis_path '/odor_evoked_' timenow_str '_' file_str '_vs_yaw.fig' ]);
saveas(f, [analysis_path '/odor_evoked_' timenow_str '_' file_str '_vs_yaw.png' ]);

Nbins = 30;

f2 = figure;
hold on;
h1 = histogram(yaw_per_trial{1}, Nbins, 'FaceColor',  ac.LEFT_CLR ,'EdgeColor', 'none');
hold on;
h2 = histogram(yaw_per_trial{2}, Nbins, 'FaceColor',  ac.RIGHT_CLR ,'EdgeColor', 'none');
xlabel('Yaw vel (deg/s)');
ylabel('Count');
legend([h1,h2],'Left Odor','Right Odor');
saveas(f2, [analysis_path '/odor_evoked_' timenow_str '_yaw_hist.fig' ]);
saveas(f2, [analysis_path '/odor_evoked_' timenow_str '_yaw_hist.png' ]);

f3 = figure;
hold on;
h1 = histogram(ephys_per_trial{1}, Nbins, 'FaceColor',  ac.LEFT_CLR ,'EdgeColor', 'none');
hold on;
h2 = histogram(ephys_per_trial{2}, Nbins, 'FaceColor',  ac.RIGHT_CLR ,'EdgeColor', 'none');
 
if(CALC_PSTH == 1)
    xlabel('Firing rate (spikes/s)');
else
    xlabel('delta Vm (mV)');
end

ylabel('Count');
legend([h1,h2],'Left Odor','Right Odor');
saveas(f3, [analysis_path '/odor_evoked_' timenow_str '_ephys_hist.fig' ]);
saveas(f3, [analysis_path '/odor_evoked_' timenow_str '_ephys_hist.png' ]);

end