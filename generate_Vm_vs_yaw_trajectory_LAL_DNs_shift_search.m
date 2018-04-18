function [ r_sqr_value, shift_amount_in_sec ] = generate_Vm_vs_yaw_trajectory_LAL_DNs_shift_search(shift, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH, analysis_path_type)

file_str = 'Vm';
WINDOW_EPHYS_AVG = 0.2;
if (CALC_PSTH == 1)
    volts_PSTH = calculate_PSTH_for_LAL_DN( t_volts, t_vel_all, volts_all, analysis_path_type );
    file_str = 'PSTH';
    WINDOW_EPHYS_AVG = 0.3;
end

OVERLAY_ODOR_EVOKED_DATA = 0;

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

SHOW_BOTH = 0;

if(shift == -1)
    USE_RANDOM_YAW = 1;
else
    USE_RANDOM_YAW = 0;    
end

FWD_CUTOFF = 1.0;

ac = get_analysis_constants();

cur_trial_cnt = 1;

cm = colormap(jet(575));

WIN_BEFORE_STIM_ONSET = 1.0;

% Image X: [0 250], Y: [-2000 2000], 
if (CALC_PSTH == 1)
    image_cnts = zeros(250, 4000);
else
    image_cnts = zeros(28, 4000); % Voltage range: [-12 16]
end

yaw_data = [];
ephys_data = [];

for tt = 1:length(t_volts)

    yaw_stim_data{tt} = [];
    ephys_stim_data{tt} = [];

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
    
    tt_trial_idx = 1;
    
    for trial = 1:size(t_volts{tt},1)

        % Throw away trials where the fly is not moving during the
        % stimulation period        
        cur_fwd = forward_vel_all{tt}(trial,:);
        [ ~, fwd_avg_win_t ] = find( (t_vel_all{tt} >= stim_start_t-WIN_BEFORE_STIM_ONSET) & (t_vel_all{tt} <= stim_start_t) );        
        
        % Try 1. Just averate the yaw in the stim window 
        %avg_fwd = mean( cur_fwd(fwd_avg_win_t) );
        avg_fwd = mean( cur_fwd );
        
        fwd_vels(end+1) = avg_fwd;
        
        if (avg_fwd < FWD_CUTOFF)
            continue;
        end        
        
        DT = 2;
        TIME_PER_YAW_DATAPOINT_IN_SECONDS = 0.01;
        
        shift_amount_in_sec = DT * TIME_PER_YAW_DATAPOINT_IN_SECONDS;
        
        % downsample yaw 
        if( USE_RANDOM_YAW == 1 )
            rand_tt    = ceil(rand * length(t_volts));
            rand_trial = ceil(rand * size(t_volts{rand_tt},1));
            cur_yaw = yaw_vel_all{rand_tt}(rand_trial,:);
        else            
            cur_yaw = yaw_vel_all{tt}(trial,:);        
        end
        
        yaw_t_down = squeeze(mean(reshape(t_vel_all{tt}(1,:), [DT, length(cur_yaw)/DT]),1));
        cur_yaw_down = squeeze(mean(reshape(cur_yaw, [DT, length(cur_yaw)/DT]),1));
        cur_fwd_down = squeeze(mean(reshape(cur_fwd, [DT, length(cur_fwd)/DT]),1));
        
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
            
            cur_volt_down = squeeze(mean(reshape(cur_volt, [DT_v, length(cur_volt)/DT_v]),1));            
        end     
        
        if( USE_RANDOM_YAW == 1 )            
            shift_factor = 1;
        else
            shift_factor = shift;            
        end
                
        if(CALC_PSTH == 1)            
            cur_ephys_plot = cur_psth_down(1:end-shift_factor+1);
        else
            cur_ephys_plot = cur_volt_down(1:end-shift_factor+1);
        end
        
        yaw_keep = [];
        ephys_keep = [];
        t_keep = [];
        
        cur_fwd_plot  = cur_fwd_down(shift_factor:end);
        cur_yaw_plot  = cur_yaw_down(shift_factor:end);
        cur_t_vel     = yaw_t_down(shift_factor:end);
        
        for ii = 1:length(cur_fwd_plot)
            if (cur_fwd_plot(ii) < FWD_CUTOFF)
                continue;
            end
            
            yaw_keep(end+1)  = cur_yaw_plot( ii );
            ephys_keep(end+1) = cur_ephys_plot( ii );
            t_keep(end+1)    = cur_t_vel( ii );
        end
                
        yaw_data = horzcat(yaw_data, yaw_keep);
        ephys_data = horzcat(ephys_data, ephys_keep);

        %hold on;        
        %plot(yaw_keep, ephys_keep, 'o', 'MarkerSize', 3, 'Color', rgb('Maroon') );
                        
        cur_trial_cnt = cur_trial_cnt + 1;        
        tt_trial_idx = tt_trial_idx + 1;
    end
end

xlabel('Yaw (deg/s)');


%mdl = fitlm(yaw_data, ephys_data);
%hold on;
%plot(mdl);
%r_sqr_value = mdl.RSquared.Ordinary;

[curvefit,gof,output] = fit(yaw_data', ephys_data','exp2');
r_sqr_value = gof.rsquare;
hold on;
plot( curvefit );


title(['R^2: ' num2str(r_sqr_value)]);
xlim([-1500 1500]);
ylim([0 200]);

set(gca, 'FontSize', 16);

end