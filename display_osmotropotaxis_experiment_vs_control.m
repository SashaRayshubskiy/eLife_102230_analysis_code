function display_osmotropotaxis_experiment_vs_control( t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_dir )

% Format of inputs: cell_array 
% { experiment type, first arg exp=1,control=2 }{directory}{trial_type}(trial,data)
% 

first_stim_t = 3.0;
last_stim_t =  3.5;

ac = get_analysis_constants();

f = figure;

ax1 = subplot(3,1,1);
hold on;
t = squeeze(t_vel_all{1}{1}{1}(1,:));

for tt = 1:3
    
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
        
    for e_type = 1:2
        all_forward_vel = [];
        
        for d_idx = 1:length(t_vel_all{e_type})
            all_forward_vel = vertcat(all_forward_vel, forward_vel_all{e_type}{d_idx}{tt}(:,:));
        end
    
        avg_vel_fwd = squeeze(mean(all_forward_vel));
        
        std_fwd_trace = squeeze(std(all_forward_vel,1 ))/ sqrt(size(all_forward_vel,1));
        
        fh = fill([ t fliplr(t) ], [(avg_vel_fwd+std_fwd_trace) fliplr((avg_vel_fwd-std_fwd_trace)) ], cur_color_single );
        set(fh, 'EdgeColor', 'None');
        
        lineStyle = '-';
        if(e_type == 2)
            lineStyle = '--';
        end
        
        pp(e_type, tt) = plot(t,avg_vel_fwd, 'color', cur_color_avg, 'LineWidth', 2, 'LineStyle', lineStyle );
        trial_cnts(e_type, tt) = size(all_forward_vel,1);
    end
    
    xlim([0 t(end)]);
    ylabel('Fwd vel (mm/s)');
end

legend( [ pp(1,1), pp(1,2), pp(1,3), pp(2,1), pp(2,2), pp(2,3) ], ...
        ['Exp Left trials (' num2str(trial_cnts(1, 1)) ')'], ...
        ['Exp Right trials (' num2str(trial_cnts(1, 2)) ')'], ...
        ['Exp Both trials (' num2str(trial_cnts(1, 3)) ')'], ...
        ['Control Left trials (' num2str(trial_cnts(2, 1)) ')'], ...
        ['Control Right trials (' num2str(trial_cnts(2, 2)) ')'], ...
        ['Control Both trials (' num2str(trial_cnts(2, 3)) ')'] );  
   
set(gca,'fontsize',14)

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

ax2 = subplot(3,1,2);
hold on;
for tt = 1:3
    
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
    
    for e_type = 1:2
        all_yaw_vel = [];
        
        for d_idx = 1:length(t_vel_all{e_type})
            all_yaw_vel = vertcat(all_yaw_vel, yaw_vel_all{e_type}{d_idx}{tt}(:,:));
        end
        
        avg_vel_yaw = squeeze(mean(all_yaw_vel));
        
        std_yaw_trace = squeeze(std(all_yaw_vel,1 ))/ sqrt(size(all_yaw_vel,1));
        
        fh = fill([ t fliplr(t) ], [(avg_vel_yaw+std_yaw_trace) fliplr((avg_vel_yaw-std_yaw_trace)) ], cur_color_single );
        set(fh, 'EdgeColor', 'None');
        
        lineStyle = '-';
        if(e_type == 2)
            lineStyle = '--';
        end
        
        pp(e_type, tt) = plot(t,avg_vel_yaw, 'color', cur_color_avg, 'LineWidth', 2, 'LineStyle', lineStyle );
        trial_cnts(e_type, tt) = size(all_yaw_vel,1);
    end
    
    xlim([0 t(end)]);
    ylabel('Yaw vel (deg/s)');
end
set(gca,'fontsize',14)

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

ax3 = subplot(3,1,3);

t = squeeze(t_volts{1}{1}{1}(1,:));

if 1
hold on;
for tt = 1:3
    
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
    
    for e_type = 1:2
        all_volts = [];
        
        for d_idx = 1:length(t_vel_all{e_type})
            all_volts = vertcat(all_volts, volts_all{e_type}{d_idx}{tt}(:,:));
        end
        
        avg_vols = squeeze(mean(all_volts));
        baseline = mean(avg_vols);
        avg_vols_baseline_corr = avg_vols - repmat(baseline, [1 length(avg_vols)]);
        
        std_volt_trace = squeeze(std( all_volts,1 ))/ sqrt(size(all_volts,1));                
               
        fh = fill([ t fliplr(t) ], [(avg_vols_baseline_corr+std_volt_trace) fliplr((avg_vols_baseline_corr-std_volt_trace)) ], cur_color_single );
        set(fh, 'EdgeColor', 'None');
        
        lineStyle = '-';
        if(e_type == 2)
            lineStyle = '--';
        end
        
        pp(e_type, tt) = plot(t, avg_vols_baseline_corr, 'color', cur_color_avg, 'LineWidth', 2, 'LineStyle', lineStyle );
    end

    xlim([0 t(end)]);
    ylim([-2 3])
    ylabel('delta volts (mV)');
    xlabel('Time (s)');
end
set(gca,'fontsize',14)

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');
end

linkaxes([ax1 ax2 ax3], 'x');

saveas(f, [analysis_dir '/osmotropotaxis_vs_control_behavior_and_ephys_avg.fig']);
saveas(f, [analysis_dir '/osmotropotaxis_vs_control_behavior_and_ephys_avg.png']);




end

