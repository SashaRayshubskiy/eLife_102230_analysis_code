function display_thermotaxis_experiment_by_fly( t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_dir ) 

first_stim_t = 3.0;
last_stim_t =  3.5;

ac = get_analysis_constants();

f = figure;

ax1 = subplot(3,1,1);
hold on;
t = squeeze(t_vel_all{1}{1}{1}(1,:));

for tt = 1:2
            
    for e_type = 1:1
        
        if(tt == ac.LEFT)
            cur_color_avg = ac.LEFT_CLR;
            cur_color_single = ac.LEFT_CLR_SINGLE;
        elseif(tt == ac.RIGHT)
            cur_color_avg = ac.RIGHT_CLR;
            cur_color_single = ac.RIGHT_CLR_SINGLE;
        end
                
        all_forward_vel = [];
        
        hold on;
        for d_idx = 1:length(t_vel_all{e_type})
            %all_forward_vel = vertcat(all_forward_vel, forward_vel_all{e_type}{d_idx}{tt}(:,:));
            all_forward_vel(d_idx,:) = mean(forward_vel_all{e_type}{d_idx}{tt}(:,:));
            plot(t, squeeze(all_forward_vel(d_idx,:)), 'color', cur_color_single );
        end
    
        avg_vel_fwd = squeeze(mean(all_forward_vel));
        
        std_fwd_trace = squeeze(std(all_forward_vel,1 ))/ sqrt(size(all_forward_vel,1));
        
%         fh = fill([ t fliplr(t) ], [(avg_vel_fwd+std_fwd_trace) fliplr((avg_vel_fwd-std_fwd_trace)) ], cur_color_single );
%         set(fh, 'EdgeColor', 'None');
        
        lineStyle = '-';
        if(e_type == 2)
            lineStyle = '--';
        end
        
        pp(e_type, tt) = plot(t,avg_vel_fwd, 'color', cur_color_avg, 'LineWidth', 2, 'LineStyle', lineStyle );
        trial_cnts(e_type, tt) = size(all_forward_vel,1);
    end
    
    xlim([0 t(end)]);
    ylim([-2 20]);
    ylabel('Fwd vel (mm/s)');
end

legend( [ pp(1,1), pp(1,2) ], ...
        ['Exp Left trials (' num2str(trial_cnts(1, 1)) ')'], ...
        ['Exp Right trials (' num2str(trial_cnts(1, 2)) ')'] );  
   
set(gca,'fontsize',14)

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

ax2 = subplot(3,1,2);
hold on;
for tt = 1:2
    
    for e_type = 1:1

        if(tt == ac.LEFT)
            cur_color_avg = ac.LEFT_CLR;
            cur_color_single = ac.LEFT_CLR_SINGLE;
        elseif(tt == ac.RIGHT)
            cur_color_avg = ac.RIGHT_CLR;
            cur_color_single = ac.RIGHT_CLR_SINGLE;
        end
      
        all_yaw_vel = [];
        
        for d_idx = 1:length(t_vel_all{e_type})
            %all_yaw_vel = vertcat(all_yaw_vel, yaw_vel_all{e_type}{d_idx}{tt}(:,:));
            tb = find((t > 2.5) & (t < 3.0) );
            
            baseline = squeeze(mean(mean( yaw_vel_all{e_type}{d_idx}{tt}(:,tb))));
            yavg = mean( yaw_vel_all{e_type}{d_idx}{tt}(:,:) );
            all_yaw_vel(d_idx,:) = yavg - baseline;
            
            plot(t, squeeze(all_yaw_vel(d_idx,:)), 'color', cur_color_single );
        end
        
        avg_vel_yaw = squeeze(mean(all_yaw_vel));
        
        std_yaw_trace = squeeze(std(all_yaw_vel,1 ))/ sqrt(size(all_yaw_vel,1));
        
%         fh = fill([ t fliplr(t) ], [(avg_vel_yaw+std_yaw_trace) fliplr((avg_vel_yaw-std_yaw_trace)) ], cur_color_single );
%         set(fh, 'EdgeColor', 'None');
        
        lineStyle = '-';
        if(e_type == 2)
            lineStyle = '--';
        end
        
        pp(e_type, tt) = plot(t,avg_vel_yaw, 'color', cur_color_avg, 'LineWidth', 2, 'LineStyle', lineStyle );
        trial_cnts(e_type, tt) = size(all_yaw_vel,1);
    end
    
    xlim([0 t(end)]);
    ylim([-70 70]);
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
for tt = 1:2
        
    for e_type = 1:1
      
        if(tt == ac.LEFT)
            cur_color_avg = ac.LEFT_CLR;
            cur_color_single = ac.LEFT_CLR_SINGLE;
        elseif(tt == ac.RIGHT)
            cur_color_avg = ac.RIGHT_CLR;
            cur_color_single = ac.RIGHT_CLR_SINGLE;
        end
        
        all_volts = [];
        
        for d_idx = 1:length(t_vel_all{e_type})
            %all_volts = vertcat(all_volts, volts_all{e_type}{d_idx}{tt}(:,:));
            tb = find((t > 2.5) & (t < 3.0) );
            baseline = squeeze(mean(mean( volts_all{e_type}{d_idx}{tt}(:, tb) )));

            all_volts(d_idx,:) = mean(volts_all{e_type}{d_idx}{tt}(:,:)) - baseline;
            
            plot(t, squeeze(all_volts(d_idx,:)), 'color', cur_color_single );
        end
        
         avg_vols = squeeze(mean(all_volts));
%         baseline = mean(avg_vols);
%         avg_vols_baseline_corr = avg_vols - repmat(baseline, [1 length(avg_vols)]);
        
        std_volt_trace = squeeze(std( all_volts,1 ))/ sqrt(size(all_volts,1));                
               
%         fh = fill([ t fliplr(t) ], [(avg_vols+std_volt_trace) fliplr((avg_vols-std_volt_trace)) ], cur_color_single );
%         set(fh, 'EdgeColor', 'None');
        
        lineStyle = '-';
        if(e_type == 2)
            lineStyle = '--';
        end
        
        pp(e_type, tt) = plot(t, avg_vols, 'color', cur_color_avg, 'LineWidth', 2, 'LineStyle', lineStyle );
    end

    xlim([0 t(end)]);
    ylim([-1.5 5])
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
xlim([1.5 5.0]);

saveas(f, [analysis_dir '/thermotaxis_behavior_and_ephys_by_fly.fig']);
saveas(f, [analysis_dir '/thermotaxis_behavior_and_ephys_by_fly.png']);

end

