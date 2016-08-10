function [ output_args ] = display_avg_velocity_left_right_only( sid, bdata_raw, bdata_vel, bdata_vel_time, analysis_path, with_single_trials )

ac = get_analysis_constants;
settings = sensor_settings;

one_trial_bdata = squeeze(bdata_raw{ ac.LEFT }(1,:,:));
left_odor_stim  = squeeze(one_trial_bdata(:,7));
first_stim = find((left_odor_stim > 0.5), 1, 'first') ./ settings.sampRate;
last_stim = find((left_odor_stim > 0.5), 1, 'last') ./ settings.sampRate;

first_stim = settings.pre_stim;
last_stim = settings.pre_stim+settings.stim;

f = figure('units','normalized','outerposition',[0 0 1 1]);

start_trial = size(bdata_vel{ac.LEFT},1)/2;
end_left_trial = size(bdata_vel{ac.LEFT},1);
end_right_trial = size(bdata_vel{ac.RIGHT},1);

subplot(3,1,1);
hold on;

total_trial_cnt = size(bdata_vel{ac.LEFT}(start_trial:end_left_trial,:,:),1) + size(bdata_vel{ac.RIGHT}(start_trial:end_right_trial,:,:),1);
mean_left_vel_fwd = squeeze(mean(bdata_vel{ac.LEFT}(start_trial:end_left_trial,ac.VEL_FWD,:)));
mean_right_vel_fwd = squeeze(mean(bdata_vel{ac.RIGHT}(start_trial:end_right_trial,ac.VEL_FWD,:)));

SHOW_SEM = 1;
if( SHOW_SEM )
    sem_left_vel_fwd = squeeze(std(bdata_vel{ac.LEFT}(start_trial:end_left_trial,ac.VEL_FWD,:),1)) ./ sqrt(size(bdata_vel{ ac.LEFT }(start_trial:end_left_trial,ac.VEL_FWD,:),1));
    sem_right_vel_fwd = squeeze(std(bdata_vel{ac.RIGHT}(start_trial:end_right_trial,ac.VEL_FWD,:),1)) ./ sqrt(size(bdata_vel{ ac.RIGHT }(start_trial:end_right_trial,ac.VEL_FWD,:),1));
else
    sem_left_vel_fwd = squeeze(std(bdata_vel{ac.LEFT}(start_trial:end_left_trial,ac.VEL_FWD,:),1));
    sem_right_vel_fwd = squeeze(std(bdata_vel{ac.RIGHT}(start_trial:end_right_trial,ac.VEL_FWD,:),1));
end

fh = fill( [bdata_vel_time, fliplr(bdata_vel_time)], ... 
        [(mean_left_vel_fwd+sem_left_vel_fwd)' fliplr((mean_left_vel_fwd-sem_left_vel_fwd)')], ...
        rgb('Salmon'));
set(fh, 'EdgeColor', 'None');

fh = fill( [bdata_vel_time, fliplr(bdata_vel_time)], ... 
        [(mean_right_vel_fwd+sem_right_vel_fwd)' fliplr((mean_right_vel_fwd-sem_right_vel_fwd)')], ...
        rgb('PaleGreen'));
set(fh, 'EdgeColor', 'None');


p1 = plot(bdata_vel_time, mean_left_vel_fwd, 'color', ac.LEFT_CLR);
p2 = plot(bdata_vel_time, mean_right_vel_fwd, 'color', ac.RIGHT_CLR);
ylim([-0.01 0.05]);

legend([p1,p2], 'Left Odor w/ SEM', 'Right Odor w/ SEM');

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ first_stim first_stim last_stim last_stim ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

xlim([0, bdata_vel_time(end)]);
xlabel('Time (s)');
ylabel('Velocity (au/s)');
title(['Forward velocity. Trial count: ' num2str(total_trial_cnt)]);

subplot(3,1,2);
hold on;

mean_left_vel_side = squeeze(mean(bdata_vel{ac.LEFT}(start_trial:end_left_trial,ac.VEL_SIDE,:)));
mean_right_vel_side = squeeze(mean(bdata_vel{ac.RIGHT}(start_trial:end_right_trial,ac.VEL_SIDE,:)));

sem_left_vel_side = squeeze(std(bdata_vel{ac.LEFT}(start_trial:end_left_trial,ac.VEL_SIDE,:),1)) ./ sqrt(size(bdata_vel{ ac.LEFT }(1:end_left_trial,ac.VEL_SIDE,:),1));
sem_right_vel_side = squeeze(std(bdata_vel{ac.RIGHT}(start_trial:end_right_trial,ac.VEL_SIDE,:),1)) ./ sqrt(size(bdata_vel{ ac.RIGHT }(1:end_right_trial,ac.VEL_SIDE,:),1));

fh = fill( [bdata_vel_time, fliplr(bdata_vel_time)], ... 
        [(mean_left_vel_side+sem_left_vel_side)' fliplr((mean_left_vel_side-sem_left_vel_side)')], ...
        rgb('Salmon'));
set(fh, 'EdgeColor', 'None');

fh = fill( [bdata_vel_time, fliplr(bdata_vel_time)], ... 
        [(mean_right_vel_side+sem_right_vel_side)' fliplr((mean_right_vel_side-sem_right_vel_side)')], ...
        rgb('PaleGreen'));
set(fh, 'EdgeColor', 'None');

p1 = plot(bdata_vel_time, mean_left_vel_side, 'color', ac.LEFT_CLR);
p2 = plot(bdata_vel_time, mean_right_vel_side, 'color', ac.RIGHT_CLR);
ylim([-0.02 0.02]);

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ first_stim first_stim last_stim last_stim ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

xlim([0, bdata_vel_time(end)]);
xlabel('Time (s)');
ylabel('Velocity (au/s)');
title('Side velocity');

subplot(3,1,3);
hold on;

mean_left_vel_yaw = squeeze(mean(bdata_vel{ac.LEFT}(start_trial:end_left_trial,ac.VEL_YAW,:)));
mean_right_vel_yaw = squeeze(mean(bdata_vel{ac.RIGHT}(start_trial:end_right_trial,ac.VEL_YAW,:)));

sem_left_vel_yaw = squeeze(std(bdata_vel{ac.LEFT}(start_trial:end_left_trial,ac.VEL_YAW,:),1)) ./ sqrt(size(bdata_vel{ ac.LEFT }(start_trial:end_left_trial,ac.VEL_YAW,:),1));
sem_right_vel_yaw = squeeze(std(bdata_vel{ac.RIGHT}(start_trial:end_right_trial,ac.VEL_YAW,:),1)) ./ sqrt(size(bdata_vel{ ac.RIGHT }(start_trial:end_right_trial,ac.VEL_YAW,:),1));

fh = fill( [bdata_vel_time, fliplr(bdata_vel_time)], ... 
        [(mean_left_vel_yaw+sem_left_vel_yaw)' fliplr((mean_left_vel_yaw-sem_left_vel_yaw)')], ...
        rgb('Salmon'));
set(fh, 'EdgeColor', 'None');

fh = fill( [bdata_vel_time, fliplr(bdata_vel_time)], ... 
        [(mean_right_vel_yaw+sem_right_vel_yaw)' fliplr((mean_right_vel_yaw-sem_right_vel_yaw)')], ...
        rgb('PaleGreen'));
set(fh, 'EdgeColor', 'None');

p1 = plot(bdata_vel_time, mean_left_vel_yaw, 'color', ac.LEFT_CLR);
p2 = plot(bdata_vel_time, mean_right_vel_yaw, 'color', ac.RIGHT_CLR);
ylim([-0.1 0.1]);

xlim([0, bdata_vel_time(end)]);
xlabel('Time (s)');
ylabel('Velocity (au/s)');
title('Yaw velocity');

yy = ylim;
y_min = yy(1); y_max = yy(2);
hh = fill([ first_stim first_stim last_stim last_stim ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

saveas(f, [analysis_path '/avg_velocity_sid_' num2str( sid ) '.fig']);
saveas(f, [analysis_path '/avg_velocity_sid_' num2str( sid ) '.png']);

end

