function display_ratio_curve_ephys_PSTH_based(data_to_analyze, t_volts, t_vel_all, volts_all, analysis_path)

first_stim_t = 0.5;
last_stim_t =  1.0;

l_r_diff  = {};
l_r_curve = {};

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

psth_dt_samples = ephys_SR/ball_SR;
psth_dt = psth_dt_samples / (1.0*ephys_SR);

SPIKE_THRESHOLD_LAL_DN = 0.5; % For A2

for d = 1:length(t_volts)
    l_r_diff{d}  = [];
    l_r_curve{d} = [];
    
    for rat_index = 2:10
        [l,r] = get_ratio_for_index(rat_index);
        
        trial_cnt = size(volts_all{d}{rat_index},1);
        psths = [];
        cur_volts = volts_all{d}{rat_index};
        for tt = 1:trial_cnt
            cur_psth = calculate_psth( t_volts, t_vel_all{d}, cur_volts(tt,:), ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples );
            psths(tt,:) = cur_psth;
        end
            
        t     = squeeze( t_volts{d}{rat_index}(1,:) );
        t_vel = squeeze( t_vel_all{d}{rat_index}(1,:) );            
        
        [~, locs] = find( (t_vel > first_stim_t ) & (t_vel < last_stim_t));
        psth_during_stim = squeeze(mean(squeeze(mean(psths(:,locs)))));
        
        l_r_d = l-r;
        %l_r_diff(rat_index) = l_r_d;
        l_r_curve{d}(rat_index-1,:) = [ l_r_d, psth_during_stim];
    end
end

f1 = figure;
hold on;

for d = 1:length(t_volts)
    total_trials = data_to_analyze{d}{3};
    
    %norm_fact = (squeeze(l_r_curve{d}(1,2)));    % {1,0} and {0,1}
    norm_fact = (squeeze(l_r_curve{d}(1,2)) + squeeze(l_r_curve{d}(end,2)))/2.0;    % {1,0} and {0,1}
    l_r_curve{d}(:,2) = l_r_curve{d}(:,2) ./ norm_fact;
    
    l_r_curve_s(d,:,:) = sortrows(l_r_curve{d},1);
    pl = plot( squeeze(l_r_curve_s(d,:,1)), squeeze(l_r_curve_s(d,:,2)), 'DisplayName', ['Fly: ' num2str(d) ' trials per ratio: ' num2str(total_trials/10)]);
    pl.Marker = '*';
end

plot(squeeze(l_r_curve_s(1,:,1)), squeeze(mean(l_r_curve_s(:,:,2))), 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Avg');

ylabel('Normalized PSTH during stim (au)');
xlabel('L-R difference');
legend('show', 'Location', 'northwest');

set(gca,'FontSize', 16);

saveas(f1, [analysis_path '/left_right_diff_curve_psth_based.fig']);
saveas(f1, [analysis_path '/left_right_diff_curve_psth_based.png']);

end

