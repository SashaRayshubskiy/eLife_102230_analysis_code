function display_ratio_curve_behavior_only(data_to_analyze, t_vel_all, yaw_vel_all, fwd_vel_all, analysis_path)

before_light_on_t = 0.5;
stim_window_start_t = 0.60;
stim_window_end_t =  1.0;

l_r_diff  = {};
l_r_curve = {};

for d = 1:length(t_vel_all)
    l_r_diff{d}  = [];
    l_r_curve{d} = [];
    
    for rat_index = 2:10
        [l,r] = get_ratio_for_index(rat_index);
        
        cur_yaw = yaw_vel_all{d}{rat_index};
        cur_fwd = fwd_vel_all{d}{rat_index};
        
        t = squeeze(t_vel_all{d}{rat_index}(1,:));
        
        [~, locs] = find( t <stim_window_start_t);
        baseline_yaw = mean(cur_yaw(:,locs),2);                
        yaw_baseline_corr = cur_yaw - repmat(baseline_yaw, [1 size(cur_yaw,2)]);
        
        avg_yaw = squeeze(mean(yaw_baseline_corr));      
        
        [~, locs] = find( (t > stim_window_start_t ) & (t<stim_window_end_t));
        yaw_during_stim = squeeze(mean(avg_yaw(locs)));
        
        [~, locs] = find( (t < before_light_on_t ) );
        avg_fwd = mean(cur_fwd);
        avg_fwd_before_stim = squeeze(mean(avg_fwd(locs)));
        
        
        l_r_d = l-r;
        %l_r_diff(rat_index) = l_r_d;
        l_r_curve{d}(rat_index-1,:) = [ l_r_d, yaw_during_stim, avg_fwd_before_stim];
        
        disp(['l_r_d: ' num2str(l_r_d) '  yaw_during_stim: ' num2str(yaw_during_stim) ]);
    end
end

f1 = figure;

subplot(2,1,1);
hold on;
for d = 1:length(t_vel_all)
    total_trials = data_to_analyze{d}{3};
    l_r_curve_s(d,:,:) = sortrows(l_r_curve{d},1);
    pl = plot( squeeze(l_r_curve_s(d,:,1)), squeeze(l_r_curve_s(d,:,3)), 'DisplayName', ['Fly: ' num2str(d) ' trials per ratio: ' num2str(total_trials/10)]);
    pl.Marker = '*';
end

plot(squeeze(l_r_curve_s(1,:,1)), squeeze(mean(l_r_curve_s(:,:,3))), 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Avg');
ylabel('Fwd (mm/s)');

subplot(2,1,2);

hold on;
for d = 1:length(t_vel_all)
    total_trials = data_to_analyze{d}{3};
    
    % Normalize yaw
    norm_fact = (squeeze(l_r_curve_s(d,1,2)) + squeeze(l_r_curve_s(d,end,2)))/2;    
    l_r_curve{d}(:,2) = l_r_curve{d}(:,2) ./ norm_fact;
        
    l_r_curve_s(d,:,:) = sortrows(l_r_curve{d},1);
    pl = plot( squeeze(l_r_curve_s(d,:,1)), squeeze(l_r_curve_s(d,:,2)), 'DisplayName', ['Fly: ' num2str(d) ' trials per ratio: ' num2str(total_trials/10)]);
    pl.Marker = '*';
end

plot(squeeze(l_r_curve_s(1,:,1)), squeeze(mean(l_r_curve_s(:,:,2))), 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Avg');

ylabel('Yaw vel normalized (au)');
xlabel('L-R difference');

legend('show');

saveas(f1, [analysis_path '/behavior_only_left_right_diff_curve.fig']);
saveas(f1, [analysis_path '/behavior_only_left_right_diff_curve.png']);

end

