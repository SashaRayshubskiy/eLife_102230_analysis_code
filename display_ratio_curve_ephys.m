function display_ratio_curve_ephys(data_to_analyze, sids, t_volts, volts_all, analysis_path)

first_stim_t = 0.5;
last_stim_t =  0.6;

l_r_diff  = {};
l_r_curve = {};

for d = 1:length(t_volts)
    l_r_diff{d}  = [];
    l_r_curve{d} = [];
    
    for rat_index = 2:10
        [l,r] = get_ratio_for_index(rat_index);
        
        cur_volts = volts_all{d}{rat_index};
        
        t = squeeze(t_volts{d}{rat_index}(1,:));
        
        [~, locs] = find( t <first_stim_t);
        baseline_volts = mean(cur_volts(:,locs),2);
        
        volts_baseline_corr = cur_volts - repmat(baseline_volts, [1 size(cur_volts,2)]);
        
        avg_vols = squeeze(mean(volts_baseline_corr));
        
        [~, locs] = find( (t > first_stim_t ) & (t<last_stim_t));
        volt_during_stim = squeeze(mean(avg_vols(locs)));
        
        l_r_d = l-r;
        %l_r_diff(rat_index) = l_r_d;
        l_r_curve{d}(rat_index-1,:) = [ l_r_d, volt_during_stim];
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

ylabel('Normalized voltage diff during stim (au)');
xlabel('L-R difference');
legend('show', 'Location', 'northwest');


saveas(f1, [analysis_path '/left_right_diff_curve_' num2str(sids) '.fig']);
saveas(f1, [analysis_path '/left_right_diff_curve_' num2str(sids) '.png']);

end

