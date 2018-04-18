function generate_random_shift( t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH, analysis_path_type )

f = figure;
rsqs = [];

ITER_CNT = 100;
iter_idx = 1;
while(iter_idx <= ITER_CNT)
    hold on;
    [ r_val, shift_amount_in_sec ] = generate_Vm_vs_yaw_trajectory_LAL_DNs_shift_search(-1, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH, analysis_path_type);
    
    disp( [ 'Random R^2: ' num2str( r_val ) ] );
    
    rsqs(iter_idx) = r_val;
    iter_idx = iter_idx + 1;
end

saveas(f, [ analysis_path '/fit_curves_shift_factor.fig' ]);
saveas(f, [ analysis_path '/fit_curves_shift_factor.png' ]);

f = figure;

plot( rsqs, '-*');
title(['R^2 for random shifts']);
xlabel('Iteration index (au)');
ylabel('R^2 value');
set(gca, 'FontSize', 16);

saveas(f, [ analysis_path '/ramdon_shifts_r_sqr_value.fig' ]);
saveas(f, [ analysis_path '/ramdon_shifts_r_sqr_value.png' ]);

save([ analysis_path '/ramdon_shifts_r_sqr_value.mat' ], 'rsqs');

disp(['Mean random R^2 from: (' num2str(ITER_CNT) ') :: ' num2str(mean(rsqs))]);

end

