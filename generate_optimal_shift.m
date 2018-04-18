function generate_optimal_shift( t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH, analysis_path_type )

SHIFT_MAX = 100;
SHIFT_MIN = 1;

datafile = [ analysis_path '/shift_vs_r_sqr_value.mat' ];

shifts = [];
rsqs = [];

if(~exist(datafile,'file'))
    f = figure;
    
    for shift = [ SHIFT_MIN : SHIFT_MAX ]
        hold on;
        [ r_val, shift_amount_in_sec ] = generate_Vm_vs_yaw_trajectory_LAL_DNs_shift_search(shift, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path, CALC_PSTH, analysis_path_type);
        
        disp( [ 'shift in sec: ' num2str(shift*shift_amount_in_sec) '  R^2: ' num2str( r_val ) ] );
        
        shifts(end+1) = shift;
        rsqs(end+1) = r_val;
    end
    
    saveas(f, [ analysis_path '/fit_curves_shift_factor.fig' ]);
    saveas(f, [ analysis_path '/fit_curves_shift_factor.png' ]);
    
    save(datafile, 'shifts', 'rsqs', 'shift_amount_in_sec');
else
    dd                  = load(datafile);
    shifts              = dd.shifts;
    rsqs                = dd.rsqs;
    shift_amount_in_sec = dd.shift_amount_in_sec;
end

f = figure;
hold on;
shifts_in_sec = shifts*shift_amount_in_sec;
plot(shifts_in_sec, rsqs, '-*', 'DisplayName', 'Shift');

% Mean random R^2 from: 100:: 0.00019233
% calculated separately, generate_random_shift
random_rsq = 0.00019233;
random_shifts = repmat(random_rsq, [1 length(shifts)]);
plot(shifts_in_sec, random_shifts, '--', 'DisplayName', 'Shuffle');

title(['Shift amount in seconds: ' num2str(shift_amount_in_sec)]);
xlabel('Shift Time (s)');
ylabel('R^2 value');
set(gca, 'FontSize', 16);
legend('show');

saveas(f, [ analysis_path '/shift_vs_r_sqr_value.fig' ]);
saveas(f, [ analysis_path '/shift_vs_r_sqr_value.png' ]);

end

