%%

dd = load('/tmp/bump_debug.mat');

[ smoothed_bump, cur_bump_tc, cur_bump_tc_unwrapped, vect_strength_check ] = get_radial_weighted_avg_bump_pos_vect_strengh_check_v2( dd.bump_debug );

figure;

hold on;
plot( cur_bump_tc );

baseline_vals = cur_bump_tc( dd.bump_baseline_idx );
baseline_non_nan = baseline_vals(~isnan(baseline_vals));
bump_delta_tc = medfilt1( cur_bump_tc, 3, 'truncate' );

plot( bump_delta_tc, 'g' );