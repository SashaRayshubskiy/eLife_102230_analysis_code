function [ct_list] = generate_initial_turn_trials( bdata_vel_time, bdata_vel, analysis_path )

ct_list = cell(2,1);
ct_list{1} = [];
ct_list{2} = [];


turn_metadata = generate_turning_metadata( 0, bdata_vel_time, bdata_vel, analysis_path );

[ condition_trials, condition_trials_str, condition_str ] = generate_expected_vs_ignore_trial_list_v3( 0, bdata_vel_time, bdata_vel, turn_metadata, analysis_path );

avg_cond_btrace_trace_filepath = [ analysis_path '/' condition_str '_asid_' num2str( 0 ) '_sid_' num2str( 0 ) ];

with_single_trials = 0;
display_two_condition_trials_avg( condition_trials, condition_trials_str, bdata_vel_time, bdata_vel, avg_cond_btrace_trace_filepath, with_single_trials );

% for t_type = 1:2
%     ct_list{t_type} = condition_trials{t_type, 1};
% end

ct_list = condition_trials;

end

