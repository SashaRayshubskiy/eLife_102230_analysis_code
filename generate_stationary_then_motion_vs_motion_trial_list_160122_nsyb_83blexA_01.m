function [ condition_trials, condition_trials_str, condition_str ] = generate_stationary_then_motion_vs_motion_trial_list_for_160122_nsyb_83blexA_01( bdata_vel_time, bdata_vel, external_trial_id_to_internal_ordinal_map )
% Data below is for 160122_nsyb_83blexA_01

ac = get_analysis_constants;
condition_trials_str = {'stationary_then_walking_after_stim', 'walking_full_time' };
condition_str = 'stationary_vs_walking_turning';

ext_condition_trials = cell(3,2);

STW = 1;
FW = 2;

ext_condition_trials{ ac.BOTH, STW } = [21,43,68,71,116,128,177,182,183,199];
ext_condition_trials{ ac.BOTH, FW } = [22,27,32,51,58,81,88,90,114,109];

ext_condition_trials{ ac.LEFT, STW } = [14,25,31,65,74,85,119,122,193,217];
ext_condition_trials{ ac.LEFT, FW } = [13,15,56,72,86,91,201,212,280,320];

ext_condition_trials{ ac.RIGHT, STW } = [30,53,57,59,63,92,112,129,173,180,181];
ext_condition_trials{ ac.RIGHT, FW } = [33,41,64,66,89,107,124,130,203,211];

condition_trials = cell(3,2);

for trial_type = 1:3
    for cond_ord = 1:2
        for trial = 1:length(ext_condition_trials{trial_type, cond_ord})
            ext_trial = ext_condition_trials{trial_type, cond_ord}(trial);
            condition_trials{trial_type, cond_ord}(end+1) = external_trial_id_to_internal_ordinal_map{trial_type}( ext_trial );
        end
    end
end
end

