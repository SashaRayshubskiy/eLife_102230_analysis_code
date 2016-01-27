function [ is_excluded ] = is_excluded_trial( trial_ord, trial_exclusion_list, internal_to_external_trial_id_map )

    external_id = internal_to_external_trial_id_map(trial_ord,2);
    
    is_excluded = 0;
    for i=1:size(trial_exclusion_list,2)
        if(external_id == trial_exclusion_list(i))
            is_excluded = 1;
            break;
        end
    end
end

