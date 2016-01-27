function [ external_to_internal_trial_id_map ] = get_external_trial_id_to_internal_ordinal_map( btrial_meta )

external_to_internal_trial_id_map = { containers.Map(), containers.Map(), containers.Map() };

for trial_type = 1:size(btrial_meta,2)    
    keys = btrial_meta{trial_type}(:,2);
    values = [1:size(btrial_meta{trial_type},1)];
        
    external_to_internal_trial_id_map{ trial_type } = containers.Map(keys, values);
end

end

