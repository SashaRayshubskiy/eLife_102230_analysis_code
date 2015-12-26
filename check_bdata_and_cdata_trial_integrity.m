function [] = check_bdata_and_cdata_trial_integrity( btrial_metadata, ctrial_metadata )

trial_type_cnt = size(btrial_metadata,2);

for i=1:trial_type_cnt
    
    % First check that the trial metadata doesn't have duplicates.
    btrial_dups = get_duplicates(btrial_metadata{i});
    if(~isempty(btrial_dups))
        disp(['ERROR: Detected duplicates in behavioral data, trial type: ' num2str(i) ' trials: ' mat2str(btrial_dups) ]);
        return;
    end
    
    ctrial_dups = get_duplicates(ctrial_metadata{i});
    if(~isempty(ctrial_dups))
        disp(['ERROR: Detected duplicates in behavioral data, trial type: ' num2str(i) ' trials: ' mat2str(ctrial_dups) ]);
        return;
    end
    
    % Check that the behavioral trials are equal to the calcium trials
    if(~isequal(ctrial_metadata{i}, btrial_metadata{i}))
        disp(['ERROR: calcium trial metadata are not equal to behavioral trial metadata for trial type: ' num2str(i)]);
        return;
    end    
end

end
