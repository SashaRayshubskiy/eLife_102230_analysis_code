function [] = check_bdata_and_cdata_trial_integrity( btrial_metadata, ctrial_metadata )

% First check that the trial metadata doesn't have duplicates.
btrial_dups = get_duplicates(btrial_metadata);
if(~isempty(btrial_dups))    
    disp(['ERROR: Detected duplicates in behavioral data, trials: ' mat2str(btrial_dups) ]);
    return;
end

ctrial_dups = get_duplicates(ctrial_metadata);
if(~isempty(cbtrial_dups))    
    disp(['ERROR: Detected duplicates in calcium data, trials: ' mat2str(ctrial_dups) ]);
    return;
end

% Check that the behavioral trials are equal to the calcium trials
if(~isequal(ctrial_metadata, btrial_metadata))
    disp(['ERROR: calcium trial metadata are not equal to behavioral trial metadata']);
    return;
end

end

