function [ duplicated_trials ] = get_duplicates( trial_metadata )

% Inputs:
% [analysis_trial_idx, sid, tid]

tm_unique = unique(trial_metadata);
countofU = hist(trial_metadata, tm_unique);
ii = (countofU > 1);

duplicated_trials = tm_unique(ii(:,2));

end

