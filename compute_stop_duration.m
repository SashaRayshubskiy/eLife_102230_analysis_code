function [stop_dur] = compute_stop_duration( fwd_vel, time_between_samples, stopping_threshold )

stopped_idx = find( fwd_vel <= stopping_threshold );

continuous_cnts = [];
cur_cont_run_idx = 1;
continuous_cnts( cur_cont_run_idx ) = 0;
% find longest continuous run of stopping
for i = [1:length(stopped_idx)-1]
    if(stopped_idx(i+1)-stopped_idx(i) == 1)
        continuous_cnts( cur_cont_run_idx ) = continuous_cnts( cur_cont_run_idx ) + 1;
    else
        cur_cont_run_idx = cur_cont_run_idx + 1;
        continuous_cnts( cur_cont_run_idx ) = 0;
    end
end

longest_run = max( continuous_cnts );

stop_dur = longest_run * time_between_samples;

end
