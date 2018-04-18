function [ locs ] = find_thresholds_in_PSTH( t_vel_all, psth, FR_THRESHOLD )

% 
EXCLUSIVE_WINDOW = 0.1; % seconds
dt = t_vel_all(2) - t_vel_all(1);

prev = psth(1);
time_since_last_crossing = 0;
locs = [];
for i = 2:length(psth)
    
    cur = psth(i);
    
    if( (cur > prev) && (cur >= FR_THRESHOLD) && (prev < FR_THRESHOLD) && (time_since_last_crossing > EXCLUSIVE_WINDOW))
        locs(end+1) = i;
        time_since_last_crossing = 0;
    else
        time_since_last_crossing = time_since_last_crossing + dt;
    end

    prev = cur;
end

end

