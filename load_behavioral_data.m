function [ bdata ] = load_behavioral_data( datapath, sids )

% bdata = { analysis trial id, { vel_fwd, vel_side, vel_yaw, vel_time}, {frame_clock, stim, meta_time} } 
% analysis trial id maps 1-to-1 behavior and imaging data for this analysis
% session

global slash;

for sid = sids
    files = dir([datapath slash '*_sid_' num2str(sid) '_*.mat'])

    for f = 1:length(files)
        raw_data = load([datapath slash files(f).name]);
        
            
    end  
end

