function [ b_rawdata, b_metadata ] = load_behavioral_data( datapath, sids )

% bdata = { analysis trial id, { vel_fwd, vel_side, vel_yaw, vel_time}, {frame_clock, stim, meta_time} } 
% analysis trial id maps 1-to-1 behavior and imaging data for this analysis
% session

global slash;

trial_type_cnt = 3;

for i=1:trial_type_cnt
    bdata_not_sorted(i,1) = {};
    trial_type_idx_counter(i) = 1;
end

i=1;

for sid = sids
    files = dir([datapath slash '*_sid_' num2str(sid) '_*.mat'])

    for f = 1:length(files)
       
        filename = files(f).name;
        filename_split = strsplit(filename, '_');
        trial_type = filename_split(2);
        sid = filename_split(6);
        tid = filename_split(8);
        
        trial_type_idx = -1;
        if(strcmp(trial_type, 'BothOdor') == 1)
            trial_type_idx = 1;
        elseif(strcmp(trial_type, 'LeftOdor') == 1)
            trial_type_idx = 2;        
        elseif(strcmp(trial_type, 'RightOdor') == 1)
            trial_type_idx = 3;
        else
            disp(['ERROR: Trial type: ' trial_type ' is not recognized']);
        end
        
        raw_data = load([ datapath slash filename ]);
        
        bdata_not_sorted{ trial_type_idx }( end+1 ) = { sid, tid, raw_data };
    end  
end

% Sort 
for i=1:trial_type_cnt
    bdata_sorted{ i } = sortrows( bdata_not_sorted(i), [1 2] );
end

for i=1:trial_type_cnt
    cur_trial_cnt = size(bdata_sorted{i},1);
    cur_raw_data_size = size(bdata_sorted{i}{3},1);
    
    cur_bdata = zeros([cur_trial_cnt, cur_raw_data_size], 'double');

    for j=1:cur_trial_cnt
        cur_bdata(j,:,:) = bdata_sorted{i}{3};
        b_metadata{i}(j,:) = [ bdata_sorted{i}{1}, bdata_sorted{i}{2} ];
    end
    bdata_raw{j} = cur_bdata;
end

end