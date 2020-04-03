function [ bdata_raw, bdata_time, trial_metadata ] = load_behavioral_data( sids, datapath, trial_type_cnt)

% Inputs:
% sids - array or a single SID of the trials to load
% Path to data directory for desired experiment 

% Outputs:
% bdata_raw = { trial_type } [ analysis trial id, number of samples, daq channel }
% bdata_time = the daq time base, time stamp for each sample point 
% trial_metadata = { trial_type } [ analysis trial id, sid, tid ]
%                  sid is an identifier of the session, tid is the trial id
%                  within the session. These should map 1-to-1 on imaging
%                  trial data.

% analysis trial id maps 1-to-1 behavior and imaging data for this analysis
% session

global slash;

for i=1:trial_type_cnt
    bdata_not_sorted{i} = [];
    bdata_not_sorted_helper{i} = {};
    bdata_idx_per_type(i) = 1;
end

i=1;

for sid = sids
    files = dir([datapath slash '*_sid_' num2str(sid) '_*.mat']);

    for f = 1:length(files)
       
        filename = files(f).name;
        filename_split = strsplit(filename, '_');
        trial_type = filename_split(2);
        
        sid = str2num(char(filename_split(6)));
        
        tid_str = filename_split(8);
        tid_split = strsplit(char(tid_str), '.');
        tid = str2num(char(tid_split(1)));
        
        trial_type_idx = -1;

        if((strcmp(trial_type, 'LeftOdor') == 1) | (strcmp(trial_type, '2pStim') == 1) | (strcmp(trial_type, 'OdorLeftWind') == 1) | (strcmp(trial_type, 'PicoPump') == 1) )
            trial_type_idx = 1;        
        elseif( (strcmp(trial_type, 'RightOdor') == 1) | (strcmp(trial_type, 'OdorRightWind') == 1) )
            trial_type_idx = 2;
        elseif((strcmp(trial_type, 'BothOdor') == 1) | (strcmp(trial_type, 'NaturalOdor') == 1) | (strcmp(trial_type, 'OdorNoWind') == 1))
            trial_type_idx = 3;
        else
            disp(['ERROR: Trial type: ' trial_type ' is not recognized']);
            return;
        end
        
        load_path = [ datapath slash filename ];
        raw_data = load( load_path );
        % disp(['Loaded file: ' load_path]);
        
        rid = bdata_idx_per_type(trial_type_idx);

        bdata_not_sorted{ trial_type_idx }( end+1, : ) = [ sid, tid, rid ];
        bdata_not_sorted_helper{ trial_type_idx }( end+1, : ) = { rid, raw_data };
        
        bdata_idx_per_type(trial_type_idx) = bdata_idx_per_type(trial_type_idx) + 1;
    end  
end

% Sort 
for i=1:trial_type_cnt
    bdata_sorted{ i } = sortrows( bdata_not_sorted{i}, [1 2] );
end

for i=1:trial_type_cnt
    
    cur_trial_cnt = size(bdata_sorted{i},1);
    [num_samples, num_channels] = size(bdata_not_sorted_helper{i}{1,2}.trial_bdata);
    
    cur_bdata = zeros([cur_trial_cnt, num_samples, num_channels], 'double');

    for j=1:cur_trial_cnt
        rid = bdata_sorted{i}(j,3);
        rid_ptr = bdata_not_sorted_helper{i}{rid,1};
        
        if(rid ~= rid_ptr)
            disp(['ERROR: rids did not match: [ ' num2str(rid) ' , ' num2str(rid_ptr) ' ]']);
            return;
        end
               
        cur_bdata( j, :, : ) = bdata_not_sorted_helper{i}{rid,2}.trial_bdata;
        
        % [sid, tid]
        trial_metadata{ i }( j, : ) = bdata_sorted{i}(j,1:2);
    end
    
    bdata_raw{i} = cur_bdata;
end

bdata_time = bdata_not_sorted_helper{1}{1,2}.trial_time;

end