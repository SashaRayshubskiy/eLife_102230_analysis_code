function [ bdata_raw, bdata_time, trial_metadata ] = load_natural_odor_daq_data( sids, datapath )

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

bdata_not_sorted = [];
bdata_not_sorted_helper = {};
bdata_idx_per_type = 1;

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
        if((strcmp(trial_type, 'NaturalOdor') == 1))
        else
            disp(['ERROR: Trial type: ' trial_type ' is not recognized']);
            return;
        end
        
        load_path = [ datapath slash filename ];
        raw_data = load( load_path );
        disp(['Loaded file: ' load_path]);
        
        rid = bdata_idx_per_type;

        bdata_not_sorted( end+1, : ) = [ sid, tid, rid ];
        bdata_not_sorted_helper( end+1, : ) = { rid, raw_data };
        
        bdata_idx_per_type = bdata_idx_per_type + 1;
    end  
end

% Sort 
bdata_sorted = sortrows( bdata_not_sorted, [1 2] );

cur_trial_cnt = size(bdata_sorted,1);
[num_samples, num_channels] = size(bdata_not_sorted_helper{1,2}.trial_bdata);

cur_bdata = zeros([cur_trial_cnt, num_samples, num_channels], 'double');

for j=1:cur_trial_cnt
    rid = bdata_sorted(j,3);
    rid_ptr = bdata_not_sorted_helper{rid,1};
    
    if(rid ~= rid_ptr)
        disp(['ERROR: rids did not match: [ ' num2str(rid) ' , ' num2str(rid_ptr) ' ]']);
        return;
    end
    
    cur_bdata( j, :, : ) = bdata_not_sorted_helper{rid,2}.trial_bdata;
    trial_metadata( j, : ) = bdata_sorted(j,1:2);
end

bdata_raw = cur_bdata;

bdata_time = bdata_not_sorted_helper{1,2}.trial_time;

end