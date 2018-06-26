function [ cdata_raw, cdata_meta, trial_metadata ] = load_imaging_data( sids, datapath, trial_type_cnt, dx, dy, dt )

% Inputs:
% sids - array or a single SID of the trials to load
% Path to data directory for desired experiment 

% Outputs:
% cdata_raw = { trial_type } [ analysis trial id, Lines, Pixels, Channels, Planes, Volumes ]
% 
% cdata_meta.frame_rate = frames_per_second
% cdata_meta.volume_rate = volumes per second
% 
% trial_metadata = { trial_type } [ analysis trial id, sid, tid ]
%                  sid is an identifier of the session, tid is the trial id
%                  within the session. These should map 1-to-1 on
%                  behavioral trial data.

% analysis trial id maps 1-to-1 behavior and imaging data for this analysis
% session

global slash;

for i=1:trial_type_cnt
    cdata_not_sorted{i} = [];
    cdata_not_sorted_helper{i} = {};
    cdata_idx_per_type(i) = 1;
end

i=1;

for sid = sids
    files = dir([datapath slash '*_sid_' num2str(sid) '_*.tif']);

    for f = 1:length(files)
       
        filename = files(f).name;
        filename_split = strsplit(filename, '_');
        trial_type = filename_split(2);
        
        sid = str2num(char(filename_split(6)));
        
        tid_str = filename_split(8);
        tid = str2num(char(tid_str));
        
        trial_type_idx = -1;
        if( (strcmp(trial_type, 'LeftOdor') == 1) | (strcmp(trial_type, 'OdorLeftWind') == 1) | (strcmp(trial_type, 'OpenLoopLeft') == 1) | (strcmp(trial_type, 'PicoPump') == 1) )
            trial_type_idx = 1;        
        elseif((strcmp(trial_type, 'RightOdor') == 1) | (strcmp(trial_type, 'OdorRightWind') == 1) | (strcmp(trial_type, 'OpenLoopRight') == 1) )
            trial_type_idx = 2;
        elseif((strcmp(trial_type, 'BothOdor') == 1) | (strcmp(trial_type, 'NaturalOdor') == 1) | (strcmp(trial_type, 'OdorNoWind') == 1))
            trial_type_idx = 3;
        else
            disp(['ERROR: Trial type: ' trial_type ' is not recognized']);
            return;
        end
        
        load_path = [ datapath slash filename ];
        raw_data = open_tif_fast( load_path );
 
        x_size = size(raw_data, 1);
        y_size = size(raw_data, 2);
        t_size = size(raw_data, 5);

        raw_d_down = squeeze(mean(mean(mean(reshape(raw_data, [dx, x_size/dx, dy, y_size/dy, size( raw_data, 3), size(raw_data,4), dt, t_size/dt ]), 7),3),1));
        
        disp(['Loaded file: ' load_path]);
        
        rid = cdata_idx_per_type(trial_type_idx);

        cdata_not_sorted{ trial_type_idx }( end+1, : ) = [ sid, tid, rid ];
        cdata_not_sorted_helper{ trial_type_idx }( end+1, : ) = { rid, raw_d_down };
        
        cdata_idx_per_type(trial_type_idx) = cdata_idx_per_type(trial_type_idx) + 1;
    end  
end

% Sort 
for i=1:trial_type_cnt
    cdata_sorted{ i } = sortrows( cdata_not_sorted{i}, [1 2] );
end

for i=1:trial_type_cnt    
    cur_trial_cnt = size(cdata_sorted{i},1);
    
    % [Lines, Pixels, Channels, Planes, Volumes]
    
    [clines, cpixels, cchannels, cplanes, cvolumes] = size(cdata_not_sorted_helper{i}{1,2});
    
    cur_bdata = zeros([cur_trial_cnt, clines, cpixels, cchannels, cplanes, cvolumes], 'double');

    for j=1:cur_trial_cnt
        rid = cdata_sorted{i}(j,3);
        rid_ptr = cdata_not_sorted_helper{i}{rid,1};
        
        if(rid ~= rid_ptr)
            disp(['ERROR: rids did not match: [ ' num2str(rid) ' , ' num2str(rid_ptr) ' ]']);
            return;
        end
               
        cur_bdata( j, :, :, :, :, : ) = cdata_not_sorted_helper{i}{rid,2};
        trial_metadata{ i }( j, : ) = cdata_sorted{i}(j,1:2);
    end
    
    cdata_raw{i} = cur_bdata;
end

% Set some useful metadata about the imaging session
if(length(sids) == 1)
    files = dir([datapath slash '*_sid_' num2str(sids) '_*.tif']);
else
    files = dir([datapath slash '*_sid_' num2str(sids(1)) '_*.tif']);    
end

tifpath = [ datapath slash files(1).name ];
tifObj = Tiff(tifpath,'r');

frameString = tifObj.getTag('ImageDescription');
cdata_meta.frame_rate = si51_frame_string_get_value_for_key(frameString, 'hRoiManager.scanFrameRate');
cdata_meta.volume_rate = 1.0/si51_frame_string_get_value_for_key(frameString, 'hFastZ.period');

tifObj.close();

end

