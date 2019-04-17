function remove_stim_ids(analysis_path, cur_sid, stim_ids_to_remove )

old_stim_filepath = [analysis_path '/stim_window_data_sid_' num2str(cur_sid) '.mat'];
new_stim_filepath = [analysis_path '/stim_window_data_sid_' num2str(cur_sid) '_new.mat'];
    
cur_data = load( old_stim_filepath );   
  
old_bump_in_window  = cur_data.bump_in_window;
old_yaw_in_window   = cur_data.yaw_in_window;
old_fwd_in_window   = cur_data.fwd_in_window;
old_ephys_in_window = cur_data.ephys_in_window;

VPS             = cur_data.VPS;
t_bump_w        = cur_data.t_bump_w;
t_yaw_w         = cur_data.t_yaw_w;
t_ephys_w       = cur_data.t_ephys_w;

bump_in_window  = [];
yaw_in_window   = [];
fwd_in_window   = [];
ephys_in_window = [];


for i = 1:size(old_bump_in_window,1)
    if( length(find(i==stim_ids_to_remove)) == 0 )
        bump_in_window(end+1,:,:)  = old_bump_in_window(i,:,:);
        yaw_in_window(end+1,:)   = old_yaw_in_window(i,:);
        fwd_in_window(end+1,:)   = old_fwd_in_window(i,:);
        ephys_in_window(end+1,:) = old_ephys_in_window(i,:);        
    end
end

save( new_stim_filepath, ...
     'fwd_in_window', 'yaw_in_window', 'ephys_in_window', 'bump_in_window', ...
     'VPS', 't_bump_w', 't_yaw_w', 't_ephys_w' );
end

