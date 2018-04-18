function merge_sids( datapath, sid_1, sid_2)

ball_path = [datapath '/ball/'];
tp_path = [datapath '/2p/'];

ball_search_str_1 = ['bdata_*_sid_' num2str(sid_1) '_*.mat'];
ball_search_str_2 = ['bdata_*_sid_' num2str(sid_2) '_*.mat'];

ball_files_1 = dirs(ball_path, ball_search_str_1);
ball_files_2 = dirs(ball_path, ball_search_str_2);


tp_search_str_1 = ['cdata_*_sid_' num2str(sid_1) '_*.tif'];
tp_search_str_2 = ['cdata_*_sid_' num2str(sid_2) '_*.tif'];

tp_files_1 = dirs(tp_path, tp_search_str_1);
tp_files_2 = dirs(tp_path, tp_search_str_2);



end

