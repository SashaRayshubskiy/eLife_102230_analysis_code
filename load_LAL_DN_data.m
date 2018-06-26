function [ t_all, t_vel_all, yaw_all, fwd_all, ephys_all_A, ephys_all_B ] = load_LAL_DN_data( working_dir, directories_to_analyze, ephys_SR, ball_SR, get_B )

dataset_cnt = length(directories_to_analyze);

t_all = cell(1,dataset_cnt);
t_vel_all = cell(1,dataset_cnt);
yaw_all = cell(1,dataset_cnt);
fwd_all = cell(1,dataset_cnt);
ephys_all_A = cell(1,dataset_cnt);
ephys_all_B = cell(1,dataset_cnt);

for d = 1:dataset_cnt
    
    datapath = [ working_dir directories_to_analyze{ d }{ 1 } ];
    
    sids = [ directories_to_analyze{ d }{ 2 }];
    
    disp(['Loading directory: ' datapath '   sid:' num2str(sids) ]);
    
    TIME_CUTOFF_SECONDS = directories_to_analyze{ d }{ 3 }; % Analysis for: 161208_lexAOpCsChrimson_gfp_83blexA_ss730_05, sid=0
    
    dd = {};
    for s = 1:length(sids)
                
        dd_tmp = dir([datapath '/*_sid_' num2str(sids(s)) '_*.mat']);
        
        % Make sure that the data is in sorted order
        num_trials = length(dd_tmp);
        sort_vector1 = zeros(num_trials, 3);
        
        for i=1:num_trials
            cur_name = dd_tmp(i).name;
            
            cur_name_token = strsplit(cur_name, '_');
            sid = str2double(cur_name_token{6});
            tid_str = strsplit(cur_name_token{8}, '.');
            
            tid = str2double(tid_str{1});
            
            sort_vector1(i, 1) = sid;
            sort_vector1(i, 2) = tid;
            sort_vector1(i, 3) = i;
        end
        
        sort_result_1 = sortrows(sort_vector1);
        
        dd_sorted = cell(1,num_trials);
        
        for i=1:num_trials
            cur_idx = sort_result_1(i,3);
            dd_sorted{i} = dd_tmp(cur_idx).name;
        end
        
        dd = horzcat(dd, dd_sorted);
    end       
    
    tmp_t_all = [];
    tmp_t_vel_all = [];
    tmp_yaw_all = [];
    tmp_fwd_all = [];
    tmp_ephys_all_A = [];
    tmp_ephys_all_B = [];

    % Load in all the DAQ data
    for i=1:length(dd)
        filepath = [datapath '/' dd{i}];
        
        cur_data = load(filepath);
        
        D = cur_data.trial_bdata;
        
        [currentA, voltageA, currentB, voltageB] = get_dual_scaled_voltage_and_current( D, get_B );
        
        t = cur_data.trial_time;
        
        [ t_vel, vel_forward, vel_side, vel_yaw ] = get_velocity_ephysrig(t, D, datapath, 1);
        
        tmp_t_all = vertcat( tmp_t_all, t );
        tmp_t_vel_all = vertcat( tmp_t_vel_all, t_vel' );
        
        vel_yaw = vel_yaw - mean(vel_yaw);
        
        tmp_yaw_all = vertcat( tmp_yaw_all, vel_yaw' );
        
        tmp_fwd_all = vertcat( tmp_fwd_all, vel_forward' );

%         voltageA = voltageA - mean(voltageA);
%         voltageB = voltageB - mean(voltageB);
        
        tmp_ephys_all_A = vertcat( tmp_ephys_all_A, voltageA );
        tmp_ephys_all_B = vertcat( tmp_ephys_all_B, voltageB );
    end
        
    if( TIME_CUTOFF_SECONDS == -1 )
        figure;
        subplot(2,1,1);
        plot(tmp_t_all, tmp_ephys_all_A);
        xlabel('Time (s)');
        ylabel('Voltage (mV)');
        title('Patch A');

        if (length(tmp_ephys_all_B) > 0 )
            subplot(2,1,2);
            plot(tmp_t_all, tmp_ephys_all_B);
            xlabel('Time (s)');
            ylabel('Voltage (mV)');
            title('Patch B');
        end
        
        waitforbuttonpress;
    else        
        TIME_CUTOFF = TIME_CUTOFF_SECONDS * ephys_SR;
        TIME_CUTOFF_VEL = TIME_CUTOFF_SECONDS * ball_SR;
        
        t_all{d}       = tmp_t_all(1:TIME_CUTOFF);
        ephys_all_A{d} = tmp_ephys_all_A(1:TIME_CUTOFF);
        
        if (length(tmp_ephys_all_B) > 0 )
            ephys_all_B{d} = tmp_ephys_all_B(1:TIME_CUTOFF);    
        end
        
        yaw_all{d}     = tmp_yaw_all(1:TIME_CUTOFF_VEL);
        fwd_all{d}     = tmp_fwd_all(1:TIME_CUTOFF_VEL);
        t_vel_all{d}   = tmp_t_vel_all(1:TIME_CUTOFF_VEL);
    end
end

end

