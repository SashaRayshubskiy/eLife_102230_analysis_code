%% Load data
clear all;

global slash;
slash = '/';

working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/descending_neurons/';
%working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/osmotropotaxis_behavior_only/';

% Ephys data
data_type = 'ephys';
data_to_analyze = { {'170715_lexAOpCsChrimson_gfp_83blexA_ss730_ratio_optostim_02',0, 1000}, ...
    {'170718_lexAOpCsChrimson_gfp_83blexA_ss730_ratio_optostim_05', 2, 500 }, ...
    {'170719_lexAOpCsChrimson_gfp_83blexA_ss730_ratio_optostim_06', 0, 900}, ...
    {'170729_lexAOpCsChrimson_gfp_83blexA_ss730_ratio_optostim_09', 0, 1000}, ...
    {'170731_lexAOpCsChrimson_gfp_83blexA_ss730_ratio_optostim_10', 0, 1000}, ...
    {'170801_lexAOpCsChrimson_gfp_83blexA_ss730_ratio_optostim_12',0, 1000} };

% data_to_analyze = { {'170801_lexAOpCsChrimson_gfp_83blexA_ss730_ratio_optostim_12',0, 1000} };

% Behavioral only data
     % {'170725_lexAOpCsChrimson_gfp_83blexA_ss731_ratio_optostim_behavior_only_06', 0, 1000}, ...

% data_type = 'behavior';
% data_to_analyze = { {'170724_lexAOpCsChrimson_gfp_83blexA_ss730_ratio_optostim_behavior_only_02',0, 1000}, ...
%      {'170724_lexAOpCsChrimson_gfp_83blexA_ss731_ratio_optostim_behavior_only_04', 0, 600 }, ...
%      {'170724_lexAOpCsChrimson_gfp_83blexA_ss731_ratio_optostim_behavior_only_05', 0, 480}, ...
%      {'170725_lexAOpCsChrimson_gfp_83blexA_ss731_ratio_optostim_behavior_only_08', 0, 1000}, ... 
%      {'170725_lexAOpCsChrimson_gfp_83blexA_ss731_ratio_optostim_behavior_only_09', 0, 1000}, ...
%      {'170725_lexAOpCsChrimson_gfp_83blexA_ss731_ratio_optostim_behavior_only_10', 0, 1000}
%      };

%data_to_analyze = { {'170724_lexAOpCsChrimson_gfp_83blexA_ss731_ratio_optostim_behavior_only_05', 0, 480} };


analysis_path = [working_dir '/ratio_analysis/'];

if(~exist(analysis_path, 'dir'))
    mkdir(analysis_path);
end

trial_type_cnt = 10;

first_stim_t = 0.5;
last_stim_t =  1.0;

% L: { 1, 1, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, 0 }
% R: { 1, 0, 7/8, 3/4, 5/8, 1/2, 3/8, 1/4, 1/8, 1 }
da_cnt = length(data_to_analyze);

data_dump_path = [analysis_path '/' data_type '_data_to_analyze.mat'];
if( ~exist(data_dump_path, 'file') )
    
    t_volts = cell(1,da_cnt);
    volts_all = cell(1,da_cnt);
    t_vel_all = cell(1,da_cnt);
    yaw_vel_all = cell(1,da_cnt);
    forward_vel_all = cell(1,da_cnt);
    
    for d = 1:da_cnt
        t_volts{d} = cell(1,trial_type_cnt);
        volts_all{d} = cell(1,trial_type_cnt);
        t_vel_all{d} = cell(1,trial_type_cnt);
        yaw_vel_all{d} = cell(1,trial_type_cnt);
        forward_vel_all{d} = cell(1,trial_type_cnt);
    end
    
    for d = 1:da_cnt
        
        cur_dir = data_to_analyze{d}{1};
        sids = data_to_analyze{d}{2};
        TRIAL_CNT_MAX = data_to_analyze{d}{3};
        
        % Good data:
        % 170715_lexAOpCsChrimson_gfp_83blexA_ss730_ratio_optostim_02 sid = 0;
        % 170718_lexAOpCsChrimson_gfp_83blexA_ss730_ratio_optostim_05 sid = 0;
        
        datapath = [ working_dir cur_dir ];
        
        [ bdata_raw, bdata_time, trial_metadata ] = load_behavioral_data_ratio_optostim( sids, datapath, trial_type_cnt );
        
        TIME_OFFSET = bdata_time(1);
        
        % Transform raw data into velocity
        
        for tt = 1:trial_type_cnt
            for trial = 1:size(bdata_raw{tt},1)
                
                tid = trial_metadata{tt}(trial,2);
                if(tid > TRIAL_CNT_MAX)
                    continue;
                end
                
                t = bdata_time;
                D = squeeze(bdata_raw{tt}(trial,:,:));
                
                [currentA, voltageA, currentB, voltageB] = get_dual_scaled_voltage_and_current( D );
                
                
                volts_all{d}{tt}(trial, :) = voltageA;
                current_all{d}{tt}(trial, :) = currentA;
                %volts_all{tt}(trial, :) = filtfilt(d,D(:,2));
                t_volts{d}{tt}(trial,:) = t-TIME_OFFSET;
                
                [ t_vel, vel_forward, vel_side, vel_yaw ] = get_velocity_ephysrig(t, D, datapath, 1 );
                
                t_vel_all{d}{tt}( trial, : ) = t_vel-TIME_OFFSET;
                
                yaw_vel_all{d}{tt}( trial, : ) = vel_yaw;
                forward_vel_all{d}{tt}( trial, : ) = vel_forward + 5.0;
            end
        end
    end

    save( data_dump_path, 't_volts', 'volts_all', 't_vel_all', 'yaw_vel_all', 'forward_vel_all', 'data_to_analyze' );    
else    
    ld_data = load(data_dump_path);
    
    t_volts         = ld_data.t_volts;
    volts_all       = ld_data.volts_all;
    t_vel_all       = ld_data.t_vel_all;
    yaw_vel_all     = ld_data.yaw_vel_all;
    forward_vel_all = ld_data.forward_vel_all;      
    
    data_to_analyze = ld_data.data_to_analyze;      
end

%% Vm based

for d = 1:da_cnt
    cur_dir = data_to_analyze{d}{1};
    cur_analysis_path = [working_dir '/' cur_dir '/analysis' ];
    
    if(~exist(cur_analysis_path, 'dir'))
        mkdir(cur_analysis_path);
    end
    
    sids = data_to_analyze{d}{2};
    display_avg_data_3_LEFT_RIGHT_BOTH_ratio_optostim(t_volts{d}, volts_all{d}, t_vel_all{d}, yaw_vel_all{d}, forward_vel_all{d}, cur_analysis_path);    
end

%% PSTH based

for d = 1:da_cnt
    cur_dir = data_to_analyze{d}{1};
    cur_analysis_path = [working_dir '/' cur_dir '/analysis' ];
    
    if(~exist(cur_analysis_path, 'dir'))
        mkdir(cur_analysis_path);
    end
    
    sids = data_to_analyze{d}{2};
    display_avg_data_3_LEFT_RIGHT_BOTH_ratio_optostim_PSTH(t_volts{d}, volts_all{d}, t_vel_all{d}, yaw_vel_all{d}, forward_vel_all{d}, cur_analysis_path);    
end

%% Get ratio curve for ephys
display_ratio_curve_ephys(data_to_analyze, t_volts, volts_all, analysis_path);

%% Get ratio curve for ephys
display_ratio_curve_ephys_PSTH_based(data_to_analyze, t_volts, t_vel_all, volts_all, analysis_path);

%% Get ratio curve for ephys
display_ratio_curve_behavior_only(data_to_analyze, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path);

%%
display_ratio_curve_behavior_only_sigmoidal_fit(data_to_analyze, t_vel_all, yaw_vel_all, forward_vel_all, analysis_path);


