%% Load data
clear all;

global slash;
slash = '/';

working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/descending_neurons/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A2 data
% { data_directory, sid, number_of_usable_trials }
% LAL_DN_Type = 'A2';
% experiment_data = { { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_03', 0, 165 }, ...
%                     { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_04', 1, 204 }, ...
%                     { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_05', 0, 150 }, ...
%                     { '161210_lexAOpCsChrimson_gfp_83blexA_ss730_06', 0, 163 } };
% 
% control_data = { { '170808_gfp_ss730_08', 0, 300 }, ... 
%                  { '170808_gfp_ss730_07', 0, 200 }, ...
%                  { '170806_gfp_ss730_05', 0, 300 }, ...
%                  { '170804_gfp_ss730_01', 0, 300 } };
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A1 data
% { data_directory, sid, number_of_usable_trials }
LAL_DN_Type = 'A1';
experiment_data = { { '170727_lexAOpCsChrimson_gfp_83blexA_ss731_12', 0, 240 }, ...  
                    { '170726_lexAOpCsChrimson_gfp_83blexA_ss731_08', 0, 135 }, ...
                    { '161218_lexAOpCsChrimson_gfp_83blexA_ss731_02', 0, 258 } };

control_data = { { '170802_gfp_ss731_03', 0, 300 }, ... 
                 { '170803_gfp_ss731_06', 0, 300 }, ...
                 { '170803_gfp_ss731_07', 0, 240 }, ...
                 { '170804_gfp_ss731_08', 0, 100 } };
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             

all_data = { experiment_data, control_data };             

analysis_dir = [working_dir '/summary_analysis/' LAL_DN_Type '_osmotropotaxis_vs_control/'];

data_file_mat = [analysis_dir 'all_data.mat' ];

first_stim_t = 3.0;
last_stim_t =  3.5;

trial_type_cnt = 3;

if(~exist(data_file_mat))
    
    t_volts_all = {};
    volts_all = {};
    
    t_vel_all = {};
    yaw_vel_all = {};
    forward_vel_all = {};
    
    for e_type = 1:2
        experiments_for_condition_cnt = length(all_data{e_type});
        t_volts_all{e_type}     = cell(1,experiments_for_condition_cnt);
        volts_all{e_type}       = cell(1,experiments_for_condition_cnt);
        
        t_vel_all{e_type}       = cell(1,experiments_for_condition_cnt);
        yaw_vel_all{e_type}     = cell(1,experiments_for_condition_cnt);
        forward_vel_all{e_type} = cell(1,experiments_for_condition_cnt);
    end
    
    for e_type = 1:2
        
        cur_experiment_type = all_data{e_type};
        for dir_idx = 1:length(cur_experiment_type)
            cur_filename        = cur_experiment_type{dir_idx}{ 1 };
            cur_sid             = cur_experiment_type{dir_idx}{ 2 };
            cur_max_trial_count = cur_experiment_type{dir_idx}{ 3 };
            
            cur_datapath = [ working_dir cur_filename ];
            
            TRIAL_CNT_MAX = cur_max_trial_count;
            [ bdata_raw, bdata_time, trial_metadata ] = load_behavioral_data( cur_sid, cur_datapath, trial_type_cnt );
            
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
                    
                    volts_all{e_type}{dir_idx}{tt}(trial, :) = voltageA;
                    current_all{e_type}{dir_idx}{tt}(trial, :) = currentA;
                    %volts_all{tt}(trial, :) = filtfilt(d,D(:,2));
                    t_volts{e_type}{dir_idx}{tt}(trial,:) = t-TIME_OFFSET;
                    
                    [ t_vel, vel_forward, vel_side, vel_yaw ] = get_velocity_ephysrig(t, D, cur_datapath, 1 );
                    
                    t_vel_all{e_type}{dir_idx}{tt}( trial, : ) = t_vel-TIME_OFFSET;
                    
                    yaw_vel_all{e_type}{dir_idx}{tt}( trial, : ) = vel_yaw;
                    forward_vel_all{e_type}{dir_idx}{tt}( trial, : ) = vel_forward;
                end
            end
        end
    end
    
    save(data_file_mat,'volts_all', 'current_all', 't_volts', 't_vel_all', 'yaw_vel_all', 'forward_vel_all', 'all_data');
else
    dd = load(data_file_mat);
    
    volts_all        = dd.volts_all;
    current_all      = dd.current_all;
    t_volts          = dd.t_volts;
    t_vel_all        = dd.t_vel_all;
    yaw_vel_all      = dd.yaw_vel_all;
    forward_vel_all  = dd.forward_vel_all;
    all_data         = dd.all_data;
end

%% 
display_osmotropotaxis_experiment_vs_control( t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_dir );

%%
display_osmotropotaxis_experiment_vs_control_PSTH( LAL_DN_Type, t_volts, volts_all, t_vel_all, yaw_vel_all, forward_vel_all, analysis_dir );
