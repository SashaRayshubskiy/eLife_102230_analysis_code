% Load all data from all the animals to be included in this analysis

clear all;
close all;

working_dir = '/data/drive1/sasha/';

settings = sensor_settings;
%ephys_SR = settings.sampRate;
ephys_SR = 4000;
ball_SR = settings.sensorPollFreq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directories_to_analyze =  { { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_03', 0, 1900 }, ...                          
                            { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_04', 1, 2346 }, ...                           
                            { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_05', 0, 1740 }, ... 
                            { '161210_lexAOpCsChrimson_gfp_83blexA_ss730_06', 0, 1882 } };

%directories_to_analyze = { { '161208_lexAOpCsChrimson_gfp_83blexA_ss730_05', 0, 1740 } };                        
%directories_to_analyze = { { '161218_lexAOpCsChrimson_gfp_83blexA_ss731_02', 0, 3000 } }; 

% directories_to_analyze = { { '161218_lexAOpCsChrimson_gfp_83blexA_ss731_02', 0, 3000 }, ...                          
%                            { '170721_lexAOpCsChrimson_gfp_83blexA_ss731_05', 2, 2829 }, ...                           
%                            { '170726_lexAOpCsChrimson_gfp_83blexA_ss731_08', 0, 1656 }, ...
%                            { '170727_lexAOpCsChrimson_gfp_83blexA_ss731_12', 0, 2800 }};

pre_stim_t = 3.0;
stim_t =  3.5;
post_stim_t = 3.0;
inter_trial_t = 5.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timenow_str = datestr(datetime, 'yymmdd_HHMMSS');
summary_analysis_path = [working_dir '/summary_analysis/ephys_vs_yaw_analysis/triggered_analysis_' timenow_str];

if(~exist(summary_analysis_path, 'dir'))
    mkdir(summary_analysis_path);
end

[t_all, t_vel_all, yaw_all, fwd_all, ephys_all_A, ephys_all_B] = load_LAL_DN_data( working_dir, directories_to_analyze, ephys_SR, ball_SR );

EVENT_FIND_PEAKS = 10;
EVENT_FR_THRESHOLD_CROSSING = 11;

EVENT_FINDING_METHOD = EVENT_FR_THRESHOLD_CROSSING;