%% Load imaging and behavioral data

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sid 0, 2
% datapath = '/data/drive2/sasha/181203_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16/'; 

% sid 0
% datapath = '/data/drive2/sasha/181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17/'; 

% sid 1
% datapath = '/data/drive2/sasha/181211_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18/'; 

% sid 0
% datapath = '/data/drive2/sasha/190131_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19/'; 

% sid 0
% datapath = '/data/drive2/sasha/190204_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22/'; 

% No bump return up trials
%                     { '181211_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18', 1 }, ...
% Fly poorly moving { '181203_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16', 0, 0.2, 0.5, 0.5, 0.5 }, ...

% { path to data, sid, A2 FR threshold, JUMP_THRESHOLD, BUMP_RETURN_STABILITY_WINDOW, BUMP_SPEED_THRESHOLD }
basedir = '/data/drive2/sasha/';
exp_directories = { { '181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17', 0, 0.2, 0.5, 0.5, 0.5 }, ...
                    { '190131_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19', 0, 0.15, 0.5, 0.5, 0.5 }, ...
                    { '190204_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22', 0, 0.15, 0.5, 0.5, 0.5 }, ...
                    { '190205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23', 1, 0.15, 0.25, 2.5, 0.1 } };                           
                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sid 0
% datapath = '/data/drive2/sasha/181206_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_01/';

% sid 0
% datapath = '/data/drive2/sasha/190208_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_02/';

% sid 1
% datapath = '/data/drive2/sasha/190211_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_03/';

% sid 0
% datapath = '/data/drive2/sasha/190212_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_04/';

ctrl_directories = { { '181206_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_01', 0, 0.4 }, ... 
                     { '190208_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_02', 0, 0.4 }, ... 
                     { '190211_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_03', 1, 0.4 }, ... 
                     { '190212_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_04', 0, 0.4 } };

%%
extract_key_variables( basedir, exp_directories );                 
extract_key_variables( basedir, ctrl_directories );                 

% saves: raw_glom_ball_ephys.mat

                 
%% Extract all valid stimulation, this assumes key variables were extracted

extract_stim_windows( basedir, exp_directories );    
extract_stim_windows( basedir, ctrl_directories );    

% saves: stim_window_data_sid_<N>.mat


%% Experiment: Filter stims for clear bump returns for exp, calculate peak bump velocity
% and align it to ball and ephys 

RUN_EXPERIMENT = 55;
RUN_CONTROL    = 56;
RUN_TEST       = 57;

ANALYSIS_TYPE = RUN_EXPERIMENT; experiment_type_str = 'experiment';
% ANALYSIS_TYPE = RUN_CONTROL; experiment_type_str = 'control';
% ANALYSIS_TYPE = RUN_TEST; experiment_type_str = 'test';

if( ANALYSIS_TYPE == RUN_EXPERIMENT )
    cur_dirs = exp_directories;
elseif( ANALYSIS_TYPE == RUN_CONTROL )
    cur_dirs = ctrl_directories;
elseif( ANALYSIS_TYPE == RUN_TEST )    
%      cur_dirs = { { '181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17', 0, 0.2, 0.5, 0.5, 0.5 } };
%      cur_dirs = { { '190131_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19', 0, 0.15, 0.5, 0.5, 0.5 } };
%      cur_dirs = { { '190204_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22', 0, 0.15, 0.5, 0.5, 0.5 } };
      cur_dirs = { { '190205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23', 1, 0.15, 0.25, 2.5, 0.1 } }; 
end

%%
[ bump_jumps_up_returns_down, bump_jumps_down_returns_up, no_response ] = filter_bump_returns_experiment_v3( basedir, cur_dirs );
% bump_conditions = { bump_jumps_up_returns_down, bump_jumps_down_returns_up, no_response };
% bump_conditions_str = { 'bump_jumps_up_returns_down', 'bum....p_jumps_down_returns_up', 'no_response' };
%bump_conditions = { bump_jumps_down_returns_up, no_response };
%bump_conditions_str = { 'bump_jumps_down_returns_up', 'no_response' };
bump_conditions = { bump_jumps_down_returns_up, bump_jumps_up_returns_down };
bump_conditions_str = { 'bump_jumps_down_returns_up', 'bump_jumps_up_returns_down' };

[ bump_win_all, yaw_win_all, fwd_win_all, ephys_win_all, timebase_bump, timebase_yaw, timebase_ephys  ] = align_by_bump_velocity_v2( basedir, cur_dirs, bump_conditions, bump_conditions_str );

%% Classify by bump returns up or down only 

[ bump_returns_up, bump_returns_down, no_response ] = filter_bump_returns_experiment_return_up_down_post_jump( basedir, cur_dirs );

% Classify by bump returns up or down only
bump_conditions = { bump_returns_up, bump_returns_down };
bump_conditions_str = { 'bump_returns_up', 'bump_returns_down' };

% 
[ bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, ephys_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys  ] = align_by_bump_velocity_with_PSTH( basedir, cur_dirs, bump_conditions, bump_conditions_str );

%% 
% display_CX_summary_all_flies( bump_conditions, bump_conditions_str, bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, ephys_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, experiment_type_str );

display_CX_summary_all_flies_workaround_3( bump_conditions, bump_conditions_str, bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, ephys_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, experiment_type_str );

%% Display summary bar plots
% display_CX_summary_bar_plots_all_flies( bump_conditions, bump_conditions_str, bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, ephys_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, experiment_type_str );

display_CX_summary_bar_plots_all_flies_use_mean( bump_conditions, bump_conditions_str, bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, ephys_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, experiment_type_str );

%%
ANALYSIS_TYPE = RUN_EXPERIMENT; 
% ANALYSIS_TYPE = RUN_CONTROL; experiment_type_str = 'control';
% ANALYSIS_TYPE = RUN_TEST; experiment_type_str = 'test';

experiment_type_str = 'experiment';
cur_dirs = exp_directories;
[ exp_bump_response, exp_fwd_response, exp_yaw_response, exp_ephys_response ] = calculate_A2_post_stim_response( basedir, cur_dirs, experiment_type_str );

experiment_type_str = 'control';
cur_dirs = ctrl_directories;
[ ctl_bump_response, ctl_fwd_response, ctl_yaw_response, ctl_ephys_response ] = calculate_A2_post_stim_response( basedir, cur_dirs, experiment_type_str );




%% Display experiment vs. control 

control_data = load('/data/drive2/sasha/CX_summary/control_data.mat');
experiment_data = load('/data/drive2/sasha/CX_summary/experiment_data.mat');

NUM_CONDITIONS = length( bump_conditions );
display_CX_summary_all_flies_experiment_vs_control( experiment_data, control_data, NUM_CONDITIONS, '/data/drive2/sasha/CX_summary/');


%% Control: Calculate velocity of spontaneous bump motion along with yaw and ephys. 
calculate_spontaneous_post_stim_bump_velocity( basedir, ctrl_directories );





