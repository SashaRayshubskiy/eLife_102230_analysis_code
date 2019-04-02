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

basedir = '/data/drive2/sasha/';
exp_directories = { { '181203_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16', 0 }, ...
                    { '181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17', 0 }, ...
                    { '181211_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18', 1 }, ...
                    { '190131_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19', 0 }, ...
                    { '190204_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22', 0 } };                           
                             
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

ctrl_directories = { { '181206_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_01', 0 }, ... 
                     { '190208_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_02', 0 }, ... 
                     { '190211_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_03', 1 }, ... 
                     { '190212_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_04', 0 } };

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


% [ bump_jumps_up_returns_down, bump_jumps_down_returns_up, no_response ] = filter_bump_returns_experiment( basedir, test_dir );

RUN_EXPERIMENT = 55;
RUN_CONTROL    = 56;
RUN_TEST       = 57;

% ANALYSIS_TYPE = RUN_EXPERIMENT; experiment_type_str = 'experiment';
ANALYSIS_TYPE = RUN_CONTROL; experiment_type_str = 'control';
% ANALYSIS_TYPE = RUN_TEST; experiment_type_str = 'test';

if( ANALYSIS_TYPE == RUN_EXPERIMENT )
    cur_dirs = exp_directories;
elseif( ANALYSIS_TYPE == RUN_CONTROL )
    cur_dirs = ctrl_directories;
elseif( ANALYSIS_TYPE == RUN_TEST )
    % cur_dirs = { { '181203_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16', 0 } }; 
    cur_dirs = { { '181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17', 0 } }; 
end

[ bump_jumps_up_returns_down, bump_jumps_down_returns_up, no_response ] = filter_bump_returns_experiment_v2( basedir, cur_dirs );

%% 
% bump_conditions = { bump_jumps_up_returns_down, bump_jumps_down_returns_up, no_response };
% bump_conditions_str = { 'bump_jumps_up_returns_down', 'bum....p_jumps_down_returns_up', 'no_response' };
%bump_conditions = { bump_jumps_down_returns_up, no_response };
%bump_conditions_str = { 'bump_jumps_down_returns_up', 'no_response' };
bump_conditions = { bump_jumps_down_returns_up, bump_jumps_up_returns_down };
bump_conditions_str = { 'bump_jumps_down_returns_up', 'bump_jumps_up_returns_down' };

[ bump_win_all, yaw_win_all, fwd_win_all, ephys_win_all, timebase_bump, timebase_yaw, timebase_ephys  ] = align_by_bump_velocity( basedir, cur_dirs, bump_conditions, bump_conditions_str );

%% 
display_CX_summary_all_flies( bump_conditions, bump_conditions_str, bump_win_all, yaw_win_all, fwd_win_all, ephys_win_all, timebase_bump, timebase_yaw, timebase_ephys, experiment_type_str );


%% Control: Calculate velocity of spontaneous bump motion along with yaw and ephys. 
calculate_spontaneous_post_stim_bump_velocity( basedir, ctrl_directories );





