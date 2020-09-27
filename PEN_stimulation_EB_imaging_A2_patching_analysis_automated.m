%% Load imaging and behavioral data

set( 0, 'DefaultFigureRenderer', 'painters' );
set( 0, 'DefaultAxesColor', 'none');

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

ctrl_directories = { { '181206_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_01', 0, 0.4, 0.5, 0.5 }, ... 
                     { '190208_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_02', 0, 0.4, 0.5, 0.5 }, ... 
                     { '190211_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_03', 1, 0.4, 0.5, 0.5 }, ... 
                     { '190212_Lex_6f_60D05_Gal4_P2X2_P2X2_recomb_04', 0, 0.4, 0.5, 0.5 } };
                 
excluded_flies = { {'181203_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16', 0 }, ...
                   {'181211_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_18', 1 }, ...
                   {'190218_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_25', 0 }, ...
                   {'190220_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_27', 1 }, ...
                   {'190221_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_28', 1 } };

%%
examine_fly_activity( basedir, excluded_flies );             
examine_fly_activity( basedir, exp_directories );             

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
     cur_dirs = { { '181205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_17', 0, 0.2, 0.5, 0.5, 0.5 } };
%      cur_dirs = { { '190131_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_19', 0, 0.15, 0.5, 0.5, 0.5 } };
%      cur_dirs = { { '190204_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_22', 0, 0.15, 0.5, 0.5, 0.5 } };
%       cur_dirs = { { '190205_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_23', 1, 0.15, 0.25, 2.5, 0.1 } }; 
end

%% Classify by bump returns up or down only 
[ bump_returns_up, bump_returns_down, no_response ] = filter_bump_returns_experiment_return_up_down_post_jump( basedir, cur_dirs );

%% For panel comparing bump return vs. not return with color coded panel by yaw
[ bump_returns_up, bump_returns_down, no_response, bump_tc_FF ] = filter_bump_returns_for_figure( basedir, cur_dirs );

%% DO NOT TOUCH: Classify by bump returns up or down only
bump_conditions = { bump_returns_up, bump_returns_down };
bump_conditions_str = { 'bump_returns_up', 'bump_returns_down' };

[ bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, Vm_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, ephys_win_all ] = align_by_bump_velocity_with_PSTH( basedir, cur_dirs, bump_conditions, bump_conditions_str );

%%
%%%%%%%%%%%%%%%%%%%%%%%%
% This unit is to explore bump magnitude with subsequent turn
%%%%%%%%%%%%%%%%%%%%%%%%%
bump_conditions = { bump_returns_up, bump_returns_down };
bump_conditions_str = { 'bump_returns_up', 'bump_returns_down' };

[ bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, Vm_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, ephys_win_all, full_fwd, full_yaw, full_bump, full_bump_t, full_ephys, full_yaw_t, full_ephys_t ] = align_by_bump_velocity_with_PSTH_with_full_trial_data( basedir, cur_dirs, bump_conditions, bump_conditions_str );

%% Examine parameters that were asked by Nature Neuroscience reviewers
examine_bump_turn_dynamics_NN( full_fwd, full_yaw, full_bump, full_bump_t, full_ephys, full_yaw_t, full_ephys_t, yaw_win_all, timebase_yaw );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Experimenting for figure
bump_conditions = { bump_returns_up, bump_returns_down };
bump_conditions_str = { 'bump_returns_up', 'bump_returns_down' };

[ bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, Vm_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, ephys_win_all, final_stims_passed_all_checks ] = align_by_bump_velocity_with_PSTH_for_figure( basedir, cur_dirs, bump_conditions, bump_conditions_str, bump_tc_FF );

%% 
% display_CX_summary_all_flies( bump_conditions, bump_conditions_str, bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, ephys_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, experiment_type_str );

display_CX_summary_all_flies_workaround_3( bump_conditions, bump_conditions_str, bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, Vm_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, experiment_type_str );

%% Display summary bar plots
% display_CX_summary_bar_plots_all_flies( bump_conditions, bump_conditions_str, bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, ephys_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, experiment_type_str );

display_CX_summary_bar_plots_all_flies_use_mean( bump_conditions, bump_conditions_str, bump_pos_win_all, bump_vel_win_all, yaw_win_all, fwd_win_all, Vm_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, experiment_type_str );


%% display spontaneous vs. CX-evoked turning scatter plot

%% Vm
display_spontaneous_vs_cx_evoked_turning_scatter( basedir, cur_dirs, yaw_win_all, ephys_win_all, timebase_yaw, timebase_ephys );

%%  FR
display_spontaneous_vs_cx_evoked_turning_scatter_FR( basedir, cur_dirs, yaw_win_all, ephys_win_all, timebase_yaw, timebase_ephys );

%% Breakout trials by jumps, and color code by weather the bump return or not
display_bump_jump_colorcode_returns( basedir, cur_dirs, final_stims_passed_all_checks );


%% Calculate experiment statistics

num_stims = [];

for d = 1:length( cur_dirs )
    
    cur_datapath = cur_dirs{d}{1};
    cur_sid      = cur_dirs{d}{2};
    
    datapath = [ basedir '/' cur_datapath ]; 
    analysis_path = [datapath '/analysis/'];
    stim_filepath = [analysis_path '/stim_window_data_sid_' num2str(cur_sid) '.mat'];
    
    cur_data = load( stim_filepath );    
    bump_data   = cur_data.bump_in_window;

    num_stims(d) = size( bump_data, 1);
    disp(['Stims for fly: ' num2str(d) ' = ' num2str(num_stims(d))]);
end

disp(['Avg: ' num2str( mean( num_stims )) ]);
disp(['SEM: ' num2str( get_sem( num_stims, 2 )) ]);

%% Controls 1: Display percent of bump jumps
experiment_type_str = 'experiment';
cur_dirs = exp_directories;
[exp_bump_jump_percent_per_fly] = calculate_percent_of_bump_jumps_v2( basedir, cur_dirs, experiment_type_str );

% Controls 1: Display percent of bump jumps
experiment_type_str = 'control';
cur_dirs = ctrl_directories;
[ctl_bump_jump_percent_per_fly] = calculate_percent_of_bump_jumps_v2( basedir, cur_dirs, experiment_type_str );

%% Controls 1: Display percent of bump jumps

%
f = figure;
SEM_DIM = 2;

%color_per_fly = { rgb('Violet'), rgb('Red'), rgb('Black'), rgb('Blue') };
color_per_fly = { rgb('Black'), rgb('Black'), rgb('Black'), rgb('Black') };

for cond = [1 2]
    
    if( cond == 1 )
        cur_data   = ctl_bump_jump_percent_per_fly * 100;
    else
        cur_data   = exp_bump_jump_percent_per_fly * 100;
    end
    
    for d = 1:size( cur_data, 2 )
        subplot( 1, 1, 1 );
        hold on;
        plot( cond, cur_data( d ), 'o', 'color', color_per_fly{d} )        
    end
    
    subplot( 1, 1, 1 );
    hold on;
    data_avg = mean( cur_data, SEM_DIM );
    plot( [cond-0.25, cond+0.25], [ data_avg, data_avg ]  );
    ylim([0 100]);
    ylabel('Percent jump');
    set(gca, 'XTick', [1 2]);
    set(gca, 'XTickLabel', {'Ctl', 'Exp'});
end

title('Percent jump');
filename = [ '/data/drive2/sasha/CX_summary' ];
saveas( f, [ filename '/percent_jump_summary.fig'] );
saveas( f, [ filename '/percent_jump_summary.png'] );

%% Controls 2: Display exp vs. control for post stim A2 response: examine the effect of ATP stimulation on A2
experiment_type_str = 'experiment';
cur_dirs = exp_directories;
[ exp_bump_response, exp_fwd_response, exp_yaw_response, exp_ephys_response ] = calculate_A2_post_stim_response( basedir, cur_dirs, experiment_type_str );

%
experiment_type_str = 'control';
cur_dirs = ctrl_directories;
[ ctl_bump_response, ctl_fwd_response, ctl_yaw_response, ctl_ephys_response ] = calculate_A2_post_stim_response( basedir, cur_dirs, experiment_type_str );

%%
f = figure;
SEM_DIM = 2;

%color_per_fly = { rgb('Violet'), rgb('Red'), rgb('Black'), rgb('Blue') };
color_per_fly = { rgb('Black'), rgb('Black'), rgb('Black'), rgb('Black') };

for cond = [1 2]
    
    if( cond == 1 )
        cur_yaw   = exp_yaw_response;
        cur_ephys = exp_ephys_response;
    else
        cur_yaw   = ctl_yaw_response;
        cur_ephys = ctl_ephys_response;
    end
    
    for d = 1:size( cur_yaw, 2 )
        subplot( 1, 2, 1 );
        hold on;
        plot( cond, cur_yaw( d ), 'o', 'color', color_per_fly{d} )
        
        subplot( 1, 2, 2 );
        hold on;
        plot( cond, cur_ephys( d ), 'o', 'color', color_per_fly{d} )
    end
    
    subplot( 1, 2, 1 );
    hold on;
    yaw_avg = mean( cur_yaw, SEM_DIM );
    yaw_sem = get_sem( cur_ephys, SEM_DIM );
    plot( [cond-0.25, cond+0.25], [ yaw_avg, yaw_avg ], 'DisplayName', [ 'avg: ' num2str( yaw_avg ) ' sem: ' num2str( yaw_sem ) ]  );
    ylim([-100 100]);
    ylabel('Yaw (deg/s)');    
    set(gca, 'XTickLabel', {'0', 'Exp', 'Ctrl', '3'})    
    
    subplot( 1, 2, 2 );
    hold on;
    ephys_avg = mean( cur_ephys, SEM_DIM );
    ephys_sem = get_sem( cur_ephys, SEM_DIM );
    plot( [cond-0.25, cond+0.25], [ ephys_avg, ephys_avg ], 'DisplayName', [ 'avg: ' num2str( ephys_avg ) ' sem: ' num2str( ephys_sem ) ] );  
    ylim([-2 2]);
    set(gca, 'XTickLabel', {'0', 'Exp', 'Ctrl', '3'})
    ylabel('Vm (mv)');        
end

subplot(1,2,1);
legend;clos
subplot(1,2,2);
legend;

title('Post ATP injection response');
filename = [ '/data/drive2/sasha/CX_summary' ];
saveas( f, [ filename '/post_stim_response_summary.fig'] );
saveas( f, [ filename '/post_stim_response_summary.png'] );

%% Test bump time courses

experiment_type_str = 'control';
cur_dirs = ctrl_directories;
[ctl_bump_jump_percent_per_fly] = calculate_percent_of_bump_jumps_v2( basedir, cur_dirs, experiment_type_str );




% ATTIC
%% Display experiment vs. control 

control_data = load('/data/drive2/sasha/CX_summary/control_data.mat');
experiment_data = load('/data/drive2/sasha/CX_summary/experiment_data.mat');

NUM_CONDITIONS = length( bump_conditions );
display_CX_summary_all_flies_experiment_vs_control( experiment_data, control_data, NUM_CONDITIONS, '/data/drive2/sasha/CX_summary/');


%% Control: Calculate velocity of spontaneous bump motion along with yaw and ephys. 
calculate_spontaneous_post_stim_bump_velocity( basedir, ctrl_directories );


%%
[ bump_jumps_up_returns_down, bump_jumps_down_returns_up, no_response ] = filter_bump_returns_experiment_v3( basedir, cur_dirs );
% bump_conditions = { bump_jumps_up_returns_down, bump_jumps_down_returns_up, no_response };
% bump_conditions_str = { 'bump_jumps_up_returns_down', 'bum....p_jumps_down_returns_up', 'no_response' };
%bump_conditions = { bump_jumps_down_returns_up, no_response };
%bump_conditions_str = { 'bump_jumps_down_returns_up', 'no_response' };
bump_conditions = { bump_jumps_down_returns_up, bump_jumps_up_returns_down };
bump_conditions_str = { 'bump_jumps_down_returns_up', 'bump_jumps_up_returns_down' };

[ bump_win_all, yaw_win_all, fwd_win_all, ephys_win_all, timebase_bump, timebase_yaw, timebase_ephys  ] = align_by_bump_velocity_v2( basedir, cur_dirs, bump_conditions, bump_conditions_str );

