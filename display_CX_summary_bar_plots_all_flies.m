function display_CX_summary_bar_plots_all_flies( bump_conditions, bump_conditions_str,  bump_pos_win_all, bump_win_all, ...
                                                 yaw_win_all, fwd_win_all, Vm_win_all, PSTH_win_all, ...
                                                 timebase_bump, timebase_yaw, timebase_Vm, experiment_type_str )

save_params.bump_conditions = bump_conditions;
save_params.bump_conditions_str = bump_conditions_str;

save_params.timebase_bump = timebase_bump;
save_params.bump_pos_win_all = bump_pos_win_all;
save_params.bump_win_all = bump_win_all;

save_params.timebase_yaw = timebase_yaw;
save_params.fwd_win_all = fwd_win_all;
save_params.yaw_win_all = yaw_win_all;

save_params.timebase_Vm = timebase_Vm;
save_params.Vm_win_all = Vm_win_all;
save_params.PSTH_win_all = PSTH_win_all;

CX_summary_path = '/data/drive2/sasha/CX_summary/';
save([ CX_summary_path '/' experiment_type_str '_data.mat'], '-struct', 'save_params' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BUMP_VEL_RANGE_MIN = -0.5;
BUMP_VEL_RANGE_MAX =  0.5;

YAW_VEL_RANGE_MIN = -0.5;
YAW_VEL_RANGE_MAX =  0;

EPHYS_RANGE_MIN = -0.5;
EPHYS_RANGE_MAX =  0;

EPHYS_BASELINE_RANGE_MIN = -1.0;
EPHYS_BASELINE_RANGE_MAX = -0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure('units','normalized','outerposition',[0 0 1 1]);

n_per_cond = zeros(1, size( bump_win_all, 1 ));

bump_vel_param_in_range = cell( 1, size( bump_win_all, 1 ) );
yaw_vel_param_in_range  = cell( 1, size( bump_win_all, 1 ) );

Vm_param_in_range       = cell( 1, size( bump_win_all, 1 ) );
FR_param_in_range       = cell( 1, size( bump_win_all, 1 ) );

for cond = 1:size( bump_win_all, 1 )

    cur_cond_str = bump_conditions_str{ cond };
      
    if( ( strcmp(cur_cond_str, 'bump_jumps_up_returns_down') == 1 ) || ( strcmp(cur_cond_str, 'bump_returns_down') == 1 ) )
        sem_clr = rgb('PaleGreen');
        avg_clr = rgb('SeaGreen');        
    elseif( ( strcmp(cur_cond_str, 'bump_jumps_down_returns_up') == 1 ) || ( strcmp(cur_cond_str, 'bump_returns_up') == 1 ) )
        sem_clr = rgb('Violet');
        avg_clr = rgb('DarkMagenta');    
    elseif( strcmp(cur_cond_str, 'no_response') == 1 )  
        sem_clr = rgb('Silver');
        avg_clr = rgb('Black');
    end
    
    bump_pos_all = [];
    bump_vel_all = [];
    fwd_all = [];
    yaw_all = [];
    Vm_all = [];
    PSTH_all = [];

    bump_vel_param_in_range{ cond } = [];
    yaw_vel_param_in_range{ cond } = [];
    
    Vm_param_in_range{ cond } = [];
    FR_param_in_range{ cond } = [];   

    for d = 1:size( bump_win_all, 2 )
        cur_bump_vel = bump_win_all{ cond, d };
        
        if( length(cur_bump_vel) == 0 )
            continue;
        else
            n_per_cond( cond ) = n_per_cond( cond ) + 1;
            bump_t  = timebase_bump{ cond, d };
            yaw_t   = timebase_yaw{ cond, d };
            Vm_t = timebase_Vm{ cond, d };
            
            bump_vel_in_range_idx = find( ( bump_t >= BUMP_VEL_RANGE_MIN ) &  (  bump_t <= BUMP_VEL_RANGE_MAX ) );
            yaw_vel_in_range_idx = find( ( yaw_t >= YAW_VEL_RANGE_MIN ) &  (  yaw_t <= YAW_VEL_RANGE_MAX ) );
            yaw_vel_baseline_in_range_idx = find( ( yaw_t >= EPHYS_BASELINE_RANGE_MIN ) &  (  yaw_t <= EPHYS_BASELINE_RANGE_MAX ) );
            ephys_in_range_idx = find( ( Vm_t >= EPHYS_RANGE_MIN ) &  (  Vm_t <= EPHYS_RANGE_MAX ) );
            ephys_baseline_in_range_idx = find( ( Vm_t >= EPHYS_BASELINE_RANGE_MIN ) &  (  Vm_t <= EPHYS_BASELINE_RANGE_MAX ) );
        end
        
        % Bump position
        cur_bump_pos = bump_pos_win_all{ cond, d };
        if( size(cur_bump_pos,1) == 1 )
            bump_pos_all(end+1,:) = cur_bump_pos;
        else
            bump_pos_all(end+1,:) = mean( cur_bump_pos );
        end
        
        % Bump velocity
        if( size(cur_bump_vel,1) == 1 )
            bump_vel_all(end+1,:) = cur_bump_vel;
        else
            bump_vel_all(end+1,:) = mean( cur_bump_vel );
        end
        
        avg_bump_vel_for_cur_fly = squeeze( bump_vel_all(end,:) );
        
        if( cond == 1 )
            bump_vel_param_in_range{cond}(end+1) = max( avg_bump_vel_for_cur_fly( bump_vel_in_range_idx ));
        else
            bump_vel_param_in_range{cond}(end+1) = min( avg_bump_vel_for_cur_fly( bump_vel_in_range_idx ));
        end
               
        % Fwd velocity
        cur_fwd = fwd_win_all{ cond, d };
        
        if( size(cur_fwd, 1 ) == 1)
            fwd_all(end+1,:) = cur_fwd;
        else
            fwd_all(end+1,:) = mean( cur_fwd );            
        end
        
        % Yaw velocity
        cur_yaw = convert_yaw_to_degrees( yaw_win_all{ cond, d } );
        
        if( size(cur_yaw, 1 ) == 1)
            yaw_all(end+1,:) = cur_yaw;
        else
            yaw_all(end+1,:) = mean( cur_yaw );            
        end        
        
        avg_yaw_for_cur_fly = squeeze( yaw_all(end,:) );
        
        if( cond == 1 )
            yaw_vel_param_in_range{cond}(end+1) = min( avg_yaw_for_cur_fly( yaw_vel_in_range_idx ));
        else
            yaw_vel_param_in_range{cond}(end+1) = max( avg_yaw_for_cur_fly( yaw_vel_in_range_idx ));
        end
                
        % Vm
        cur_Vm = Vm_win_all{ cond, d };        
        if( size(cur_Vm, 1 ) == 1)
            Vm_all(end+1,:) = cur_Vm;
        else
            Vm_all(end+1,:) = mean( cur_Vm );            
        end
        avg_Vm_for_cur_fly = squeeze( Vm_all(end,:) );
        Vm_param_in_range{cond}(end+1) = max( avg_Vm_for_cur_fly( ephys_in_range_idx )) - mean( avg_Vm_for_cur_fly( ephys_baseline_in_range_idx ));
        
        % Firing rate
        cur_PSTH = PSTH_win_all{ cond, d };
        if( size(cur_PSTH, 1 ) == 1)
            PSTH_all(end+1,:) = cur_PSTH;
        else
            PSTH_all(end+1,:) = mean( cur_PSTH );            
        end
        avg_PSTH_for_cur_fly = squeeze( PSTH_all(end,:) );
        FR_param_in_range{cond}(end+1) = max( avg_PSTH_for_cur_fly( yaw_vel_in_range_idx )) - mean( avg_PSTH_for_cur_fly( yaw_vel_baseline_in_range_idx ));
    end
    
    if( length( bump_vel_all ) == 0 )
        continue;
    end
    
    ax1(1) = subplot(1,4,1);

    bump_sem = get_sem( bump_vel_param_in_range{cond}, 1 );
    bump_avg = mean( bump_vel_param_in_range{cond} );
    
    ERROR_BAR_CLR = rgb('Black');
    hold on;
        
    pl(cond) = bar( cond, bump_avg, 'FaceColor', avg_clr );
    errorbar( cond, bump_avg, bump_sem, 'color', ERROR_BAR_CLR );
    ylabel('EB bump velocity (wed/s)');
    
    ax1(2) = subplot(1,4,2);
    
    yaw_sem = get_sem( yaw_vel_param_in_range{cond}, 1 );
    yaw_avg = mean( yaw_vel_param_in_range{cond} );
    
    hold on;
    
    bar( cond, yaw_avg, 'FaceColor', avg_clr );
    errorbar( cond, yaw_avg, yaw_sem, 'color', ERROR_BAR_CLR );
    ylabel('Yaw velocity (deg/s)');

    ax1(3) = subplot(1,4,3);
    Vm_sem = get_sem( Vm_param_in_range{ cond }, 1 );
    Vm_avg = mean( Vm_param_in_range{ cond } );
    
    hold on;
    
    bar( cond, Vm_avg, 'FaceColor', avg_clr );
    errorbar( cond, Vm_avg, Vm_sem, 'color', ERROR_BAR_CLR );
    ylabel('Delta Vm (mV)');
    
    ax1(4) = subplot(1,4,4);
    FR_sem = get_sem( FR_param_in_range{ cond }, 1 );
    FR_avg = mean( FR_param_in_range{ cond } );
    
    hold on;
    
    bar( cond, FR_avg, 'FaceColor', avg_clr );
    errorbar( cond, FR_avg, FR_sem, 'color', ERROR_BAR_CLR );
    ylabel('Delta firing rate (spikes/s)');
end

if( length( pl ) == 1 )
    ll = legend([pl(1)], [ bump_conditions_str{1}, ' (' num2str(n_per_cond( 1 )) ')'] );    
elseif( length( pl ) == 2 )
    ll = legend([pl(1), pl(2)], [ bump_conditions_str{1}, ' (' num2str(n_per_cond( 1 )) ')'], [ bump_conditions_str{2}, ' (' num2str(n_per_cond( 2 )) ')'] );
elseif( length( pl ) == 3 )
    ll = legend([pl(1), pl(2), pl(3)], [ bump_conditions_str{1}, ' (' num2str(n_per_cond( 1 )) ')'], [ bump_conditions_str{2}, ' (' num2str(n_per_cond( 2 )) ')'], [ bump_conditions_str{3}, ' (' num2str(n_per_cond( 3 )) ')'] );
end
set(ll, 'Interpreter', 'none');

saveas(f, [ CX_summary_path '/' experiment_type_str '_bump_conditions_all_bar_plots_flies.fig' ]);
saveas(f, [ CX_summary_path '/' experiment_type_str '_bump_conditions_all_bar_plots_flies.png' ]);

end

