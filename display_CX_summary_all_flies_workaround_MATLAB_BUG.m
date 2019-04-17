function display_CX_summary_all_flies_workaround_MATLAB_BUG( bump_conditions, bump_conditions_str,  bump_pos_win_all, bump_win_all, ...
                                                             yaw_win_all, fwd_win_all, Vm_win_all, PSTH_win_all, ...
                                                             timebase_bump, timebase_yaw, timebase_Vm, experiment_type_str )
SHOW_FIRING_RATE = 11;
SHOW_VM          = 12;

% DISPAY_TYPE = SHOW_FIRING_RATE;
DISPLAY_TYPE = SHOW_VM;
                                                         
                                                    
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

f = figure('units','normalized','outerposition',[0 0 1 1]);

n_per_cond = zeros(1, size( bump_win_all, 1 ));

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

    for d = 1:size( bump_win_all, 2 )
        cur_bump_vel = bump_win_all{ cond, d };
        
        if( length(cur_bump_vel) == 0 )
            continue;
        else
            bump_t  = timebase_bump{ cond, d };
            yaw_t   = timebase_yaw{ cond, d };
            Vm_t = timebase_Vm{ cond, d };
            n_per_cond( cond ) = n_per_cond( cond ) + 1;
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

        cur_fwd = fwd_win_all{ cond, d };
        
        if( size(cur_fwd, 1 ) == 1)
            fwd_all(end+1,:) = cur_fwd;
        else
            fwd_all(end+1,:) = mean( cur_fwd );            
        end
        
        cur_yaw = yaw_win_all{ cond, d };
        
        if( size(cur_yaw, 1 ) == 1)
            yaw_all(end+1,:) = cur_yaw;
        else
            yaw_all(end+1,:) = mean( cur_yaw );            
        end

        cur_Vm = Vm_win_all{ cond, d };        
        if( size(cur_Vm, 1 ) == 1)
            Vm_all(end+1,:) = cur_Vm;
        else
            Vm_all(end+1,:) = mean( cur_Vm );            
        end
        
        cur_PSTH = PSTH_win_all{ cond, d };
        if( size(cur_PSTH, 1 ) == 1)
            PSTH_all(end+1,:) = cur_PSTH;
        else
            PSTH_all(end+1,:) = mean( cur_PSTH );            
        end
    end
    
    if( length( bump_vel_all ) == 0 )
        continue;
    end
    
    ax1(1) = subplot(5,1,1);

    bump_sem = get_sem(bump_pos_all,1);
    bump_avg = mean( bump_pos_all );        
    
    hold on;
    
    fh = fill( [bump_t, fliplr(bump_t)], ...
    [(bump_avg-bump_sem) fliplr((bump_avg+bump_sem))], ...
    sem_clr );
    set(fh, 'EdgeColor', 'None');
    pl(cond) = plot( bump_t, bump_avg, 'color', avg_clr );
    ylabel('EB bump position (wedge id)');
    
    ax1(2) = subplot(5,1,2);
    
    bump_sem = get_sem(bump_vel_all,1);
    bump_avg = mean( bump_vel_all );        
    
    hold on;
    
    fh = fill( [bump_t, fliplr(bump_t)], ...
    [(bump_avg-bump_sem) fliplr((bump_avg+bump_sem))], ...
    sem_clr );
    set(fh, 'EdgeColor', 'None');
    pl(cond) = plot( bump_t, bump_avg, 'color', avg_clr );
    ylabel('EB bump velocity (wedge/s)');

    ax1(3) = subplot(5,1,3);
    fwd_sem = get_sem( fwd_all,1);
    fwd_avg = mean( fwd_all );        
    
    hold on;
    
    fh = fill( [yaw_t, fliplr(yaw_t)], ...
               [(fwd_avg-fwd_sem) fliplr((fwd_avg+fwd_sem))], ...
               sem_clr );

    set(fh, 'EdgeColor', 'None');
    plot( yaw_t, fwd_avg, 'color', avg_clr );
    ylabel('Fwd (au/s)');

    ax1(4) = subplot(5,1,4);
    yaw_sem = get_sem( yaw_all,1);
    yaw_avg = mean( yaw_all );        
    
    hold on;
    
    fh = fill( [ yaw_t, fliplr(yaw_t)], ...
    [(yaw_avg-yaw_sem) fliplr((yaw_avg+yaw_sem))], ...
    sem_clr );
    set(fh, 'EdgeColor', 'None');
    plot( yaw_t, yaw_avg, 'color', avg_clr );
    ylabel('Yaw (au/s)');


    if( DISPLAY_TYPE == SHOW_VM )
        ax1(5) = subplot(5,1,5);
        Vm_sem = get_sem(Vm_all,1);
        Vm_avg = mean( Vm_all );
        
        hold on;
        
        fh = fill( [Vm_t, fliplr(Vm_t)], ...
            [(Vm_avg-Vm_sem) fliplr((Vm_avg+Vm_sem))], ...
            sem_clr );
        set(fh, 'EdgeColor', 'None');
        plot( Vm_t, Vm_avg, 'color', avg_clr );
        ylabel('Vm (mV)');
    else
        ax1(5) = subplot(5,1,5);
        PSTH_sem = get_sem(PSTH_all,1);
        PSTH_avg = mean( PSTH_all );
        
        hold on;
        
        fh = fill( [yaw_t, fliplr(yaw_t)], ...
            [(PSTH_avg-PSTH_sem) fliplr((PSTH_avg+PSTH_sem))], ...
            sem_clr );
        set(fh, 'EdgeColor', 'None');
        plot( yaw_t, PSTH_avg, 'color', avg_clr );
        ylabel('PSTH (spikes/s)');
    end
    
    xlabel('Time from max bump return vel (s)');
    xlim([-2 1.5]);
    linkaxes(ax1, 'x');
end

if( length( pl ) == 1 )
    ll = legend([pl(1)], [ bump_conditions_str{1}, ' (' num2str(n_per_cond( 1 )) ')'] );    
elseif( length( pl ) == 2 )
    ll = legend([pl(1), pl(2)], [ bump_conditions_str{1}, ' (' num2str(n_per_cond( 1 )) ')'], [ bump_conditions_str{2}, ' (' num2str(n_per_cond( 2 )) ')'] );
elseif( length( pl ) == 3 )
    ll = legend([pl(1), pl(2), pl(3)], [ bump_conditions_str{1}, ' (' num2str(n_per_cond( 1 )) ')'], [ bump_conditions_str{2}, ' (' num2str(n_per_cond( 2 )) ')'], [ bump_conditions_str{3}, ' (' num2str(n_per_cond( 3 )) ')'] );
end
set(ll, 'Interpreter', 'none');

if( DISPLAY_TYPE == SHOW_VM )
    analysis_str = 'Vm';
else
    analysis_str = 'FR';
end

saveas(f, [ CX_summary_path '/' experiment_type_str '_bump_conditions_all_flies_' analysis_str '.fig' ]);
saveas(f, [ CX_summary_path '/' experiment_type_str '_bump_conditions_all_flies_' analysis_str '.png' ]);

end

