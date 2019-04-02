function display_CX_summary_all_flies( bump_conditions, bump_conditions_str, bump_win_all, yaw_win_all, fwd_win_all, ephys_win_all, timebase_bump, timebase_yaw, timebase_ephys, experiment_type_str )

f = figure('units','normalized','outerposition',[0 0 1 1]);

n_per_cond = zeros(1, size( bump_win_all, 1 ));

for cond = 1:size( bump_win_all, 1 )

    cur_cond_str = bump_conditions_str{ cond };
    
    if( strcmp(cur_cond_str, 'bump_jumps_up_returns_down') == 1 ) 
        sem_clr = rgb('PaleGreen');
        avg_clr = rgb('SeaGreen');        
    elseif( strcmp(cur_cond_str, 'bump_jumps_down_returns_up') == 1 )  
        sem_clr = rgb('Violet');
        avg_clr = rgb('DarkMagenta');    
    elseif( strcmp(cur_cond_str, 'no_response') == 1 )  
        sem_clr = rgb('Silver');
        avg_clr = rgb('Black');    
    end
    
    bump_all = [];
    fwd_all = [];
    yaw_all = [];
    ephys_all = [];

    for d = 1:size( bump_win_all, 2 )
        cur_bump = bump_win_all{ cond, d };
        
        if( length(cur_bump) == 0 )
            continue;
        else
            bump_t  = timebase_bump{ cond, d };
            yaw_t   = timebase_yaw{ cond, d };
            ephys_t = timebase_ephys{ cond, d };
            n_per_cond( cond ) = n_per_cond( cond ) + 1;
        end
        
        bump_all(end+1,:) = mean( cur_bump );

        cur_fwd = fwd_win_all{ cond, d };
        fwd_all(end+1,:) = mean( cur_fwd );

        cur_yaw = yaw_win_all{ cond, d };
        yaw_all(end+1,:) = mean( cur_yaw );

        cur_ephys = ephys_win_all{ cond, d };
        ephys_all(end+1,:) = mean( cur_ephys );
    end
    
    if( length( bump_all ) == 0 )
        continue;
    end
    
    ax1(1) = subplot(4,1,1);
    
    bump_sem = get_sem(bump_all,1);
    bump_avg = mean( bump_all );        
    
    hold on;
    
    fh = fill( [bump_t, fliplr(bump_t)], ...
    [(bump_avg-bump_sem) fliplr((bump_avg+bump_sem))], ...
    sem_clr );
    set(fh, 'EdgeColor', 'None');
    pl(cond) = plot( bump_t, bump_avg, 'color', avg_clr );
    ylabel('EB bump velocity (wedge/s)');

    ax1(2) = subplot(4,1,2);
    fwd_sem = get_sem( fwd_all,1);
    fwd_avg = mean( fwd_all );        
    
    hold on;
    
    fh = fill( [yaw_t, fliplr(yaw_t)], ...
               [(fwd_avg-fwd_sem) fliplr((fwd_avg+fwd_sem))], ...
               sem_clr );

    set(fh, 'EdgeColor', 'None');
    plot( yaw_t, fwd_avg, 'color', avg_clr );
    ylabel('Fwd (au/s)');

    ax1(3) = subplot(4,1,3);
    yaw_sem = get_sem( yaw_all,1);
    yaw_avg = mean( yaw_all );        
    
    hold on;
    
    fh = fill( [ yaw_t, fliplr(yaw_t)], ...
    [(yaw_avg-yaw_sem) fliplr((yaw_avg+yaw_sem))], ...
    sem_clr );
    set(fh, 'EdgeColor', 'None');
    plot( yaw_t, yaw_avg, 'color', avg_clr );
    ylabel('Yaw (au/s)');

    ax1(4) = subplot(4,1,4);
    ephys_sem = get_sem(ephys_all,1);
    ephys_avg = mean( ephys_all );        
    
    hold on;
    
    fh = fill( [ephys_t, fliplr(ephys_t)], ...
    [(ephys_avg-ephys_sem) fliplr((ephys_avg+ephys_sem))], ...
    sem_clr );
    set(fh, 'EdgeColor', 'None');
    plot( ephys_t, ephys_avg, 'color', avg_clr );
    ylabel('Vm (mV)');
    xlabel('Time from max bump return vel (s)');
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

CX_summary_path = '/data/drive2/sasha/CX_summary/';

saveas(f, [ CX_summary_path '/' experiment_type_str '_bump_conditions_all_flies.fig' ]);
saveas(f, [ CX_summary_path '/' experiment_type_str '_bump_conditions_all_flies.png' ]);

end

