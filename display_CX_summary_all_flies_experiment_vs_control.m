function display_CX_summary_all_flies_experiment_vs_control( exp_data, control_data, NUM_COND, analysis_path )

f = figure('units','normalized','outerposition',[0 0 1 1]);

all_data = { exp_data, control_data };

n_per_cond = zeros(2, NUM_COND);

for ed = 1:2
    
    cur_data = all_data{ ed };
    
    %for cond = 1:size( cur_data.bump_win_all, 1 )
    for cond = 1:1
        
        cur_cond_str = cur_data.bump_conditions_str{ cond };
        
        if(ed == 1)
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
        else
            if( ( strcmp(cur_cond_str, 'bump_jumps_up_returns_down') == 1 ) || ( strcmp(cur_cond_str, 'bump_returns_down') == 1 ) )
                sem_clr = rgb('Moccasin');
                avg_clr = rgb('Gold');
            elseif( ( strcmp(cur_cond_str, 'bump_jumps_down_returns_up') == 1 ) || ( strcmp(cur_cond_str, 'bump_returns_up') == 1 ) )
                sem_clr = rgb('LightSalmon');
                avg_clr = rgb('Red');
            elseif( strcmp(cur_cond_str, 'no_response') == 1 )
                sem_clr = rgb('SkyBlue');
                avg_clr = rgb('Navy');
            end
        end        
        
        bump_all = [];
        fwd_all = [];
        yaw_all = [];
        ephys_all = [];
        
        for d = 1:size( cur_data.bump_win_all, 2 )
            cur_bump = cur_data.bump_win_all{ cond, d };
            
            if( length(cur_bump) == 0 )
                continue;
            else
                bump_t  = cur_data.timebase_bump{ cond, d };
                yaw_t   = cur_data.timebase_yaw{ cond, d };
                ephys_t = cur_data.timebase_ephys{ cond, d };
                n_per_cond( ed, cond ) = n_per_cond( ed, cond ) + 1;
            end
            
            if( size(cur_bump,1) == 1 )
                bump_all(end+1,:) = cur_bump;
            else
                bump_all(end+1,:) = mean( cur_bump );
            end
            
            cur_fwd = cur_data.fwd_win_all{ cond, d };
            
            if( size(cur_fwd, 1 ) == 1)
                fwd_all(end+1,:) = cur_fwd;
            else
                fwd_all(end+1,:) = mean( cur_fwd );
            end
            
            cur_yaw = cur_data.yaw_win_all{ cond, d };
            
            if( size(cur_yaw, 1 ) == 1)
                yaw_all(end+1,:) = cur_yaw;
            else
                yaw_all(end+1,:) = mean( cur_yaw );
            end
            
            cur_ephys = cur_data.ephys_win_all{ cond, d };
            
            if( size(cur_ephys, 1 ) == 1)
                ephys_all(end+1,:) = cur_ephys;
            else
                ephys_all(end+1,:) = mean( cur_ephys );
            end
        end
        
        if( length( bump_all ) == 0 )
            continue;
        end
        
        ax1(1) = subplot(4,1,1);
        
        bump_sem = get_sem( bump_all, 1 );
        bump_avg = mean( bump_all );
        
        hold on;
        
        fh = fill( [bump_t, fliplr(bump_t)], ...
            [(bump_avg-bump_sem) fliplr((bump_avg+bump_sem))], ...
            sem_clr );
        set(fh, 'EdgeColor', 'None');
        
        pl(ed,cond) = plot( bump_t, bump_avg, 'color', avg_clr );
        
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
end

subplot(4,1,1);
if( size( pl,2 ) == 1 )
    ll = legend([pl(1,1), pl(2,1)], [ 'Exp: ' cur_data.bump_conditions_str{1}, ' (' num2str(n_per_cond( 1, 1 )) ')'], [ 'Cntrl: ' cur_data.bump_conditions_str{1}, ' (' num2str(n_per_cond( 2, 1 )) ')'] );
elseif( size( pl,2 ) == 2 )
    ll = legend([pl(1,1), pl(2,1), pl(1,2), pl(2,2)], ...
        [ 'Exp: ' cur_data.bump_conditions_str{1}, ' (' num2str(n_per_cond( 1, 1 )) ')'], [ 'Cntrl: ' cur_data.bump_conditions_str{1}, ' (' num2str(n_per_cond( 2, 1 )) ')'], ...
        [ 'Exp: ' cur_data.bump_conditions_str{2}, ' (' num2str(n_per_cond( 1, 2 )) ')'], [ 'Cntrl: ' cur_data.bump_conditions_str{2}, ' (' num2str(n_per_cond( 2, 2 )) ')'] );

elseif( size( pl,2 ) == 3 )
    ll = legend([pl(1,1), pl(2,1), pl(1,2), pl(2,2), pl(1,3), pl(2,3)], ...
        [ 'Exp: ' cur_data.bump_conditions_str{1}, ' (' num2str(n_per_cond( 1, 1 )) ')'], [ 'Cntrl: ' cur_data.bump_conditions_str{1}, ' (' num2str(n_per_cond( 2, 1 )) ')'], ...
        [ 'Exp: ' cur_data.bump_conditions_str{2}, ' (' num2str(n_per_cond( 1, 2 )) ')'], [ 'Cntrl: ' cur_data.bump_conditions_str{2}, ' (' num2str(n_per_cond( 2, 2 )) ')'], ... 
        [ 'Exp: ' cur_data.bump_conditions_str{3}, ' (' num2str(n_per_cond( 1, 3 )) ')'], [ 'Cntrl: ' cur_data.bump_conditions_str{3}, ' (' num2str(n_per_cond( 2, 3 )) ')'] );
end
set(ll, 'Interpreter', 'none');

saveas(f, [ analysis_path '/experiment_vs_control_all_flies_returns_up_only.fig' ]);
saveas(f, [ analysis_path '/experiment_vs_control_all_flies_returns_up_only.png' ]);

end

