function display_spontaneous_vs_cx_evoked_turning_scatter_FR( basedir, cur_dirs, yaw_win_all, ephys_win_all, timebase_yaw, timebase_ephys )

set(0, 'DefaultFigureRenderer', 'painters');

FILT_FACTOR = 0.04;
ac = get_analysis_constants;

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;
BIN_SIZE = 0.050; % s
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

SHIFT_FACTOR = 3;

psth_dt_samples = ephys_SR/ball_SR;
psth_dt = psth_dt_samples / (1.0*ephys_SR);

PRE_PEAK_BUMP_RETURN_VEL = 0.75;
POST_PEAK_BUMP_RETURN_VEL = 0.0;

num_flies = length( cur_dirs );

YAW_DATA_IDX = 1;
FR_DATA_IDX = 2;
spontaneous_turning = cell( num_flies, 2 );
CX_turning          = cell( num_flies, 2 );

for i=1:num_flies
    for j=1:2
        spontaneous_turning{ i, j } = [];
        CX_turning{ i, j }          = [];
    end
end

analysis_path = [basedir '/CX_summary/'];

for d = 1:num_flies

    cur_datapath = cur_dirs{ d }{ 1 };
    cur_sid      = cur_dirs{ d }{ 2 };
    
    % datapath = [ basedir '/' cur_datapath ];
    % analysis_path = [datapath '/analysis/'];
    
    FR_THRESHOLD =  cur_dirs{ d }{ 3 };
    
    f = figure;
    
    for cond = 1:size( ephys_win_all,1 )
                
        
        % stim_filepath = [analysis_path '/stim_window_data_sid_' num2str( cur_sid ) '.mat'];
        
        % 'fwd_in_window', 'yaw_in_window', 'ephys_in_window', 'bump_in_window', 'VPS', 't_bump_w', 't_yaw_w', 't_ephys_w'
        %         cur_data = load( stim_filepath );
        %         bump_data   = cur_data.bump_in_window;
        %         fwd_data    = cur_data.fwd_in_window;
        %         yaw_data    = cur_data.yaw_in_window;
        %         ephys_data  = cur_data.ephys_in_window;
        %         t_bump_w    = cur_data.t_bump_w;
        %         t_yaw_w     = cur_data.t_yaw_w;
        %         t_ephys_w   = cur_data.t_ephys_w;
        %         VPS         = cur_data.VPS;
        
        cur_ephys_all          = ephys_win_all{ cond, d }(:,1:end-1);
        cur_yaw_all            = yaw_win_all{ cond, d }(:,1:end-1);
        
        cur_yaw_t              = timebase_yaw{ cond, d }(1:end-1);
        cur_ephys_t            = timebase_ephys{ cond, d }(1:end-1);
        
        for tr = 1 : size( cur_yaw_all, 1 )
            
            cur_ephys = squeeze( cur_ephys_all( tr, : ) );
            cur_yaw   = convert_yaw_to_degrees( squeeze( cur_yaw_all( tr, : ) ) );
            
            % cur_PSTH = calculate_PSTH_for_LAL_DN( cur_ephys_t, cur_yaw_t, cur_ephys, 'A2' );
            cur_PSTH = calculate_psth_debug( cur_ephys_t, cur_yaw_t, cur_ephys, ephys_SR, FR_THRESHOLD, psth_dt_samples, 0 );
            
            A2_PSTH_down = squeeze(mean(reshape( cur_PSTH, [ DT_YAW, length( cur_PSTH ) / DT_YAW ] ), 1 ));
            
            yaw_t_down = squeeze( mean( reshape( cur_yaw_t, [ DT_YAW, length( cur_yaw_t ) / DT_YAW ] ),1));
            yaw_down = squeeze( mean( reshape( cur_yaw, [ DT_YAW, length( cur_yaw_t ) / DT_YAW ] ),1));
            
            hold on;
            
            % Assuming the data is aligned to bump return velocity.
            time_of_peak = 0.0;                       
            CX_turn_start = time_of_peak - PRE_PEAK_BUMP_RETURN_VEL;
            CX_turn_end   = time_of_peak + POST_PEAK_BUMP_RETURN_VEL;                        
            
            % Plots spontaneous data
            cur_SPON_turning_t = find( ( yaw_t_down < CX_turn_start ) | ( yaw_t_down > CX_turn_end ) );
            cur_SPON_turning_t = cur_SPON_turning_t(SHIFT_FACTOR+1:end);
            
            if( length( spontaneous_turning{ d, YAW_DATA_IDX } ) == 0 )
                spontaneous_turning{ d, YAW_DATA_IDX } = yaw_down( cur_SPON_turning_t );
                spontaneous_turning{ d, FR_DATA_IDX } = A2_PSTH_down( cur_SPON_turning_t-SHIFT_FACTOR );                
            else
                spontaneous_turning{ d, YAW_DATA_IDX } = horzcat( spontaneous_turning{ d, YAW_DATA_IDX }, yaw_down( cur_SPON_turning_t ) );
                spontaneous_turning{ d, FR_DATA_IDX } = horzcat( spontaneous_turning{ d, FR_DATA_IDX }, A2_PSTH_down( cur_SPON_turning_t-SHIFT_FACTOR ) );
            end
            
            for ii = 1:( length(cur_SPON_turning_t ) - SHIFT_FACTOR )
                cur_index_yaw   = cur_SPON_turning_t( ii );
                cur_yaw_1 = yaw_down( cur_index_yaw + SHIFT_FACTOR );
                plot( A2_PSTH_down(ii), cur_yaw_1, 'o', 'MarkerSize', 3, 'color', 'b' );
            end
                                    
            % Plot CX-evoked turning epoch            
            cur_CX_turning_t = find( ( yaw_t_down >= CX_turn_start ) & ( yaw_t_down <= CX_turn_end ) );
            cur_CX_turning_t = cur_CX_turning_t(SHIFT_FACTOR+1:end);
            if( length( CX_turning{ d, YAW_DATA_IDX } ) == 0 )
                CX_turning{ d, YAW_DATA_IDX } = yaw_down( cur_CX_turning_t );
                CX_turning{ d, FR_DATA_IDX }  = A2_PSTH_down( cur_CX_turning_t-SHIFT_FACTOR );                
            else
                CX_turning{ d, YAW_DATA_IDX } = horzcat( CX_turning{ d, YAW_DATA_IDX }, yaw_down( cur_CX_turning_t ) );
                CX_turning{ d, FR_DATA_IDX }  = horzcat( CX_turning{ d, FR_DATA_IDX }, A2_PSTH_down( cur_CX_turning_t-SHIFT_FACTOR ) );
            end
                                    
            for ii = 1: ( length( cur_CX_turning_t ) - SHIFT_FACTOR )
                cur_index_yaw   = cur_CX_turning_t( ii );
                cur_yaw_2 = yaw_down( cur_index_yaw +SHIFT_FACTOR );
                plot( A2_PSTH_down( cur_index_yaw ), cur_yaw_2, 'o', 'MarkerSize', 10, 'color', 'r' );
            end
        end
    end
    
    
    [ spont_FO, G ] = fit( spontaneous_turning{ d, FR_DATA_IDX }', spontaneous_turning{ d, YAW_DATA_IDX }', 'poly1' );
    spon_rsq = G.rsquare;
    plt1 = plot( spont_FO );    
    set(plt1, 'Color', 'b');
    set(plt1, 'DisplayName', [ 'spontaneous rsq: ' num2str( spon_rsq ) ' slope: ' num2str( spont_FO.p1 ) ]);
    slopes_spon( d ) = spont_FO.p1;
    rsq_spon( d )    = spon_rsq;
    
    [ spon_FR_base(d,:), spon_mean(d,:), spon_95ci(d,:) ] = get_95_ci_for_spont_vs_cx( spontaneous_turning{ d, FR_DATA_IDX }, spontaneous_turning{ d, YAW_DATA_IDX } );
    
    [ CX_FO, G ] = fit( CX_turning{ d, FR_DATA_IDX }', CX_turning{ d, YAW_DATA_IDX }', 'poly1' );
    CX_rsq = G.rsquare;
    plt1 = plot( CX_FO );    
    set(plt1, 'Color', 'r');
    set(plt1, 'DisplayName', [ 'CX rsq: ' num2str( CX_rsq ) ' slope: ' num2str( CX_FO.p1 ) ]);
    slopes_CX( d ) = CX_FO.p1;
    rsq_CX( d )    = CX_rsq;

    [ CX_FR_base(d,:), CX_mean(d,:), CX_95ci(d,:) ] = get_95_ci_for_spont_vs_cx( CX_turning{ d, FR_DATA_IDX }, CX_turning{ d, YAW_DATA_IDX } );    
    
    tt = title(['Fly: ' num2str( d ) ' ' cur_datapath ]);
    set( tt, 'Interpreter', 'none' );

    xlabel('Firing rate ( spikes/s )');
    ylabel('Yaw (deg/s)');
    set(gca, 'FontSize', 14);
    
    legend();
    
    saveas(f,[analysis_path '/' cur_datapath '_FR_vs_yaw_spontaneous_CX_turning.fig']);
    saveas(f,[analysis_path '/' cur_datapath '_FR_vs_yaw_spontaneous_CX_turning.png']);
    saveas(f,[analysis_path '/' cur_datapath '_FR_vs_yaw_spontaneous_CX_turning.svg']);
    % close(f);
end

f = figure;

subplot(1,2,1);
hold on;
for d = 1:num_flies
    plot( [1, 2], [ slopes_spon(d), slopes_CX(d) ], '-o', 'DisplayName', ['Fly: ' num2str(d)]);
end

avg_slope_spon = mean( slopes_spon );
avg_slope_CX   = mean( slopes_CX );

PLOT_X_LINE_WIDTH = 0.25;
plot( [ 1-PLOT_X_LINE_WIDTH, 1+PLOT_X_LINE_WIDTH ], [ avg_slope_spon, avg_slope_spon ], 'color', 'b', 'LineWidth', 2.0 );
plot( [ 2-PLOT_X_LINE_WIDTH, 2+PLOT_X_LINE_WIDTH ], [ avg_slope_CX, avg_slope_CX ], 'color', 'b', 'LineWidth', 2.0 );

xlim([0 3]);
ylabel('slope of linear fit');
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Spontaneous', 'Central Complex'});
xtickangle( 45 );
legend();

subplot(1,2,2);
hold on;

for d = 1:num_flies
    plot( [1, 2], [ rsq_spon( d ), rsq_CX( d ) ], '-o', 'DisplayName', ['Fly: ' num2str(d)]);
end

avg_rsq_spon = mean( rsq_spon );
avg_rsq_CX   = mean( rsq_CX );

PLOT_X_LINE_WIDTH = 0.25;
plot( [ 1-PLOT_X_LINE_WIDTH, 1+PLOT_X_LINE_WIDTH ], [ avg_rsq_spon, avg_rsq_spon ], 'color', 'b', 'LineWidth', 2.0 );
plot( [ 2-PLOT_X_LINE_WIDTH, 2+PLOT_X_LINE_WIDTH ], [ avg_rsq_CX, avg_rsq_CX ], 'color', 'b', 'LineWidth', 2.0 );

xlim([0 3]);
ylabel('rsq of linear fit');
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Spontaneous', 'Central Complex'});
xtickangle( 45 );
legend();

saveas( f, [ analysis_path '/FR_vs_yaw_spontaneous_CX_turning_slope_rsq_per_fly.fig' ] );
saveas( f, [ analysis_path '/FR_vs_yaw_spontaneous_CX_turning_slope_rsq_per_fly.png' ] );
saveas( f, [ analysis_path '/FR_vs_yaw_spontaneous_CX_turning_slope_rsq_per_fly.svg' ] );

f = figure;

hold on;

% Plot spontaneous
for d = 1:num_flies
    cur_high_ci = spon_mean( d,: ) + spon_95ci( d,: );
    cur_low_ci  = spon_mean( d,: ) - spon_95ci( d,: );
    plt1 = plot( spon_FR_base(d,:), cur_high_ci, 'color', rgb('Gray') );    
    plot( spon_FR_base(d,:), cur_low_ci, 'color', rgb('Gray') );    
    plot( spon_FR_base(d,:), spon_mean(d,:), '--', 'color', rgb('Gray') );    
end

% Plot CX
for d = 1:num_flies
    cur_high_ci = CX_mean( d,: ) + CX_95ci( d,: );
    cur_low_ci  = CX_mean( d,: ) - CX_95ci( d,: );
    plt2 = plot( CX_FR_base(d,:), cur_high_ci, 'color', rgb('Red') );    
    plot( CX_FR_base(d,:), cur_low_ci, 'color', rgb('Red') );    
    plot( CX_FR_base(d,:), CX_mean(d,:), '--', 'color', rgb('Red') );    
end

xlabel('firing rate (spikes/s)');
ylabel('yaw velocity (deg/s)');
legend([plt1, plt2], 'Spontaneous','Central Complex');

saveas( f, [ analysis_path '/FR_vs_yaw_spontaneous_CX_turning_95ci.fig' ] );
saveas( f, [ analysis_path '/FR_vs_yaw_spontaneous_CX_turning_95ci.png' ] );
saveas( f, [ analysis_path '/FR_vs_yaw_spontaneous_CX_turning_95ci.svg' ] );

end

