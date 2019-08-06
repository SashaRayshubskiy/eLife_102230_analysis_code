function [ bump_pos_win_all, bump_win_all, yaw_win_all, fwd_win_all, Vm_win_all, PSTH_win_all, timebase_bump, timebase_yaw, timebase_ephys, ephys_win_all, final_stims_passed_all_checks ] = align_by_bump_velocity_with_PSTH_for_figure( basedir, directories, bump_conditions, bump_conditions_str, bump_tc_FF )

set( 0, 'DefaultFigureRenderer', 'painters' );
set( 0, 'DefaultAxesColor', 'none' );

settings = sensor_settings;
BALL_FR = settings.sensorPollFreq;
EPHYS_FR = settings.sampRate;

psth_dt_samples = EPHYS_FR/BALL_FR;
SPIKE_THRESHOLD_A2 = 0.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BUMP_CONDITION_RETURNED_UP       = 33;
BUMP_CONDITION_RETURNED_DOWN     = 34;
BUMP_CONDITION_NOT_RETURNED_UP   = 35;
BUMP_CONDITION_NOT_RETURNED_DOWN = 36;
BUMP_CONDITION_NO_RESPONSE       = 37;
BUMP_CONDITION_UNDEFINED         = 38;

DEBUG_VERBOSE                    = 77;
DEBUG_OFF                        = 78;
DEBUG_LEVEL                      = DEBUG_VERBOSE;

TIME_BEFORE_EB_VEL_CHANGE        = 2.5; % Used to be 1.0 s
TIME_AFTER_EB_VEL_CHANGE         = 2.0;

% BUMP_SPEED_THRESHOLD           = 0.5; % wedges/s CHANGED ON 4/12/2019
BUMP_SPEED_THRESHOLD_MIN         = 0.1; % wedges/s
BUMP_SPEED_THRESHOLD_MAX         = 20; % wedges/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bump_pos_win_all              = cell( length(bump_conditions), length( directories ));
bump_win_all                  = cell( length(bump_conditions), length( directories ));
yaw_win_all                   = cell( length(bump_conditions), length( directories ));
fwd_win_all                   = cell( length(bump_conditions), length( directories ));
Vm_win_all                    = cell( length(bump_conditions), length( directories ));
PSTH_win_all                  = cell( length(bump_conditions), length( directories ));
timebase_bump                 = cell( length(bump_conditions), length( directories ));
timebase_yaw                  = cell( length(bump_conditions), length( directories ));
timebase_ephys                = cell( length(bump_conditions), length( directories ));
final_stims_passed_all_checks = cell( length(bump_conditions), length( directories ));

ephys_win_all          = cell( length(bump_conditions), length( directories ));

bump_and_yaw_fig                       = figure;
bump_trial_to_trial_colored_by_yaw_fig = figure;

% for cond = 1
for cond = 1:length(bump_conditions)
    cur_cond     = bump_conditions{cond};
    cur_cond_str = bump_conditions_str{cond};
    
    BUMP_CONDITION = BUMP_CONDITION_UNDEFINED;
    if( ( strcmp(cur_cond_str, 'bump_jumps_up_returns_down') == 1 ) || ( strcmp(cur_cond_str, 'bump_returns_down') == 1 ) )
        BUMP_CONDITION = BUMP_CONDITION_RETURNED_DOWN;
    elseif( ( strcmp(cur_cond_str, 'bump_jumps_down_returns_up') == 1 ) || ( strcmp(cur_cond_str, 'bump_returns_up') == 1 ) )
        BUMP_CONDITION = BUMP_CONDITION_RETURNED_UP;
    elseif( strcmp(cur_cond_str, 'no_response') == 1 )  
        BUMP_CONDITION = BUMP_CONDITION_NO_RESPONSE;
    elseif( strcmp(cur_cond_str, 'bump_not_returned_up') == 1 )  
        BUMP_CONDITION = BUMP_CONDITION_NOT_RETURNED_UP;
    elseif( strcmp(cur_cond_str, 'bump_not_returned_down') == 1 )  
        BUMP_CONDITION = BUMP_CONDITION_NOT_RETURNED_DOWN;
    end
        
    for d = 1:length( directories )
    % for d = 1
        
        BUMP_SPEED_THRESHOLD_MIN = directories{ d }{ 6 };
    
        cur_datapath = directories{ d }{ 1 };
        cur_sid      = directories{ d }{ 2 };
        
        datapath = [ basedir '/' cur_datapath ];
        analysis_path = [datapath '/analysis/'];

        t_now = datetime('now');
        t_now_analysis_path = [datapath '/analysis/' datestr(t_now)];
        mkdir(t_now_analysis_path);
        
        cur_cond_bump_tc     = cur_cond{ d, 1 };
        cur_cond_stim_ids    = cur_cond{ d, 2 };
        cur_cond_bump_tc_FF  = bump_tc_FF{cond};
        
              
        if( ( BUMP_CONDITION == BUMP_CONDITION_RETURNED_UP ) || ( BUMP_CONDITION == BUMP_CONDITION_RETURNED_DOWN ) )
            cur_cond_bump_return_idx = cur_cond{d,3};
        end
        
        stim_filepath = [analysis_path '/stim_window_data_sid_' num2str(cur_sid) '.mat'];

        % 'fwd_in_window', 'yaw_in_window', 'ephys_in_window', 'bump_in_window', 'VPS', 't_bump_w', 't_yaw_w', 't_ephys_w' 
        cur_data = load( stim_filepath );
        bump_data   = cur_data.bump_in_window;
        fwd_data    = cur_data.fwd_in_window;
        yaw_data    = cur_data.yaw_in_window;
        ephys_data  = cur_data.ephys_in_window;
        t_bump_w    = cur_data.t_bump_w;
        t_yaw_w     = cur_data.t_yaw_w;
        t_ephys_w   = cur_data.t_ephys_w;     
        VPS         = cur_data.VPS;
        
        dt_bump = t_bump_w(2) - t_bump_w(1);

        bump_vel_win = [];
        bump_pos_win = [];
        yaw_win      = [];
        fwd_win      = [];
        Vm_win       = [];
        PSTH_win     = [];
        t_bump_win   = [];
        t_yaw_win    = [];
        t_Vm_win     = [];
        
        ephys_win              = [];
        bump_return_peak_t     = [];
        t_peak     = t_bump_w;
                        
        if( DEBUG_LEVEL == DEBUG_VERBOSE )
            aligned_data_fig = figure('units','normalized','outerposition',[0 0 1 1]);
        end
        
        kk = 0;
        ccc = 0;
        stims_passed_all_checks = [];
        
        for s = 1:size(cur_cond_bump_tc,1)
            
            % Calculate bump velocity
            cur_bump_tmp = cur_cond_bump_tc(s,:);

            cur_bump_tmp_FF = cur_cond_bump_tc_FF(s,:);
            
            % if bump is a nan for large enough consequitive time points,
            % than disqualify.
            [ cur_bump ] = assess_and_fix_bump_quality( cur_bump_tmp, VPS );                       
            
            if(length(cur_bump) == 0 )
                disp(['Trial: ' num2str(s) ' has been disqualified, too many time points with no clear bump.']);
                continue;
            end
            
            [ cur_bump_FF ] = assess_and_fix_bump_quality( cur_bump_tmp_FF, VPS );
            
            if( ( BUMP_CONDITION == BUMP_CONDITION_RETURNED_UP ) || ( BUMP_CONDITION == BUMP_CONDITION_RETURNED_DOWN ) )
                cur_bump_return_idx_range = cur_cond_bump_return_idx{ s };
                
                assert( cur_bump_return_idx_range( end ) > -1 );
                
                cur_bump_tc = cur_bump( cur_bump_return_idx_range );
            else
                cur_bump_tc = cur_bump;
                cur_bump_return_idx_range = [1:length(cur_bump_tc)];
            end
                       
            cur_stim_id = cur_cond_stim_ids(s);
            
            cur_bump_rois  = squeeze(bump_data( cur_stim_id, :, : ));
            cur_fwd   = fwd_data( cur_stim_id, : );
            
            % cur_yaw_drift   = yaw_data( cur_stim_id, : );
            cur_yaw_drift   = convert_yaw_to_degrees( yaw_data( cur_stim_id, : ) );
            cur_yaw         = cur_yaw_drift - mean( cur_yaw_drift );
            
            cur_ephys = ephys_data( cur_stim_id, : );
            
            % remove spikes and drift
            FILT_FACTOR = 0.04;
            cur_Vm_w_drift = medfilt1( cur_ephys, FILT_FACTOR * EPHYS_FR, 'truncate' );
            cur_Vm = cur_Vm_w_drift - mean( cur_Vm_w_drift );                      
            
            % With PSTH 
            SPIKE_THRESHOLD_A2 = directories{ d }{ 3 };
            
            cur_debug = 0;
             if( s == 4 )
                 cur_debug = 1;
             end
             
            [ cur_PSTH ] = calculate_psth_debug( t_ephys_w(1:end-1), t_yaw_w(1:end-1), cur_ephys(1:end-1), EPHYS_FR, SPIKE_THRESHOLD_A2, psth_dt_samples, cur_debug );
            
            % Take max or min depending on condition
            bump_vel = diff( cur_bump_tc ) / dt_bump;
            
            bump_vel_all = diff( cur_bump ) / dt_bump;
            
            if( BUMP_CONDITION == BUMP_CONDITION_RETURNED_UP )
                bump_vel_to_search = bump_vel;
            elseif( BUMP_CONDITION == BUMP_CONDITION_RETURNED_DOWN )
                bump_vel_to_search = -1.0*bump_vel;            
            else
                bump_vel_to_search = bump_vel;                
            end
            
            % Remove artifactual peaks from bump wrapping
            bump_vel_to_search( find( abs(bump_vel_to_search) >= BUMP_SPEED_THRESHOLD_MAX ) ) = 0;
            
                                    
            % Check that the average bump speed is about a threshold. This
            % eliminates cases where there isn't a clear peak in bump
            % movement.
            if( mean(abs(bump_vel_to_search)) < BUMP_SPEED_THRESHOLD_MIN )
                continue;
            end
            
            % Filter out noise from            
            [~, locs] = findpeaks( bump_vel_to_search, 'NPeaks', 1, 'SortStr', 'descend' ); 
            
            if( length(locs) == 0 )
                % No peaks detected
                continue;
            end
                
            EB_bump_vel_align_idx_in_bump_vel_to_search = locs(1);
            
            if( ( BUMP_CONDITION == BUMP_CONDITION_RETURNED_UP ) || ( BUMP_CONDITION == BUMP_CONDITION_RETURNED_DOWN ) )
                EB_bump_vel_align_idx = locs(1) + cur_bump_return_idx_range(1);
            else
                EB_bump_vel_align_idx = locs(1);
            end
            
            if( DEBUG_LEVEL == DEBUG_VERBOSE )            
            %if( s == 17 ) % 17 is a nice example for fly _17
            % if( 0 )            
                f = figure('units','normalized','outerposition',[0 0 1 1]);

                ax1(1) = subplot(5, 1, 1);
                hold on;
                imagesc( t_bump_w, [1:size(cur_bump_rois,1)], cur_bump_rois );
                colormap(flipud(gray));
                axis tight;
                ylim([0 9]);
                caxis([-0.5 2]);

                tt = title( [cur_cond_str ' stim: ' num2str(s)  ] );
                set(tt, 'Interpreter', 'none');

                if( s == 12 )
                    ccc = ccc + 1;
                end
                
                BUMP_OFFSET = 6.0;
                plot( t_bump_w, cur_bump + BUMP_OFFSET );

                ax1(2) = subplot(5, 1, 2);
                hold on;
                plot( t_bump_w( cur_bump_return_idx_range(1:end-1) ), bump_vel_to_search );
                plot( t_bump_w( cur_bump_return_idx_range(1)), 0, 'o' );
                plot( t_bump_w( cur_bump_return_idx_range(end-1)), 0, 'o' );
                plot( t_bump_w( EB_bump_vel_align_idx-1 ), bump_vel_to_search( EB_bump_vel_align_idx_in_bump_vel_to_search ), 'X', 'color', 'g' );
                
                ylabel('EB bump vel (au/s)');
                
                ax1(3) = subplot(5, 1, 3);
                hold on;
                plot( t_yaw_w, cur_yaw );                
                ylabel('Yaw vel (au/s)');
                
                ax1(4) = subplot(5, 1, 4);
                hold on;
                plot( t_ephys_w, cur_ephys -    mean( cur_Vm_w_drift ) );                
                ylabel('Vm (mV)');
                xlabel('Time (s)');

                ax1(5) = subplot(5, 1, 5);
                hold on;
                plot( t_yaw_w(1:end-1), cur_PSTH );                
                ylabel('Firing rate (spikes/s)');
                xlabel('Time (s)');

                % waitforbuttonpress;
                linkaxes(ax1, 'x');
                saveas(f, [t_now_analysis_path '/eb_bump_vel_peak_t_detect_' cur_cond_str '_stim_' num2str( s ) '.fig'] );
                saveas(f, [t_now_analysis_path '/eb_bump_vel_peak_t_detect_' cur_cond_str '_stim_' num2str( s ) '.png'] );
                close(f);
            end
            
            % Align peak of bump
            EB_bump_vel_align = t_bump_w( EB_bump_vel_align_idx-1 );           
            
            EB_FRAMES_BEFORE_EB_VEL_CHANGE = floor( TIME_BEFORE_EB_VEL_CHANGE * VPS );
            EB_FRAMES_AFTER_EB_VEL_CHANGE  = floor( TIME_AFTER_EB_VEL_CHANGE * VPS );
            
            cur_EB_bump_vel_win_start  = EB_bump_vel_align_idx - EB_FRAMES_BEFORE_EB_VEL_CHANGE;
            cur_EB_bump_vel_win_end    = EB_bump_vel_align_idx + EB_FRAMES_AFTER_EB_VEL_CHANGE;
            
            if( ( cur_EB_bump_vel_win_start < 1) || ( cur_EB_bump_vel_win_end >= length( cur_bump ) ) )
                % skip this stimulus
                disp(['Stim removed: out of bounds in alignment: range: { 1, ' num2str(length( cur_bump )) ' }  cur stim: { ' num2str(cur_EB_bump_vel_win_start) ' , ' num2str(cur_EB_bump_vel_win_end) ' }']);
                continue;
            end
            
%             if( cur_EB_bump_vel_win_start < 1) 
%                 cur_EB_bump_vel_win_start = 1;
%             end            
%             
%             if( cur_EB_bump_vel_win_end >= length( cur_bump ) )
%                 cur_EB_bump_vel_win_end = length( cur_bump )-1; % Need to account for n-1 in vel vs. pos
%             end
                  
            % Align yaw
            xx = find( t_yaw_w < EB_bump_vel_align );
            yaw_align_idx = xx(end);
            
            YAW_FRAMES_BEFORE_EB_VEL_CHANGE = floor( TIME_BEFORE_EB_VEL_CHANGE * BALL_FR );
            YAW_FRAMES_AFTER_EB_VEL_CHANGE  = floor( TIME_AFTER_EB_VEL_CHANGE * BALL_FR );
            cur_yaw_win_start  = yaw_align_idx - YAW_FRAMES_BEFORE_EB_VEL_CHANGE;
            cur_yaw_win_end    = yaw_align_idx + YAW_FRAMES_AFTER_EB_VEL_CHANGE;
            
            if( ( cur_yaw_win_start < 1) || ( cur_yaw_win_end > length( cur_yaw ) ) )
                % skip this stimulus
                disp(['Stim removed: out of bounds in alignment: range: { 1, ' num2str(length( cur_yaw )) ' }  cur stim: { ' num2str(cur_yaw_win_start) ' , ' num2str(cur_yaw_win_end) ' }']);                
                continue;
            end
            
%             if( cur_yaw_win_start < 1)
%                 cur_yaw_win_start = 1;
%             end            
%             
%             if( cur_yaw_win_end > length( cur_yaw )  )
%                 cur_yaw_win_end = length( cur_yaw );
%             end
            
            % Align ephys (Vm)
            xx = find( t_ephys_w < EB_bump_vel_align );
            ephys_align_idx = xx(end);
            
            EPHYS_FRAMES_BEFORE_EB_VEL_CHANGE = floor( TIME_BEFORE_EB_VEL_CHANGE * EPHYS_FR );
            EPHYS_FRAMES_AFTER_EB_VEL_CHANGE  = floor( TIME_AFTER_EB_VEL_CHANGE * EPHYS_FR );
            cur_Vm_win_start  = ephys_align_idx - EPHYS_FRAMES_BEFORE_EB_VEL_CHANGE;
            cur_Vm_win_end    = ephys_align_idx + EPHYS_FRAMES_AFTER_EB_VEL_CHANGE;

            if( ( cur_Vm_win_start < 1) || ( cur_Vm_win_end > length( cur_ephys ) ) )
                % skip this stimulus
                disp(['Stim removed: out of bounds in alignment: range: { 1, ' num2str(length( cur_ephys )) ' }  cur stim: { ' num2str(cur_Vm_win_start) ' , ' num2str(cur_Vm_win_end) ' }']);
                continue;
            end
            
%             if( cur_Vm_win_start < 1)
%                 cur_Vm_win_start = 1;
%             end            
%             
%             if( cur_Vm_win_end > length( cur_ephys )  )
%                 cur_Vm_win_end = length( cur_ephys );
%             end            
            
            % If we got this far, then the indicies are within range
            cur_bump_win = cur_bump( cur_EB_bump_vel_win_start:cur_EB_bump_vel_win_end );
            cur_bump_pre_win = cur_bump( 1:cur_EB_bump_vel_win_start ); 
            
            cur_bump_win_FF = cur_bump_FF( cur_EB_bump_vel_win_start:cur_EB_bump_vel_win_end );
            cur_bump_pre_win_FF = cur_bump_FF( 1:cur_EB_bump_vel_win_start );             
            
            bump_pos_win(end+1,:) = cur_bump_win;            
            
            cur_EB_vel_win = bump_vel_all( cur_EB_bump_vel_win_start:cur_EB_bump_vel_win_end );
            bump_vel_win(end+1,:) = cur_EB_vel_win;
            t_bump_win = t_bump_w( cur_EB_bump_vel_win_start:cur_EB_bump_vel_win_end ) - EB_bump_vel_align;
            t_bump_pre_win = t_bump_w( 1:cur_EB_bump_vel_win_start ) - EB_bump_vel_align;
            
            % Align yaw
            cur_yaw_win = cur_yaw( cur_yaw_win_start:cur_yaw_win_end );
            yaw_win(end+1,:) = cur_yaw_win;
            t_yaw_win = t_yaw_w( cur_yaw_win_start:cur_yaw_win_end ) - t_yaw_w(yaw_align_idx);

            % Align fwd
            cur_fwd_win = cur_fwd( cur_yaw_win_start:cur_yaw_win_end );
            fwd_win(end+1,:) = cur_fwd_win;            
            
            % Align ephys (Vm)
            cur_Vm_win = cur_Vm( cur_Vm_win_start:cur_Vm_win_end );
            Vm_win(end+1,:) = cur_Vm_win;
            t_Vm_win = t_ephys_w( cur_Vm_win_start:cur_Vm_win_end ) - t_ephys_w(ephys_align_idx);
            ephys_win(end+1, :) = cur_ephys( cur_Vm_win_start:cur_Vm_win_end );

            % Align ephys (FR)
            cur_PSTH_win = cur_PSTH( cur_yaw_win_start:cur_yaw_win_end );
            PSTH_win(end+1,:) = cur_PSTH_win;
            
            stims_passed_all_checks(end+1) = s;
            
            if( DEBUG_LEVEL == DEBUG_VERBOSE )
                figure(aligned_data_fig);
                ax2(1) = subplot(6,1,1);
                hold on;
                plot( t_bump_win, cur_bump_win, '-', 'DisplayName', [ 'stim: ' num2str(s) ] );
                ylim([-5 5]);
                ylabel('PVA position (wedge loc)');
                tt = title( [cur_datapath ': ' cur_cond_str ]);
                set(tt, 'Interpreter', 'none');

                ax2(2) = subplot(6,1,2);
                hold on;
                plot( t_bump_win, cur_EB_vel_win, '-' );
                ylim([-15 15]);
                ylabel('EB Vel (au/s)');
                
                ax2(3) = subplot(6,1,3);
                hold on;
                plot( t_yaw_win, cur_fwd_win, '-' );
                ylabel('Fwd (au/s)');
                
                ax2(4) = subplot(6,1,4);
                hold on;
                plot( t_yaw_win, cur_yaw_win, '-' );
                ylabel('Yaw (au/s)');
                
                ax2(5) = subplot(6,1,5);
                hold on;
                plot( t_Vm_win, cur_Vm_win, '-' );
                ylabel('Vm (mV)');
                
                ax2(6) = subplot(6,1,6);
                hold on;
                plot( t_yaw_win, cur_PSTH_win, '-' );
                ylabel('PSTH (spikes/s)');
                xlabel('Time (s)');
                xlim([ t_yaw_win(1) t_yaw_win(end) ]);
                %linkaxes(ax1, 'x');
                
                figure( bump_and_yaw_fig );
                yyaxis left;
                hold on;
                plot( t_bump_win, cur_EB_vel_win, '-b' );
                set(gca, 'TickDir', 'Out');
                ylabel('EB bump vel (wed/s)');
                
                yyaxis right;
                plot( t_yaw_win, cur_yaw_win, '-r' );
                set(gca, 'TickDir', 'Out');
                xlabel('Time (s)');
                ylabel('vyaw (deg/s)');

                
                %%%%%%%%%
                figure( bump_trial_to_trial_colored_by_yaw_fig );

                pre_bump_max_vel_win_t = find( (t_yaw_win <= 0) & (t_yaw_win > -0.5 ) );
                
                avg_yaw = mean( cur_yaw_win( pre_bump_max_vel_win_t ) );
                
                cur_clr = 'none';
                if( avg_yaw > 0 ) cur_clr = 'b'; else cur_clr = 'g'; end
                                    
                subplot( 2, 1, 1 );
                hold on;
                plot( t_bump_win, cur_bump_win_FF, 'color', cur_clr );
                plot( t_bump_pre_win, cur_bump_pre_win_FF, 'color', cur_clr );
                set(gca, 'TickDir', 'Out');
                ylabel('EB bump position (wed)');
                
                subplot( 2, 1, 2 );
                hold on;        
                plot( t_bump_w, cur_bump_FF, 'color', cur_clr );
                set(gca, 'TickDir', 'Out');
                ylabel('EB bump position (wed)');
                
%                 yyaxis right;
%                 plot( t_yaw_win, cur_yaw_win, '-r' );
%                 set(gca, 'TickDir', 'Out');
%                 xlabel('Time (s)');
%                 ylabel('vyaw (deg/s)');
                               
            end
        end
        
        bump_pos_win_all{ cond, d }       = bump_pos_win;
        bump_win_all{ cond, d }           = bump_vel_win;
        yaw_win_all{ cond, d }            = yaw_win;
        fwd_win_all{ cond, d }            = fwd_win;
        Vm_win_all{ cond, d }             = Vm_win;
        PSTH_win_all{ cond, d }           = PSTH_win;
        timebase_bump{ cond, d }          = t_bump_win;
        timebase_yaw{ cond, d }           = t_yaw_win;
        timebase_ephys{ cond, d }         = t_Vm_win;
        ephys_win_all{ cond, d }          = ephys_win;
        final_stims_passed_all_checks{ cond, d } = stims_passed_all_checks;
        
        if( DEBUG_LEVEL == DEBUG_VERBOSE )
            
            if( length( bump_vel_win ) == 0 )
                continue;
            end
            
            figure(aligned_data_fig);
            ax2(1) = subplot(6,1,1);
            hold on;
            
            if( size(bump_pos_win, 1) == 1 )
                avg_bump_pos = bump_pos_win;
            else
                avg_bump_pos = mean( bump_pos_win );
            end
            
            plot( t_bump_win, avg_bump_pos, '-', 'LineWidth', 2.0 );
            ylim([ -5 5 ]);
            ylabel('PVA position (wedge loc)');
            tt = title( [cur_datapath ': ' cur_cond_str ' ( ' num2str(size( bump_vel_win, 1 )) ')']);
            set(tt, 'Interpreter', 'none');

            ax2(2) = subplot(6,1,2);
            hold on;
            
            if( size(bump_vel_win, 1) == 1 )
                avg_bump = bump_vel_win;
            else
                avg_bump = mean( bump_vel_win );
            end
            
            plot( t_bump_win, avg_bump, '-', 'LineWidth', 2.0 );
            ylim([-15 15]);
            ylabel('EB Vel (au/s)');
    
            ax2(3) = subplot(6,1,3);
            hold on;
            
            if( size(fwd_win, 1) == 1 )
                avg_fwd = fwd_win;
            else
                avg_fwd = mean( fwd_win );
            end
            
            plot( t_yaw_win, avg_fwd, '-', 'LineWidth', 2.0 );
            ylabel('Fwd (au/s)');
            
            ax2(4) = subplot(6,1,4);
            hold on;
            if( size(yaw_win, 1) == 1 )
                avg_yaw = yaw_win;
            else
                avg_yaw = mean( yaw_win );
            end
            
            plot( t_yaw_win, avg_yaw, '-', 'LineWidth', 2.0 );
            ylabel('Yaw (au/s)');
            
            ax2(5) = subplot(6,1,5);
            hold on;
            
            if( size(Vm_win, 1) == 1 )
                avg_ephys = Vm_win;
            else
                avg_ephys = mean( Vm_win );
            end                
            
            plot( t_Vm_win, avg_ephys, '-', 'LineWidth', 2.0 );
            ylabel('Vm (mV)');

            ax2(6) = subplot(6,1,6);
            hold on;
            
            if( size(PSTH_win, 1) == 1 )
                avg_PSTH = PSTH_win;
            else
                avg_PSTH = mean( PSTH_win );
            end                
            
            plot( t_yaw_win, avg_PSTH, '-', 'LineWidth', 2.0 );
            ylabel('Firing rate (spike/s)');
            
            xlabel('Time (s)');
            xlim([ t_yaw_win(1) t_yaw_win(end) ]);
            linkaxes(ax2, 'x');
            
            saveas( aligned_data_fig, [analysis_path '/' cur_cond_str '_bump_return_yaw_fwd_ephys_aligned.fig']);
            saveas( aligned_data_fig, [analysis_path '/' cur_cond_str '_bump_return_yaw_fwd_ephys_aligned.png']);
            
            
            figure( bump_and_yaw_fig );
            yyaxis left;
            hold on;
            plot( t_bump_win, avg_bump, 'b', 'LineWidth', 2.5 );
            set(gca, 'TickDir', 'Out');
            ylabel('EB bump vel (wed/s)');
            ylim([-10 10]);
            
            yyaxis right;
            plot( t_yaw_win, avg_yaw, 'r', 'LineWidth', 2.5 );
            set(gca, 'TickDir', 'Out');
            xlabel('Time (s)');
            ylabel('vyaw (deg/s)');
            ylim([-800 800]);
            xlim([-1.0 1.0]);
            
            title(['Number of trials: ' num2str(size(yaw_win,1)) ]);
                        
            saveas( bump_and_yaw_fig, [analysis_path '/bump_vel_vs_yaw_trial_by_trial_for_CX_figure.fig'] );
            saveas( bump_and_yaw_fig, [analysis_path '/bump_vel_vs_yaw_trial_by_trial_for_CX_figure.png'] );
            saveas( bump_and_yaw_fig, [analysis_path '/bump_vel_vs_yaw_trial_by_trial_for_CX_figure.svg'] );
        end
    end    
end

saveas(bump_trial_to_trial_colored_by_yaw_fig, [analysis_path 'bump_vel_colorcoded_by_yaw.fig']);
saveas(bump_trial_to_trial_colored_by_yaw_fig, [analysis_path 'bump_vel_colorcoded_by_yaw.png']);
saveas(bump_trial_to_trial_colored_by_yaw_fig, [analysis_path 'bump_vel_colorcoded_by_yaw.svg']);


%
% bump_win_all   = cell( length(bump_conditions), length( directories ));
% yaw_win_all    = cell( length(bump_conditions), length( directories ));
% fwd_win_all    = cell( length(bump_conditions), length( directories ));
% Vm_win_all  = cell( length(bump_conditions), length( directories ));
% timebase_bump  = cell( length(bump_conditions), length( directories ));
% timebase_yaw   = cell( length(bump_conditions), length( directories ));
% timebase_ephys = cell( length(bump_conditions), length( directories ));

end

