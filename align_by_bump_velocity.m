function [ bump_win_all, yaw_win_all, fwd_win_all, ephys_win_all, timebase_bump, timebase_yaw, timebase_ephys ] = align_by_bump_velocity( basedir, directories, bump_conditions, bump_conditions_str )

settings = sensor_settings;
BALL_FR = settings.sensorPollFreq;
EPHYS_FR = settings.sampRate;

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
DEBUG_LEVEL                      = DEBUG_OFF;

TIME_BEFORE_EB_VEL_CHANGE        = 1.0;
TIME_AFTER_EB_VEL_CHANGE         = 2.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bump_win_all   = cell( length(bump_conditions), length( directories ));
yaw_win_all    = cell( length(bump_conditions), length( directories ));
fwd_win_all    = cell( length(bump_conditions), length( directories ));
ephys_win_all  = cell( length(bump_conditions), length( directories ));
timebase_bump  = cell( length(bump_conditions), length( directories ));
timebase_yaw   = cell( length(bump_conditions), length( directories ));
timebase_ephys = cell( length(bump_conditions), length( directories ));

for cond = 1:length(bump_conditions)
    cur_cond     = bump_conditions{cond};
    cur_cond_str = bump_conditions_str{cond};
    
    BUMP_CONDITION = BUMP_CONDITION_UNDEFINED;
    if( strcmp(cur_cond_str, 'bump_jumps_up_returns_down') == 1 ) 
        BUMP_CONDITION = BUMP_CONDITION_RETURNED_DOWN;
    elseif( strcmp(cur_cond_str, 'bump_jumps_down_returns_up') == 1 )  
        BUMP_CONDITION = BUMP_CONDITION_RETURNED_UP;
    elseif( strcmp(cur_cond_str, 'no_response') == 1 )  
        BUMP_CONDITION = BUMP_CONDITION_NO_RESPONSE;
    elseif( strcmp(cur_cond_str, 'bump_not_returned_up') == 1 )  
        BUMP_CONDITION = BUMP_CONDITION_NOT_RETURNED_UP;
    elseif( strcmp(cur_cond_str, 'bump_not_returned_down') == 1 )  
        BUMP_CONDITION = BUMP_CONDITION_NOT_RETURNED_DOWN;
    end
        
    for d = 1:length( directories )
    % for d = 9
        
        cur_datapath = directories{ d }{ 1 };
        cur_sid      = directories{ d }{ 2 };
        
        datapath = [ basedir '/' cur_datapath ];
        analysis_path = [datapath '/analysis/'];

        t_now = datetime('now');
        t_now_analysis_path = [datapath '/analysis/' datestr(t_now)];
        mkdir(t_now_analysis_path);
        
        cur_cond_bump_tc  = cur_cond{d,1};
        cur_cond_stim_ids = cur_cond{d,2};
              
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

        bump_win    = [];
        yaw_win     = [];
        fwd_win     = [];
        ephys_win   = [];
        t_bump_win  = [];
        t_yaw_win   = [];
        t_ephys_win = [];
        
        if( DEBUG_LEVEL == DEBUG_VERBOSE )
            aligned_data_fig = figure('units','normalized','outerposition',[0 0 1 1]);
        end
        
        kk = 0;
        ccc = 0;
        for s = 1:size(cur_cond_bump_tc,1)
            
            % Calculate bump velocity
            cur_bump_tmp = cur_cond_bump_tc(s,:);
            
            % if bump is a nan for large enough consequitive time points,
            % than disqualify.
            [ cur_bump ] = assess_and_fix_bump_quality( cur_bump_tmp, VPS );
            
            if(length(cur_bump) == 0 )
                disp(['Trial: ' num2str(s) ' has been disqualified, too many time points with no clear bump.']);
                continue;
            end

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
            
            cur_yaw_drift   = yaw_data( cur_stim_id, : );
            cur_yaw         = cur_yaw_drift - mean( cur_yaw_drift );
            
            cur_ephys = ephys_data( cur_stim_id, : );
            
            % remove spikes and drift
            FILT_FACTOR = 0.06;
            cur_Vm_w_drift = medfilt1( cur_ephys, FILT_FACTOR * EPHYS_FR, 'truncate' );
            cur_Vm = cur_Vm_w_drift - mean( cur_Vm_w_drift );
            
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
            % if( 0 )            
                f = figure('units','normalized','outerposition',[0 0 1 1]);

                ax1(1) = subplot(4, 1, 1);
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

                ax1(2) = subplot(4, 1, 2);
                hold on;
                plot( t_bump_w( cur_bump_return_idx_range(1:end-1) ), bump_vel_to_search );
                plot( t_bump_w( cur_bump_return_idx_range(1)), 0, 'o' );
                plot( t_bump_w( cur_bump_return_idx_range(end-1)), 0, 'o' );
                plot( t_bump_w( EB_bump_vel_align_idx-1 ), bump_vel_to_search( EB_bump_vel_align_idx_in_bump_vel_to_search ), 'X', 'color', 'g' );
                
                ylabel('EB bump vel (au/s)');
                
                ax1(3) = subplot(4, 1, 3);
                hold on;
                plot( t_yaw_w, cur_yaw );                
                ylabel('Yaw vel (au/s)');
                
                ax1(4) = subplot(4, 1, 4);
                hold on;
                plot( t_ephys_w, cur_Vm );                
                ylabel('Vm (mV)');
                xlabel('Time (s)');
                
                
                % waitforbuttonpress;
                linkaxes(ax1, 'x');
                saveas(f, [t_now_analysis_path '/eb_bump_vel_peak_detect_' cur_cond_str '_stim_' num2str( s ) '.fig'] );
                saveas(f, [t_now_analysis_path '/eb_bump_vel_peak_detect_' cur_cond_str '_stim_' num2str( s ) '.png'] );
                close(f);
            end
            
            % Align peak of bump
            EB_bump_vel_align = t_bump_w( EB_bump_vel_align_idx-1 );           
            
            EB_FRAMES_BEFORE_EB_VEL_CHANGE = floor( TIME_BEFORE_EB_VEL_CHANGE * VPS );
            EB_FRAMES_AFTER_EB_VEL_CHANGE  = floor( TIME_AFTER_EB_VEL_CHANGE * VPS );
            
            cur_EB_bump_vel_win_start  = EB_bump_vel_align_idx - EB_FRAMES_BEFORE_EB_VEL_CHANGE;
            cur_EB_bump_vel_win_end  = EB_bump_vel_align_idx + EB_FRAMES_AFTER_EB_VEL_CHANGE;
            
            if( ( cur_EB_bump_vel_win_start < 1) || ( cur_EB_bump_vel_win_end >= length( cur_bump ) ) )
                % skip this stimulus
                continue;
            end            
                  
            % Align yaw
            xx = find( t_yaw_w < EB_bump_vel_align );
            yaw_align_idx = xx(end);
            
            YAW_FRAMES_BEFORE_EB_VEL_CHANGE = floor( TIME_BEFORE_EB_VEL_CHANGE * BALL_FR );
            YAW_FRAMES_AFTER_EB_VEL_CHANGE  = floor( TIME_AFTER_EB_VEL_CHANGE * BALL_FR );
            cur_yaw_win_start  = yaw_align_idx - YAW_FRAMES_BEFORE_EB_VEL_CHANGE;
            cur_yaw_win_end  = yaw_align_idx + YAW_FRAMES_AFTER_EB_VEL_CHANGE;
            
            if( ( cur_yaw_win_start < 1) || ( cur_yaw_win_end > length( cur_yaw ) ) )
                % skip this stimulus
                continue;
            end            
            
            % Align ephys
            xx = find( t_ephys_w < EB_bump_vel_align );
            ephys_align_idx = xx(end);
            
            EPHYS_FRAMES_BEFORE_EB_VEL_CHANGE = floor( TIME_BEFORE_EB_VEL_CHANGE * EPHYS_FR );
            EPHYS_FRAMES_AFTER_EB_VEL_CHANGE  = floor( TIME_AFTER_EB_VEL_CHANGE * EPHYS_FR );
            cur_ephys_win_start  = ephys_align_idx - EPHYS_FRAMES_BEFORE_EB_VEL_CHANGE;
            cur_ephys_win_end    = ephys_align_idx + EPHYS_FRAMES_AFTER_EB_VEL_CHANGE;

            if( ( cur_ephys_win_start < 1) || ( cur_ephys_win_end > length( cur_ephys ) ) )
                % skip this stimulus
                continue;
            end                        
            
            % If we got this far, then the indicies are within range
            cur_EB_vel_win = bump_vel_all( cur_EB_bump_vel_win_start:cur_EB_bump_vel_win_end );
            bump_win(end+1,:) = cur_EB_vel_win;
            t_bump_win = t_bump_w( cur_EB_bump_vel_win_start:cur_EB_bump_vel_win_end ) - EB_bump_vel_align;

            % Align yaw
            cur_yaw_win = cur_yaw( cur_yaw_win_start:cur_yaw_win_end );
            yaw_win(end+1,:) = cur_yaw_win;
            t_yaw_win = t_yaw_w( cur_yaw_win_start:cur_yaw_win_end ) - t_yaw_w(yaw_align_idx);

            % Align fwd
            cur_fwd_win = cur_fwd( cur_yaw_win_start:cur_yaw_win_end );
            fwd_win(end+1,:) = cur_fwd_win;            
            
            % Align ephys
            cur_ephys_win = cur_Vm( cur_ephys_win_start:cur_ephys_win_end );
            ephys_win(end+1,:) = cur_ephys_win;
            t_ephys_win = t_ephys_w( cur_ephys_win_start:cur_ephys_win_end ) - t_ephys_w(ephys_align_idx);
 
            if( DEBUG_LEVEL == DEBUG_VERBOSE )
                figure(aligned_data_fig);
                ax1(1) = subplot(4,1,1);
                hold on;
                plot( t_bump_win, cur_EB_vel_win, '-' );
                ylim([-15 15]);
                ylabel('EB Vel (au/s)');
                tt = title( [cur_datapath ': ' cur_cond_str ]);
                set(tt, 'Interpreter', 'none');
                
                ax1(2) = subplot(4,1,2);
                hold on;
                plot( t_yaw_win, cur_fwd_win, '-' );
                ylabel('Fwd (au/s)');
                xlabel('Time (s)');
                
                ax1(3) = subplot(4,1,3);
                hold on;
                plot( t_yaw_win, cur_yaw_win, '-' );
                ylabel('Yaw (au/s)');
                
                ax1(4) = subplot(4,1,4);
                hold on;
                plot( t_ephys_win, cur_ephys_win, '-' );
                ylabel('Vm (mV)');
                xlabel('Time (s)');
            end
        end
        
        bump_win_all{ cond, d }   = bump_win;
        yaw_win_all{ cond, d }    = yaw_win;
        fwd_win_all{ cond, d }    = fwd_win;
        ephys_win_all{ cond, d }  = ephys_win;
        timebase_bump{ cond, d }  = t_bump_win;
        timebase_yaw{ cond, d }   = t_yaw_win;
        timebase_ephys{ cond, d } = t_ephys_win;
        
        if( DEBUG_LEVEL == DEBUG_VERBOSE )
            
            if( length( bump_win ) == 0 )
                continue;
            end
            
            figure(aligned_data_fig);
            ax1(1) = subplot(4,1,1);
            hold on;
            
            if( size(bump_win, 1) == 1 )
                avg_bump = bump_win;
            else
                avg_bump = mean( bump_win );
            end
            
            plot( t_bump_win, avg_bump, '-', 'LineWidth', 2.0 );
            ylim([-15 15]);
            ylabel('EB Vel (au/s)');
            tt = title( [cur_datapath ': ' cur_cond_str ' ( ' num2str(size( bump_win, 1 )) ')']);
            set(tt, 'Interpreter', 'none');
            
            ax1(2) = subplot(4,1,2);
            hold on;
            
            if( size(fwd_win, 1) == 1 )
                avg_fwd = fwd_win;
            else
                avg_fwd = mean( fwd_win );
            end
            
            plot( t_yaw_win, avg_fwd, '-', 'LineWidth', 2.0 );
            ylabel('Fwd (au/s)');
            
            ax1(3) = subplot(4,1,3);
            hold on;
            if( size(yaw_win, 1) == 1 )
                avg_yaw = yaw_win;
            else
                avg_yaw = mean( yaw_win );
            end
            
            plot( t_yaw_win, avg_yaw, '-', 'LineWidth', 2.0 );
            ylabel('Yaw (au/s)');
            
            ax1(4) = subplot(4,1,4);
            hold on;
            
            if( size(ephys_win, 1) == 1 )
                avg_ephys = ephys_win;
            else
                avg_ephys = mean( ephys_win );
            end                
            
            plot( t_ephys_win, avg_ephys, '-', 'LineWidth', 2.0 );
            ylabel('Vm (mV)');
            xlabel('Time (s)');
            linkaxes(ax1, 'x');
            
            saveas( aligned_data_fig, [analysis_path '/' cur_cond_str '_bump_return_yaw_fwd_ephys_aligned.fig']);
            saveas( aligned_data_fig, [analysis_path '/' cur_cond_str '_bump_return_yaw_fwd_ephys_aligned.png']);
        end
    end    
end
end

