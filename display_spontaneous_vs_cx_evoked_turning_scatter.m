function display_spontaneous_vs_cx_evoked_turning_scatter( basedir, cur_dirs, yaw_win_all, ephys_win_all, timebase_yaw, timebase_ephys )

FILT_FACTOR = 0.04;
ac = get_analysis_constants;

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;
BIN_SIZE = 0.050; % s
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

CX_TURN_WINDOW_T = 0.4;

for d = 1:length( cur_dirs )
% for d = 1
    
    f = figure;
    
    for cond = 1:size( ephys_win_all,1 )
                
        cur_datapath = cur_dirs{ d }{ 1 };
        cur_sid      = cur_dirs{ d }{ 2 };
        
        datapath = [ basedir '/' cur_datapath ];
        analysis_path = [datapath '/analysis/'];
        
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
            
            VmFilt_A2 = medfilt1( cur_ephys, FILT_FACTOR * ephys_SR, 'truncate' );
            VmFilt_A2_corr = VmFilt_A2 - mean(VmFilt_A2);
            
            t_down = squeeze(mean(reshape(cur_ephys_t, [DT_EPHYS, length(cur_ephys_t)/DT_EPHYS]), 1));
            A2_Vm_down = squeeze(mean(reshape(VmFilt_A2_corr, [ DT_EPHYS, length(VmFilt_A2)/DT_EPHYS ] ),1));
            
            yaw_t_down = squeeze( mean( reshape( cur_yaw_t, [ DT_YAW, length( cur_yaw_t ) / DT_YAW ] ),1));
            yaw_down = squeeze( mean( reshape( cur_yaw, [ DT_YAW, length( cur_yaw_t ) / DT_YAW ] ),1));
            
            hold on;
            
            % Plots spontaneous data
            SHIFT_FACTOR = 3;
            for ii = 1:(length(A2_Vm_down)-SHIFT_FACTOR)
                cur_yaw_1 = yaw_down( ii+SHIFT_FACTOR );
                plot(A2_Vm_down(ii), cur_yaw_1, 'o', 'MarkerSize', 3, 'color', 'b' );
            end
            
            % Assuming the data is aligned to bump return velocity.
            time_of_peak = 0.0;                       
            CX_turn_start = time_of_peak - CX_TURN_WINDOW_T;
            CX_turn_end   = time_of_peak + CX_TURN_WINDOW_T;                        
                        
            % Plot CX-evoked turning epoch            
            cur_t_down     = find( ( t_down >= CX_turn_start ) & ( t_down <= CX_turn_end ) );
            cur_yaw_t_down = find( ( yaw_t_down >= CX_turn_start ) & ( yaw_t_down <= CX_turn_end ) );
                                    
            for ii = 1:(length(cur_t_down)-SHIFT_FACTOR)
                cur_index_yaw   = cur_yaw_t_down( ii );
                cur_index_ephys = cur_t_down( ii );
                cur_yaw_2 = yaw_down( cur_index_yaw +SHIFT_FACTOR );
                plot(A2_Vm_down(cur_index_ephys), cur_yaw_2, 'o', 'MarkerSize', 10, 'color', 'r' );
            end            
        end
    end
    
    tt = title(['Fly: ' num2str( d ) ' ' cur_datapath ]);
    set(tt, 'Interpreter', 'none');

    xlabel('A2 Vm (mV)');
    ylabel('Yaw (deg/s)');

    saveas(f,[analysis_path '/Vm_vs_yaw_spontaneous_CX_turning_Vm.fig']);
    saveas(f,[analysis_path '/Vm_vs_yaw_spontaneous_CX_turning_Vm.png']);
    % close(f);
end


% display linear regression slope of spontaneous vs. CX-evoked (connecting
% lines within a fly)

end

