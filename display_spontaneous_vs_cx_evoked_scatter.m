function display_spontaneous_vs_cx_evoked_turning_scatter(  cur_dirs, yaw_win_all, Vm_win_all, bump_return_peak_t_all, t_peak_t )

FILT_FACTOR = 0.04;
ac = get_analysis_constants;

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;
BIN_SIZE = 0.050; % s
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

% 181022_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15
% counter_turn_epochs =  { { 2, { 11.03, 13.32 } }, ...
%                          { 3, { 16.07, 17.32 } }, ...
%                          { 4, { 24.26, 27.45 }  }, ...
%                          { 8, { 46.64, 48.28 } }, ...
%                          { 9, { 18.3, 19.52 } }, ...
%                          { 10, { 19.11, 20.83 } }};

for d = 1:length( directories )
    
    f = figure;
    
    cur_datapath = directories{ d }{ 1 };
    cur_sid      = directories{ d }{ 2 };
    
    datapath = [ basedir '/' cur_datapath ];
    analysis_path = [datapath '/analysis/'];
        
    stim_filepath = [analysis_path '/stim_window_data_sid_' num2str( cur_sid ) '.mat'];
    
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
    
    
    for tr = 1:size(ephys_data{1},1)
        
        cur_ephys = squeeze( ephys_data{1}(tr,:) ); 
        cur_yaw   = squeeze( bdata_vel{ 1 }( tr, ac.VEL_YAW, : ) );
        
        VmFilt_A2 = medfilt1( cur_ephys, FILT_FACTOR * ephys_SR, 'truncate' );
        VmFilt_A2_corr = VmFilt_A2 - mean(VmFilt_A2);
        
        t_down = squeeze(mean(reshape(ephys_time, [DT_EPHYS, length(ephys_time)/DT_EPHYS]), 1));
        A2_Vm_down = squeeze(mean(reshape(VmFilt_A2_corr, [ DT_EPHYS, length(VmFilt_A2)/DT_EPHYS ] ),1));
        
        yaw_t_down = squeeze(mean(reshape(bdata_vel_time, [DT_YAW, length(bdata_vel_time)/DT_YAW]),1));
        yaw_down = squeeze(mean(reshape(cur_yaw, [DT_YAW, length(bdata_vel_time)/DT_YAW]),1));
        
        hold on;
        
        SHIFT_FACTOR = 3;
        for ii = 1:(length(A2_Vm_down)-SHIFT_FACTOR)
            cur_yaw_1 = yaw_down( ii+SHIFT_FACTOR );
            plot(A2_Vm_down(ii), cur_yaw_1, 'o', 'MarkerSize', 3, 'color', 'b' );
        end
        
        for jj = 1:length( counter_turn_epochs )
            cur_tr = counter_turn_epochs{jj}{1};
            if( cur_tr == tr )
                cur_start_t = counter_turn_epochs{jj}{2}{1};
                cur_end_t = counter_turn_epochs{jj}{2}{2};
                
                cur_t_down = find( ( t_down >= cur_start_t ) & ( t_down <= cur_end_t ) );
                cur_yaw_t_down = find( ( yaw_t_down >= cur_start_t ) & ( yaw_t_down <= cur_end_t ) );
                
                for ii = 1:(length(cur_t_down)-SHIFT_FACTOR)
                    cur_index_yaw   = cur_yaw_t_down(ii);
                    cur_index_ephys = cur_t_down(ii);
                    cur_yaw_2 = yaw_down( cur_index_yaw +SHIFT_FACTOR );
                    plot(A2_Vm_down(cur_index_ephys), cur_yaw_2, 'o', 'MarkerSize', 10, 'color', 'r' );
                end
            end
        end
    end
end

xlabel('A2 Vm (mV)');
ylabel('Yaw (au)');

saveas(f,[analysis_path '/Vm_vs_yaw_spontaneous_CX_turning_' num2str(sid) '.fig']);
saveas(f,[analysis_path '/Vm_vs_yaw_spontaneous_CX_turning_' num2str(sid) '.png']);

% display linear regression slope of spontaneous vs. CX-evoked (connecting
% lines within a fly)

end

