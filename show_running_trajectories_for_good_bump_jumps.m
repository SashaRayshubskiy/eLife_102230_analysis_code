function show_running_trajectories_for_good_bump_jumps(  sid, bdata_vel_time, bdata_vel, trial_stim_id_map, pico_stim_data, analysis_path )

settings = sensor_settings;
ephys_FR = settings.sampRate;
ball_FR = settings.sensorPollFreq;

traj = get_single_trial_trajectories( sid, bdata_vel_time, bdata_vel );

num_trials = size(traj{1},1);

%TRIAL_MAX = 6;
TRIAL_MAX = num_trials;

EBYE = EB_yaw_ephys_data; 
% clean_bump_return_stims = EBYE.Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_16_181203_clean_bump_return(:,1);
clean_bump_return_stims = [];

f = figure;
cm = colormap(jet(TRIAL_MAX));
hold on; 
for tr = 1:TRIAL_MAX
    disp_x = squeeze(traj{1}(tr,1,:));
    disp_y = squeeze(traj{1}(tr,2,:));       
    plot(disp_x, disp_y, 'color', cm(tr,:), 'DisplayName', ['trial: ' num2str(tr)]);
        
    if( length(trial_stim_id_map) > 0 )
        cur_stim_meta = trial_stim_id_map{tr};
    end

    if( length(clean_bump_return_stims) == 0 )
        cur_stim = 10.0*squeeze(pico_stim_data{1}( tr, : ));
        stim_idx = find(diff(cur_stim) > 3.0 );

        stim_idx_yaw = floor(( stim_idx / ephys_FR ) * ball_FR);
        stim_x = disp_x( stim_idx_yaw );
        stim_y = disp_y( stim_idx_yaw );
        plot(stim_x, stim_y, 'x', 'MarkerSize', 10.0, 'color', 'r', 'LineWidth', 2.0);
    else
        
        j = 1;
        while( j <= length(cur_stim_meta) )
            
            cur_global_stim_idx = cur_stim_meta(j,1);
            % Check if the stim caused a clean bump jump/return
            if( any( clean_bump_return_stims == cur_global_stim_idx )  == 1 )
                
                cur_stim_idx_in_trial = cur_stim_meta(j,2);
                
                stim_idx_yaw = floor(( cur_stim_idx_in_trial / ephys_FR ) * ball_FR);
                
                stim_x = disp_x( stim_idx_yaw );
                stim_y = disp_y( stim_idx_yaw );
                
                plot(stim_x, stim_y, 'x', 'MarkerSize', 10.0, 'color', 'r', 'LineWidth', 2.0);
            end
            
            j = j + 1;
        end        
    end
    
    waitforbuttonpress;
    cla();    
end

legend();
xlabel('X displacement (au)');
ylabel('Y displacement (au)');

saveas(f,[analysis_path '/running_trajectories_' num2str(sid) '_max_trials_' num2str( TRIAL_MAX ) '.fig']);
saveas(f,[analysis_path '/running_trajectories_' num2str(sid) '_max_trials_' num2str( TRIAL_MAX ) '.png']);

end
