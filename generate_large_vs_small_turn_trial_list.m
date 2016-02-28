function [ condition_trials, condition_trials_str, condition_str ] = generate_large_vs_small_turn_trial_list( sid, bdata_vel_time, bdata_vel, analysis_path )

ac = get_analysis_constants();
condition_trials_str = { 'large_turns', 'small_turns' };
condition_str = 'large_vs_small_turning';

settings = sensor_settings;
prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

total_time = prestim + stim + poststim;

trial_cnt = size( bdata_vel, 2 );
condition_trials = cell(trial_cnt,2);

%clust_period = find( (bdata_vel_time > (prestim-0.5)) & (bdata_vel_time < (prestim+stim+1.5)));
clust_period = find( (bdata_vel_time > (prestim)) & (bdata_vel_time < (prestim+stim+1.5)));
t = bdata_vel_time(clust_period);

if 1
f = figure;

turn_period = find( (bdata_vel_time > (prestim)) & (bdata_vel_time < (prestim+0.8)));
counter_turn_period = find( (bdata_vel_time > (prestim+0.4)) & (bdata_vel_time < (prestim+1.4)));

turn_metadata = cell(1, trial_cnt);
yaw_freq_cutoff = 3.0; % Hz

for trial_type = 1:trial_cnt    
    
    subplot(2,1,trial_type);
        
    yaw_data = squeeze(bdata_vel{ trial_type }( :, ac.VEL_YAW, : ));
    
    turn_metadata{ trial_type } = zeros(size(yaw_data,1), 4);
    
    hold on;
    for i=1:size(yaw_data,1)
        cur_trial_data_1 = squeeze(yaw_data(i,:));
        
        cur_trial_data = fft_filter( cur_trial_data_1', yaw_freq_cutoff, settings.sensorPollFreq );
        
        % Extract key trial parameters from yaw data
        % 1. Max magnitude and time of turn
        % 2. Max magnitude and time of counter turn
        cur_turn_period = cur_trial_data( turn_period );
        cur_counter_turn_period = cur_trial_data( counter_turn_period );  

        if( trial_type == ac.LEFT )            
            [turn_mag, turn_i] = min(cur_turn_period);                        
            [counter_turn_mag, counter_turn_i] = max(cur_counter_turn_period);                        
        else
            [turn_mag, turn_i] = max(cur_turn_period);
            [counter_turn_mag, counter_turn_i] = min(cur_counter_turn_period);                        
        end
        
        turn_i = turn_i + turn_period( 1 );
        counter_turn_i = counter_turn_i + counter_turn_period( 1 );

        turn_t = bdata_vel_time( turn_i );
        counter_turn_t = bdata_vel_time( counter_turn_i );            
        turn_metadata{ trial_type }(i, :) = [ turn_t, turn_mag, counter_turn_t, counter_turn_mag ];
        
        %plot(t, cur_trial_data(clust_period), 'color', ac.single_clr_by_type{ trial_type }, 'LineWidth', 0.5);
        %plot(t, cur_trial_data_filt(clust_period), 'color', 'b', 'LineWidth', 0.5, 'LineStyle', '--');
        plot( turn_t, turn_mag, 'X', 'color', rgb('Maroon'), 'LineWidth', 0.5 );
        plot( counter_turn_t, counter_turn_mag, 'X', 'color', rgb('Indigo'), 'LineWidth', 0.5 );
        
        % waitforbuttonpress;
    end
        
    plot(t, squeeze(mean(yaw_data(:,clust_period))), 'LineWidth', 2.0, 'color', ac.clr_by_type{ trial_type } );
    
    title(ac.task_str{trial_type});
    xlim([t(1) t(end)]);
    ylim([-0.35 0.35]);
end

saveas(f, [analysis_path '/turn_meta_all_trials_to_cluster_sid_' num2str(sid) '.png']);
saveas(f, [analysis_path '/turn_meta_all_trials_to_cluster_sid_' num2str(sid) '.fig']);
%saveas(f, [analysis_path '/all_trials_to_cluster_sid_' num2str(sid) '.png']);
%saveas(f, [analysis_path '/all_trials_to_cluster_sid_' num2str(sid) '.fig']);
end

NUM_BIN = 30;

f2 = figure; 
subplot(2,4,1)
histogram(turn_metadata{1}(:,1));
title('Left: Turning times');

subplot(2,4,2)
histogram(turn_metadata{1}(:,2));
title('Left: Turning magnitude');

subplot(2,4,3)
histogram(turn_metadata{1}(:,3));
title('Left: Counter turn times');

subplot(2,4,4)
histogram(turn_metadata{1}(:,4));
title('Left: Counter turn magnitude');

subplot(2,4,5)
histogram(turn_metadata{2}(:,1));
title('Right: Turning times');

subplot(2,4,6)
histogram(turn_metadata{2}(:,2));
title('Right: Turning magnitude');

subplot(2,4,7)
histogram(turn_metadata{2}(:,3));
title('Right: Counter turn times');

subplot(2,4,8)
histogram(turn_metadata{2}(:,4));
title('Right: Counter turn magnitude');

saveas(f2, [analysis_path '/turn_metadata_histograms_sid_' num2str(sid) '.png']);
saveas(f2, [analysis_path '/turn_metadata_histograms_sid_' num2str(sid) '.fig']);

if 0
f1 = figure();

MAX_CLUST = 5;

for trial_type = 1:trial_cnt
    
    yaw_data = squeeze(bdata_vel{ trial_type }( :, ac.VEL_YAW, clust_period));
    
    if 0
        yaw_data_1 = zeros( size(bdata_vel{ trial_type },1), size(bdata_vel{ trial_type },3) );
        for tr = 1:size(bdata_vel{ trial_type },1)
            yaw_data_1(tr,:) = fft_filter( squeeze(bdata_vel{ trial_type }( tr, ac.VEL_YAW, :)), yaw_freq_cutoff, settings.sensorPollFreq );
        end
        
        yaw_data = yaw_data_1( :, clust_period );
    end
    
    if 0
        Y = pdist(yaw_data);
        Z = linkage(Y);
        dendrogram(Z,4);
    end
    
    if 1
    T = clusterdata( yaw_data, 'distance', 'seuclidean', 'linkage', 'ward', 'maxclust', MAX_CLUST );
    %T = clusterdata( yaw_data, 'linkage', 'ward', 'maxclust', MAX_CLUST );
    
    for c = 1:MAX_CLUST
        if(trial_type == 1)
            subplot( MAX_CLUST, 2, 2*(c-1)+1 );
        else
            subplot( MAX_CLUST, 2, 2*c );
        end
        
        cur_cluster = find(T == c);
        
        cur_cluster_trials = yaw_data(cur_cluster, :);
        hold on;
        for i=1:size(cur_cluster_trials,1)
            plot(t, squeeze(cur_cluster_trials(i,:)), 'color', ac.single_clr_by_type{ trial_type }, 'LineWidth', 0.5);
        end
                        
        plot(t, squeeze(mean(cur_cluster_trials)), 'LineWidth', 2.0, 'color', ac.clr_by_type{ trial_type } );
        
        title([ac.task_str{trial_type} ' ( ' num2str(size(cur_cluster_trials,1)) ' )']);
        xlim([t(1) t(end)]);
        ylim([-0.35 0.35]);
    end
    end
end

saveas(f1, [analysis_path '/clusterdata_maxclust_' num2str(MAX_CLUST) '_sid_' num2str(sid) '.png']);
saveas(f1, [analysis_path '/clusterdata_maxclust_' num2str(MAX_CLUST) '_sid_' num2str(sid) '.fig']);
end

end