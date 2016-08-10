function [ clusters ] = generate_turning_clusters( sid, bdata_vel_time, bdata_vel, analysis_path, max_clusters )

ac = get_analysis_constants();

settings = sensor_settings;
prestim = settings.pre_stim;
stim    = settings.stim;

trial_cnt = size( bdata_vel, 2 );

%clust_period = find( (bdata_vel_time > (prestim-0.5)) & (bdata_vel_time < (prestim+stim+1.5)));
clust_period = find( (bdata_vel_time > (prestim-0.5)) & (bdata_vel_time < (prestim+stim+1.5)));
t = bdata_vel_time(clust_period);

f1 = figure();

MAX_CLUST = max_clusters;

clusters = cell(1,trial_cnt);

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
    
    clusters{ trial_type } = T;
    
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