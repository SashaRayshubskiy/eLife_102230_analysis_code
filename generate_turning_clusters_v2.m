function [ clusters ] = generate_turning_clusters_v2( sid, bdata_vel_time, bdata_vel, analysis_path, max_clusters, idx )

ac = get_analysis_constants();

settings = sensor_settings;
prestim = settings.pre_stim;
stim    = settings.stim;

trial_cnt = size( bdata_vel, 2 );

%clust_period = find( (bdata_vel_time > (prestim-0.5)) & (bdata_vel_time < (prestim+stim+1.5)));
%clust_period = find( (bdata_vel_time > (prestim-0.5)) & (bdata_vel_time < (prestim+stim+1.5)));
%clust_period = find( (bdata_vel_time > (prestim-0.2)) & (bdata_vel_time < (prestim+0.0)));
clust_period = find( (bdata_vel_time > (prestim)) & (bdata_vel_time < (prestim+0.5)));
t = bdata_vel_time(clust_period);

plot_period = find( (bdata_vel_time > (prestim-1.0)) & (bdata_vel_time < (prestim+stim+1.5)));
t_p = bdata_vel_time(plot_period);

f1 = figure();

MAX_CLUST = max_clusters;

clusters = cell(1,trial_cnt);

for trial_type = 1:trial_cnt
    
    yaw_data = squeeze(bdata_vel{ trial_type }( :, ac.VEL_YAW, clust_period));
    yaw_data_p = squeeze(bdata_vel{ trial_type }( :, ac.VEL_YAW, plot_period));
    
    if 1
    T = clusterdata( yaw_data, 'distance', 'seuclidean', 'linkage', 'ward', 'maxclust', MAX_CLUST );
    %T = clusterdata( yaw_data, 'linkage', 'ward', 'maxclust', MAX_CLUST );
    
    clusters{ trial_type } = T;
    
    for c = 1:MAX_CLUST

        subplot(3,1,c);
        
        cur_cluster = find(T == c);
        
        %cur_cluster_trials = yaw_data(cur_cluster, :);
        cur_cluster_trials = yaw_data_p(cur_cluster, :);
        
        hold on;
        for i=1:size(cur_cluster_trials,1)
            plot(t_p, squeeze(cur_cluster_trials(i,:)), 'color', ac.single_clr_by_type{ trial_type }, 'LineWidth', 0.5);
        end
                        
        plot(t_p, squeeze(mean(cur_cluster_trials)), 'LineWidth', 2.0, 'color', ac.clr_by_type{ trial_type } );
        
        title([ac.task_str{trial_type} ' ( ' num2str(size(cur_cluster_trials,1)) ' )']);
        xlim([t_p(1) t_p(end)]);
        ylim([-0.45 0.45]);
        %ylim([-0.3 0.3]);
    end
    end
end

%title([ 'clust_period = ' num2str(clust_period) ] );
saveas(f1, [analysis_path '/clusterdata_v2_maxclust_' num2str(MAX_CLUST) '_sid_' num2str(sid) '_' num2str(idx) '.png']);
saveas(f1, [analysis_path '/clusterdata_v2_maxclust_' num2str(MAX_CLUST) '_sid_' num2str(sid) '_' num2str(idx) '.fig']);

end