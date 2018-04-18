%% Load data
clear all;

global slash;
slash = '/';

working_dir = '/Users/sasha/Dropbox/Wilson_lab/paper_1/data/descending_neurons/';
summary_analysis_path = [ working_dir '/summary_analysis' ];

COMPARE_ALL = 11;
COMPARE_BOTH = 12;
COMPARE_WAB320_WAB326 = 13;
ANALYSIS_TYPE = COMPARE_ALL;

if (ANALYSIS_TYPE == COMPARE_ALL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each category is a list of pairs {directory name, sid} representing
% behavior of a single fly. 
% 1. Both neurons labeled
% 2. Left neuron labeled
% 3. Right neuron labeled
% 4. None neurons labeled

LAL_DN_TYPE = 'A2';

% BOTH_LIST = {{'170209_hsFLP_ReachR_ss730_05', 3}, {'170209_hsFLP_ReachR_ss730_05', 5}, {'170210_hsFLP_ReachR_ss730_06', 3}, ...
%              {'170210_hsFLP_ReachR_ss730_06', 4}, {'170214_hsFLP_ReachR_ss730_08', 0}, {'170214_hsFLP_ReachR_ss730_08', 2}, ...
%              {'170215_hsFLP_ReachR_ss730_09', 1}, {'170215_hsFLP_ReachR_ss730_09', 3}, {'170216_hsFLP_ReachR_ss730_10', 0}, ...
%              {'170216_hsFLP_ReachR_ss730_10', 2}, {'170216_hsFLP_ReachR_ss730_10', 3}, {'170216_hsFLP_ReachR_ss730_10', 4}, ...
%              {'170217_hsFLP_ReachR_ss730_11', 1}, {'170220_hsFLP_ReachR_ss730_14', 0}, {'170220_hsFLP_ReachR_ss730_14', 1}, ...
%              {'170221_hsFLP_ReachR_ss730_15', 1}, {'170222_hsFLP_ReachR_ss730_16', 0}, {'170222_hsFLP_ReachR_ss730_16', 1}, ...
%              {'170222_hsFLP_ReachR_ss730_16', 2}, {'170222_hsFLP_ReachR_ss730_16', 3}, {'170222_hsFLP_ReachR_ss730_16', 5}, ...
%              {'170223_hsFLP_ReachR_ss730_17', 0}, {'170223_hsFLP_ReachR_ss730_17', 1}, {'170223_hsFLP_ReachR_ss730_17', 4} };

BOTH_LIST_SAME = {{'170214_hsFLP_ReachR_ss730_08', 0}, {'170215_hsFLP_ReachR_ss730_09', 1}, {'170216_hsFLP_ReachR_ss730_10', 2}, ...
                  {'170216_hsFLP_ReachR_ss730_10', 3}, {'170217_hsFLP_ReachR_ss730_11', 1}, {'170220_hsFLP_ReachR_ss730_14', 1}, ...
                  {'170222_hsFLP_ReachR_ss730_16', 5}, {'170222_hsFLP_ReachR_ss730_16', 0}, {'170222_hsFLP_ReachR_ss730_16', 3}, ...
                  {'170223_hsFLP_ReachR_ss730_17', 0}, {'170223_hsFLP_ReachR_ss730_17', 1}, {'170227_hsFLP_ReachR_ss730_20', 0}, ...
                  {'170227_hsFLP_ReachR_ss730_20', 1}, {'170227_hsFLP_ReachR_ss730_19', 2}, {'170228_hsFLP_ReachR_ss730_22', 0}, ...
                  {'170228_hsFLP_ReachR_ss730_21', 2}, {'170304_hsFLP_ReachR_ss730_25', 2}, {'170306_hsFLP_ReachR_ss730_26', 4}, ...
                  {'170306_hsFLP_ReachR_ss730_26', 1}, {'170306_hsFLP_ReachR_ss730_26', 3}, {'170306_hsFLP_ReachR_ss730_27', 3}, ...
                  {'170306_hsFLP_ReachR_ss730_27', 4}, {'170308_hsFLP_ReachR_ss730_28', 0}, {'170308_hsFLP_ReachR_ss730_29', 0}, ...
                  {'170308_hsFLP_ReachR_ss730_28', 1}, {'170308_hsFLP_ReachR_ss730_28', 2}, {'170308_hsFLP_ReachR_ss730_29', 3}, ...
                  {'170308_hsFLP_ReachR_ss730_28', 3}, {'170308_hsFLP_ReachR_ss730_29', 4}, {'170308_hsFLP_ReachR_ss730_29', 5}, ...
                  {'170308_hsFLP_ReachR_ss730_28', 7}
                  };
            
            
LEFT_LIST = {{'170208_hsFLP_ReachR_ss730_04', 0}, {'170210_hsFLP_ReachR_ss730_06', 0}, {'170220_hsFLP_ReachR_ss730_14', 2}, ...
             {'170223_hsFLP_ReachR_ss730_17', 3}, {'170228_hsFLP_ReachR_ss730_21', 0}, {'170228_hsFLP_ReachR_ss730_21', 1}, ...
             {'170306_hsFLP_ReachR_ss730_27', 1}, {'170315_hsFLP_ReachR_ss730_38', 1}, {'170322_hsFLP_ReachR_ss730_46', 0}
            };

RIGHT_LIST = {{'170215_hsFLP_ReachR_ss730_09', 0}, {'170222_hsFLP_ReachR_ss730_16', 4}, {'170227_hsFLP_ReachR_ss730_19', 0}, ...
              {'170308_hsFLP_ReachR_ss730_29', 1}, {'170313_hsFLP_ReachR_ss730_34', 2}, {'170320_hsFLP_ReachR_ss730_43', 1}, ...
              {'170330_hsFLP.PEST_ReachR_ss730_01', 0} };

NONE_LIST = {{'170211_hsFLP_ReachR_GFP_control_ss730_07', 0}, {'170217_hsFLP_ReachR_ss730_11', 3}, {'170218_hsFLP_ReachR_ss730_12', 0}, ...
                         {'170218_hsFLP_ReachR_ss730_12', 2}, {'170218_hsFLP_ReachR_ss730_12', 3}, {'170218_hsFLP_ReachR_ss730_12', 4}, ...
                         {'170219_hsFLP_ReachR_ss730_13', 1}, {'170219_hsFLP_ReachR_ss730_13', 2}                          
            };        
        
LAL_DN_labeling_categories_ss730 = { BOTH_LIST_SAME, LEFT_LIST, RIGHT_LIST, NONE_LIST };
                                 
LAL_DN_labeling_categories = LAL_DN_labeling_categories_ss730;  
string_per_label = {'both', 'left', 'right', 'none'};

elseif(ANALYSIS_TYPE == COMPARE_BOTH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAL_DN_TYPE = 'A2';
BOTH_LIST_SAME = {{'170214_hsFLP_ReachR_ss730_08', 0}, {'170215_hsFLP_ReachR_ss730_09', 1}, {'170216_hsFLP_ReachR_ss730_10', 2}, ...
                  {'170216_hsFLP_ReachR_ss730_10', 3}, {'170217_hsFLP_ReachR_ss730_11', 1}, {'170220_hsFLP_ReachR_ss730_14', 1}, ...
                  {'170222_hsFLP_ReachR_ss730_16', 5}, {'170222_hsFLP_ReachR_ss730_16', 0}, {'170222_hsFLP_ReachR_ss730_16', 3}, ...
                  {'170223_hsFLP_ReachR_ss730_17', 0}, {'170223_hsFLP_ReachR_ss730_17', 1}, {'170227_hsFLP_ReachR_ss730_20', 0}, ...
                  {'170227_hsFLP_ReachR_ss730_20', 1}, {'170227_hsFLP_ReachR_ss730_19', 2}, {'170228_hsFLP_ReachR_ss730_22', 0}
                };
BOTH_4_3 = { {'170216_hsFLP_ReachR_ss730_10', 0}, {'170215_hsFLP_ReachR_ss730_09', 3}, {'170220_hsFLP_ReachR_ss730_14', 0}, ...
             {'170222_hsFLP_ReachR_ss730_16', 2}, {'170223_hsFLP_ReachR_ss730_17', 4}
           };
       
BOTH_3_4 = { {'170214_hsFLP_ReachR_ss730_08', 2}, {'170216_hsFLP_ReachR_ss730_10', 4}, {'170221_hsFLP_ReachR_ss730_15', 1}, ...
             {'170222_hsFLP_ReachR_ss730_16', 1}
           };

LAL_DN_labeling_categories_ss730 = { BOTH_LIST_SAME, BOTH_4_3, BOTH_3_4 };
                                 
LAL_DN_labeling_categories = LAL_DN_labeling_categories_ss730;  
string_per_label = {'both_same', 'both_4_3', 'both_3_4'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(ANALYSIS_TYPE == COMPARE_WAB320_WAB326)
  
LAL_DN_TYPE = 'WAB_320_compare';

WAB_326_DATA = { {'170227_hsFLP_ReachR_ss730_20', 0}, {'170227_hsFLP_ReachR_ss730_20', 1}, {'170228_hsFLP_ReachR_ss730_22', 0}                
                };

WAB_320_RIGHT_LIST = {{'170215_hsFLP_ReachR_ss730_09', 0}, {'170222_hsFLP_ReachR_ss730_16', 4}, {'170227_hsFLP_ReachR_ss730_19', 0} };

LAL_DN_labeling_categories_ss730 = { WAB_320_RIGHT_LIST, WAB_326_DATA };
                                 
LAL_DN_labeling_categories = LAL_DN_labeling_categories_ss730;  
string_per_label = {'WAB320-right', 'WAB326', };

end

trial_type_cnt = 1;

first_stim_t = 3.0;
last_stim_t =  3.5;

TRIAL_CNT_MAX = 300;
s_settings = sensor_settings;
BEHAV_SAMPLE_RATE = s_settings.sensorPollFreq;

timenow_str = datestr(datetime, 'yymmdd_HHMMSS');

f = figure('units','normalized','outerposition',[0 0 1 1]); 
colors_per_label = {rgb('Black'), rgb('Red'), rgb('Green'), rgb('Blue') };
colors_per_label_sem = {rgb('Silver'), rgb('LightSalmon'), rgb('PaleGreen'), rgb('LightBlue') };
plt_hdl = [];
legend_text = {};
for labeling_type = 1:length(LAL_DN_labeling_categories)
    
    cur_fly_list = LAL_DN_labeling_categories{ labeling_type };  
    
    if( length( cur_fly_list ) == 0 )
        continue;
    end
    
    fly_data = {};
        
    for cur_fly = 1:length( cur_fly_list )
    
        cur_entry = cur_fly_list( cur_fly );
        cur_dir = cur_entry{1}{1};
        cur_sid = cur_entry{1}{2};
        
        datapath = [ working_dir cur_dir ];
        
        analysis_path = [datapath '/analysis/'];
        
        if(~exist(analysis_path, 'dir'))
            mkdir(analysis_path);
        end
    
        pre_generated_filepath = [analysis_path '/' cur_dir '_' num2str(cur_sid) '.mat'];
        
        fly_data{cur_fly} = load_causality_experiment_behavioral_data( cur_sid, datapath, pre_generated_filepath );
    end        
    
%     t_vel_all = [];
%     yaw_vel_all = [];
%     forward_vel_all = [];
    
    trial_cnt = 0;
    all_data_fwd = [];
    all_data_yaw = [];
    f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    for cur_fly = 1:length( cur_fly_list )
        
        cur_data = fly_data{cur_fly};
                
        cur_trial_fwd_avg = squeeze(mean( cur_data.forward_vel_all ));
        cur_trial_yaw_avg = squeeze(mean( cur_data.yaw_vel_all ));
              
        t_plot = cur_data.t_vel_all;

        subplot(2,1,1);
        hold on;
        %plot(t_plot, cur_trial_fwd_avg, 'color', colors_per_label{labeling_type});
        plot(t_plot, cur_trial_fwd_avg, 'LineWidth', 1.5 );
        ylabel('Fwd velocity (mm/s)');
        yy = ylim;
        y_min = yy(1); y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');        
        xlim([0 t_plot(end)]);
        set(gca(), 'FontSize', 16)
     
        subplot(2,1,2)
        hold on;
        %plot(t_plot, cur_trial_yaw_avg, 'color', colors_per_label{labeling_type});
        plot(t_plot, cur_trial_yaw_avg, 'LineWidth', 1.5);
        
        yy = ylim;
        y_min = yy(1); y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');        
        
        ylabel('Yaw velocity (deg/s)');
        xlabel('Time (s)');
        xlim([0 t_plot(end)]);
        set(gca(), 'FontSize', 16)
                
        all_data_fwd = vertcat( all_data_fwd, cur_data.forward_vel_all );
        
        start_baseline = (first_stim_t-0.5)*BEHAV_SAMPLE_RATE;
        end_baseline = (first_stim_t)*BEHAV_SAMPLE_RATE;
        
        yaw_vel_baseline_corr = cur_data.yaw_vel_all - repmat(mean(cur_data.yaw_vel_all(:,start_baseline:end_baseline),2), [1 size(cur_data.yaw_vel_all,2)]);
        
        %all_data_yaw = vertcat( all_data_yaw, cur_data.yaw_vel_all );
        all_data_yaw = vertcat( all_data_yaw, yaw_vel_baseline_corr );
    end
    subplot(2,1,1);
    title( ['N ' string_per_label{labeling_type} ' trials = ' num2str(size(cur_data.forward_vel_all,1))] );
    
    saveas(f1, [summary_analysis_path '/' timenow_str '_LAL_DN_' LAL_DN_TYPE '_causality_' string_per_label{ labeling_type } '.fig' ]);
    saveas(f1, [summary_analysis_path '/' timenow_str '_LAL_DN_' LAL_DN_TYPE '_causality_' string_per_label{ labeling_type } '.png' ]);
    
    close(f1);
    
    figure(f)
    subplot(2,1,1);
    all_data_fwd_avg = mean(all_data_fwd);
    all_data_fwd_sem = get_sem(all_data_fwd,1);
    
    hold on;
    t_plot = cur_data.t_vel_all(1,:);
    fh = fill( [t_plot, fliplr(t_plot)], ...
               [(all_data_fwd_avg-all_data_fwd_sem) fliplr((all_data_fwd_avg+all_data_fwd_sem))], ...
               colors_per_label_sem{labeling_type});
    set(fh, 'EdgeColor', 'None');

    plot( t_plot, all_data_fwd_avg, 'color', colors_per_label{labeling_type}, 'LineWidth', 1.5);
    xlim([0 t_plot(end)]);

    yy = ylim;
    y_min = yy(1); y_max = yy(2);
    hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
    set(gca,'children',circshift(get(gca,'children'),-1));
    set(hh, 'EdgeColor', 'None');
    
    ylabel('Fwd velocity (mm/s)');
    set(gca(), 'FontSize', 16)
    
    subplot(2,1,2);
    hold on;
    
    %all_data_yaw_avg_tmp = mean(all_data_yaw);
    %all_data_yaw_avg_baseline = mean(all_data_yaw_avg_tmp);
    %all_data_yaw_avg = all_data_yaw_avg_tmp - repmat(all_data_yaw_avg_baseline, [1 size(all_data_yaw_avg_tmp,2)]);
    
    all_data_yaw_avg = mean(all_data_yaw);
    all_data_yaw_sem = get_sem(all_data_yaw,1);    
    fh = fill( [t_plot, fliplr(t_plot)], ...
               [(all_data_yaw_avg-all_data_yaw_sem) fliplr((all_data_yaw_avg+all_data_yaw_sem))], ...
               colors_per_label_sem{labeling_type} );
    set(fh, 'EdgeColor', 'None');
    
    plt_hdl(end+1) = plot( t_plot, all_data_yaw_avg, 'color', colors_per_label{labeling_type}, 'LineWidth', 1.5);
    legend_text{end+1} = ['N ' string_per_label{labeling_type} ' trials ( ' num2str( size(all_data_fwd,1) ) ' ) from ' num2str(length( cur_fly_list)) ' flies'];
    xlim([0 t_plot(end)]);

    yy = ylim;
    y_min = yy(1); y_max = yy(2);
    hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
    set(gca,'children',circshift(get(gca,'children'),-1));
    set(hh, 'EdgeColor', 'None');

    set(gca(), 'FontSize', 16)
    ylabel('Yaw velocity (deg/s)');
    xlabel('Time (s)');
end

legend(plt_hdl, legend_text(:));

save([summary_analysis_path '/' timenow_str '_LAL_DN_causality_directories_to_analyze.mat'],'LAL_DN_labeling_categories');

saveas(f, [summary_analysis_path '/' timenow_str '_LAL_DN_' LAL_DN_TYPE '_causality.fig' ]);
saveas(f, [summary_analysis_path '/' timenow_str '_LAL_DN_' LAL_DN_TYPE '_causality.png' ]);
