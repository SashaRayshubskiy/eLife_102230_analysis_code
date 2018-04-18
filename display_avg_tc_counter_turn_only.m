function display_avg_tc_counter_turn_only(cdata_raw, bdata_vel_time, bdata_vel, VPS, analysis_dir)

ac = get_analysis_constants;

settings = sensor_settings;

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

total_time = prestim + stim + poststim;
first_stim_t = prestim;
last_stim_t = stim + prestim;

nframes = size(cdata_raw{1}, 5);

t = (([0:nframes-1]))./VPS;

x_size  = size(cdata_raw{1}, 2);
y_size  = size(cdata_raw{1}, 3);
[x, y] = meshgrid(1:y_size, 1:x_size);
colorindex = 0;

order    = ac.order;

cur_trial = squeeze(cdata_raw{1}(1,:,:,:,:));

cur_trial_single_plane = squeeze(mean(cur_trial(:,:,5:end,:),3));

%%%
% Generate trials where the fly counter turns, and does not counter turn
%%%
ct_list = generate_counter_turn_trials( bdata_vel_time, bdata_vel, analysis_dir );

%%%
% Take an avg in an ROI
%%%
ff = figure;
imagesc(mean(cur_trial_single_plane, 3));

[xv, yv] = (getline(gca, 'closed'));
inpoly = inpolygon(x,y,xv,yv);

%draw the bounding polygons and label them
currcolor    = order(1+mod(colorindex,size(order,1)),:);
hold on;
plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
%%%

f = figure;
for trial_type = 1:2
    
    if( trial_type == ac.BOTH )
        cur_color_avg = rgb('DimGray');
        cur_color_single = rgb('DarkGray');
    elseif( trial_type == ac.LEFT )
        cur_color_avg = rgb('FireBrick');
        cur_color_single = rgb('LightSalmon');
    elseif( trial_type == ac.RIGHT )
        cur_color_avg = rgb('SeaGreen');
        cur_color_single = rgb('PaleGreen');
    end
    
    % subplot(2,1,trial_type)
    
    
    for cond = 1:2
        
        if(cond == 1)
            cur_style = '-';
        else
            cur_style = '--';            
        end
        
        cur_df_f = [];
        for tr_ord = 1:length(ct_list{trial_type, cond})
            
            tr = ct_list{ trial_type, cond }( tr_ord );
            
            cur_trial = squeeze( cdata_raw{trial_type}(tr,:,:,:,:) );
            
            %cur_trial_single_plane = squeeze(mean(cur_trial(:,:,9:11,:),3));
            cur_trial_single_plane = squeeze(mean(cur_trial(:,:,5:end,:),3));
            
            cur_trial_all_plane_tc = squeeze(sum(sum(double(cur_trial_single_plane).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
            
            mean_baseline = repmat(mean(cur_trial_all_plane_tc(1:8)), [1 length(cur_trial_all_plane_tc)]);
            cur_trial_all_plane_tc_df_f = (cur_trial_all_plane_tc' - mean_baseline) ./ mean_baseline;
            
            cur_df_f(tr,:) = cur_trial_all_plane_tc_df_f;
            
            % hold on;
            % plot(cur_trial_all_plane_tc_corr);
        end
        
        
        avg_tc = squeeze(mean(cur_df_f));
        std_tc = std(cur_df_f,1);
        hold on;
        %     fh = fill( [t, fliplr(t)], ...          [(avg_tc+std_tc) fliplr(avg_tc-std_tc)], ...
        %         cur_color_single );
        %     set(fh, 'EdgeColor', 'None');
        
        plot(t, avg_tc, 'color', cur_color_avg, 'LineWidth', 2.0, 'LineStyle', cur_style);
        xlim([0 t(end)]);
    end
end

yy = ylim;
y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

saveas(f, [analysis_dir '/avg_tc.fig']);
saveas(f, [analysis_dir '/avg_tc.png']);
end

