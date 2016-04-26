function generate_tt_trial_single_plane( asid, sid, cdata_raw, bdata_vel, btrial_meta, bdata_vel_time, FPS, analysis_path, rois )

global slash;

% tbt = trial-by-trial
tbt_analysis_path = [analysis_path '/asid_' num2str(asid) '_trial_by_trial/'];

if(~exist(tbt_analysis_path, 'dir'))
    mkdir( tbt_analysis_path );
end

aconstants = get_analysis_constants;
trial_type_cnt = length(cdata_raw);
settings = sensor_settings;

ac = get_analysis_constants;
order = ac.order;

% Assume 12 ROIs
SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

for trial_type = 1:trial_type_cnt
    
    cur_trial_cnt = size( cdata_raw{ trial_type }, 1 );
    cur_trial_type_str = aconstants.task_str{trial_type};
    
    for trial = 1:cur_trial_cnt
        f = figure('units','normalized','outerposition',[0 0 1 1]);
            
        subaxis(4, 3, [4 7], 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);

        %cur_trial_cdata = squeeze(cdata_raw{ trial_type }(trial,:,:,:,:,:));
        cur_trial_bdata = squeeze(bdata_vel{ trial_type }(trial,:,:));
        cur_trial_id = squeeze(btrial_meta{ trial_type }(trial, 2));
        
        cur_trial_cdata = cdata_raw{ trial_type }(trial,:,:,:);

        ref_img = squeeze(mean(cur_trial_cdata(:,:,3:end),3));
        [xsize, ysize] = size(ref_img);
        imagesc( ref_img );
        colormap gray;
        caxis([0 1500]);
        axis image;
        axis off;
        
        colorindex = 0;
        
        cur_roi_idx = 1;
        cur_roi = rois{ 1, cur_roi_idx };
        
        while( ~isempty( cur_roi ) )
            xv = cur_roi(:,1);
            yv = cur_roi(:,2);
            
            hold on;
            currcolor    = order(1+mod(colorindex,size(order,1)),:);
            plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
            text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
            
            cur_roi_idx = cur_roi_idx + 1;
            if( cur_roi_idx > size(rois,2))
                break;
            end
            
            cur_roi = rois{1, cur_roi_idx};
            colorindex = colorindex + 1;
        end
        
        tt = title(filename_prefix);
        set(tt, 'Interpreter', 'none')
            
        % 1 2 3 
        % 4 5 6 
        % 7 8 9
        % 10 11 12 

        % Generate the time courses
        colorindex = 0;
        for sa = 1:3
            if (sa == 1)
                subaxis(4, 3, [2 3], 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
            elseif( sa == 2 )
                subaxis(4, 3, [5 6], 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
            elseif( sa == 3 )
                subaxis(4, 3, [8 9], 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN); 
            end

            for roi_ord = 1:4
                roi_idx = (sa-1)*4 + roi_ord;

                cur_roi = rois{ 1, roi_idx };
        
                while( ~isempty( cur_roi ) )
                    xv = cur_roi(:,1);
                    yv = cur_roi(:,2);
                    
                    inpoly = inpolygon(x,y,xv,yv);
                    
                    currcolor = order(1+mod(colorindex,size(order,1)),:);
                    
                    tmp = squeeze(sum(sum(cur_trial_cdata .* repmat(inpoly, [1, 1, nframes])))) / sum(inpoly(:));
                    baseline = repmat(mean(tmp(base_begin:base_end)), [1 1 size(tmp,2)]);
                    itrace = (tmp-baseline) ./ baseline;
                    
                    hold on;
                    plot(t, itrace,'color', currcolor);
                    
                    cur_roi_idx = cur_roi_idx + 1;
                    if( cur_roi_idx > size(rois,2))
                        break;
                    end
                    
                    cur_roi = rois{1, cur_roi_idx};
                    colorindex = colorindex + 1;
                end
            end
        
            yy = ylim;
            y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
            hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
            set(gca,'children',circshift(get(gca,'children'),-1));
            set(hh, 'EdgeColor', 'None');
            
            xlim([0, total_time]);
            if(mod((p-1),IMAGE_COLS) == 0 )
                ylabel('dF/F');
            else
                set(gca, 'YTickLabel', '');
            end
    
            drawnow;
        end

        % Axis for behavioral data
        subaxis( 4, 3, [11 12], 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
        
        hold on;
        p1 = plot(bdata_vel_time, cur_trial_bdata(ac.VEL_FWD,:), 'color', rgb('FireBrick'));
        p2 = plot(bdata_vel_time, cur_trial_bdata(ac.VEL_YAW,:), 'color', rgb('SeaGreen'));
        
        legend( [ p1, p2 ], 'Vel fwd', 'Vel yaw' );
        
        yy = ylim;
        y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');
        
        xlim([0, total_time]);
        xlabel('Time (s)');
        ylabel('Velocity (au/s)');
        
        drawnow;
    
        cur_tbt_filename_prefix = [ tbt_analysis_path '/time_courses_in_roi_sid_' num2str(sid) '_' cur_trial_type_str '_tid_' num2str(cur_trial_id)];        
        saveas( f, [cur_tbt_filename_prefix '.fig'] );
        saveas( f, [cur_tbt_filename_prefix '.png'] );
    end
end
end

