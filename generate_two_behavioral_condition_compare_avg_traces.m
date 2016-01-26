function [btraces_per_condition, ctraces_in_roi_per_condition] = collect_two_behavioral_condition_traces( condition_trials, condition_trials_str, atid, sid, cdata_raw, bdata_vel, btrial_meta, bdata_vel_time, VPS, filename_prefix, rois )

ac = get_analysis_constants;
settings = sensor_settings;
order = ac.order;

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

IMAGE_ROWS = 4;
IMAGE_COLS = 4;
PLANES = IMAGE_ROWS * IMAGE_COLS;

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

base_begin = 1;
base_end = floor(prestim*VPS);

total_time = prestim + stim + poststim;

first_stim_t = prestim;
last_stim_t = stim + prestim;

[ysize, xsize] = size(squeeze(cdata_raw{1}(1,:,:,1,1,1)));
[x, y] = meshgrid(1:xsize, 1:ysize);

nframes = size( squeeze(cdata_raw{1}(1,:,:,1,1,:)), 3 );
t = [1:nframes]./VPS;

for trial_type = 1:size(condition_trials{1},1)
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);    
    
    for cond_ord = 1:size(condition_trials,1)
        
        cur_condition_trial_ords = condition_trials{cond_ord, trial_type};
        
        for trial_ord_idx = 1:size(cur_condition_trial_ords,1)
        
            cur_trial_ord = cur_condition_trial_ords( trial_ord_idx );
            
            trail_cdata = squeeze(cdata_raw{trial_type}(cur_trial_ord,:,:,:,:,:));            
        
            % Collect behavioral data
            
            % Collect calcium data
            for p=1:PLANES
                
                % Axis for time courses
                subaxis(IMAGE_ROWS+1, IMAGE_COLS, p, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
                
                colorindex = 0;
                
                cur_roi_idx = 1;
                cur_roi = rois{ p, cur_roi_idx };
                
                cdata_in_plane = squeeze(trial_cdata(:,:,p,:));
                
                while( ~isempty( cur_roi ) )
                    xv = cur_roi(:,1);
                    yv = cur_roi(:,2);
                    
                    inpoly = inpolygon(x,y,xv,yv);
                    
                    currcolor = order(1+mod(colorindex,size(order,1)),:);
                    
                    tmp = squeeze(sum(sum(cdata_in_plane .* repmat(inpoly, [1, 1, nframes])))) / sum(inpoly(:));
                    baseline = repmat(mean(tmp(base_begin:base_end)), [1 1 size(tmp,2)]);
                    itrace = (tmp-baseline) ./ baseline;
                    
                    hold on;
                    plot(t, itrace,'color', currcolor);
                    
                    cur_roi_idx = cur_roi_idx + 1;
                    if( cur_roi_idx > size(rois,2))
                        break;
                    end
                    
                    cur_roi = rois{p, cur_roi_idx};
                    colorindex = colorindex + 1;
                end                
            end
        
        for c = 1:IMAGE_COLS
            
            % Axis for behavioral data
            subaxis(IMAGE_ROWS+1, IMAGE_COLS, PLANES + c, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
            
            hold on;
            p1 = plot(bdata_vel_time, trial_bdata_vel(ac.VEL_FWD,:), 'color', rgb('FireBrick'));
            p2 = plot(bdata_vel_time, trial_bdata_vel(ac.VEL_YAW,:), 'color', rgb('SeaGreen'));
            
            if( c == 1 )
                legend( [ p1, p2 ], 'Vel fwd', 'Vel yaw' );
            end
            
            yy = ylim;
            y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
            hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
            set(gca,'children',circshift(get(gca,'children'),-1));
            set(hh, 'EdgeColor', 'None');
            
            xlim([0, total_time]);
            xlabel('Time (s)');
            if( c == 1 )
                ylabel('Velocity (au/s)');
            else
                set(gca, 'YTickLabel', '');
            end
            
            drawnow;
        end
    end
    
    saveas(f, [ filename_prefix '_' ac.task_str{trial_type} '.fig' ]);
    saveas(f, [ filename_prefix '_' ac.task_str{trial_type} '.png' ]);
    close(f);
end
end

