function generate_volume_time_courses_with_rois( trial_cdata, trial_bdata_vel, bdata_vel_time, rois, VPS, filename_prefix )

ac = get_analysis_constants;
settings = sensor_settings;
order = ac.order;

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

PLANES = size(trial_cdata,3);

if(PLANES == 8)
    IMAGE_ROWS = 2;
    IMAGE_COLS = 4;
elseif(PLANES == 12)
    IMAGE_ROWS = 3;
    IMAGE_COLS = 4;
else
    IMAGE_ROWS = floor(sqrt(PLANES));
    IMAGE_COLS = IMAGE_ROWS;
end

prestim = settings.pre_stim;
stim    = settings.stim;

base_begin = 1;
base_end = floor(prestim*VPS);

poststim    = settings.post_stim;
total_time = prestim + stim + poststim;

first_stim_t = prestim;
last_stim_t = stim + prestim;

f = figure('units','normalized','outerposition',[0 0 1 1]);

[ysize, xsize] = size(squeeze(trial_cdata(:,:,1,1)));
[x, y] = meshgrid(1:xsize, 1:ysize);

nframes = size(trial_cdata,4);
t = [1:nframes]./VPS;

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
    
    if( PLANES == 16 )
        if(p>=1 & p<=8 )
            ylim([-1 2.5]);
        elseif (p>=9 & p<=12 )
            ylim([-1 2.5]);
        elseif (p>=13 & p<=16 )
            ylim([-1 2.5]);
        end
    elseif( PLANES == 8 ) 
        if(p>=1 & p<=4 )
            ylim([-0.6 2.0]);
        elseif (p>=5 & p<=8 )
            ylim([-0.6 2.0]);
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
    
    set(gca, 'XTickLabel', '');
    
    if( p == 2 )
        tt = title(filename_prefix);
        set(tt, 'Interpreter', 'none')
    end
    
    drawnow;
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

saveas(f, [filename_prefix '.fig']);
saveas(f, [filename_prefix '.png']);

close(f);

end

