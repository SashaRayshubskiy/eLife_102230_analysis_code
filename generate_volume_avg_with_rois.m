function generate_volume_avg_with_rois( trial_cdata, rois, filename_prefix )

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

f = figure('units','normalized','outerposition',[0 0 1 1]);
ac = get_analysis_constants;
order = ac.order;

for i=1:PLANES
    all_axes(i) = subaxis(IMAGE_ROWS, IMAGE_COLS, i, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
    %all_axes(i) = subplot(IMAGE_ROWS, IMAGE_COLS, i);
    
    % Throw away the first few volumes due to settling time.
    ref_img = mean(squeeze(trial_cdata(:,:,i,3:end)),3);
    [xsize, ysize] = size(ref_img);
    %imagesc(imresize(ref_img, [xsize 2*ysize]));
    imagesc( ref_img );
    colormap gray;
    caxis([0 1500]);
    axis image;
    axis off;
    
    colorindex = 0;
    
    cur_roi_idx = 1;
    cur_roi = rois{ i, cur_roi_idx };
    
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
        
        cur_roi = rois{i, cur_roi_idx};
        colorindex = colorindex + 1;
    end
    
    if( i == 2 )
        tt = title(filename_prefix);
        set(tt, 'Interpreter', 'none')
    end    

    drawnow;
end

saveas(f, [filename_prefix '.fig']);
saveas(f, [filename_prefix '.png']);

close(f);

end

