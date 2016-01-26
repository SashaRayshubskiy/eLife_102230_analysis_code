function [ new_rois ] = add_rois_from_volume( cdata, rois, plane_to_append )

new_rois = rois;

f = figure('units','normalized','outerposition',[0 0 1 1]);

% Throw away the first few volumes due to settling time.
imagesc(mean(squeeze(cdata(:,:,plane_to_append,3:end)),3));
colormap gray;
caxis([0 3000]);
axis image;
axis off;

ac = get_analysis_constants;
order = ac.order;

colorindex = 0;

cur_roi_idx = 1;
cur_roi = new_rois{ plane_to_append, cur_roi_idx };

while( ~isempty( cur_roi ) )
    xv = cur_roi(:,1);
    yv = cur_roi(:,2);
    
    hold on;
    currcolor    = order(1+mod(colorindex,size(order,1)),:);
    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    
    cur_roi_idx = cur_roi_idx + 1;
    cur_roi = new_rois{plane_to_append, cur_roi_idx};
    colorindex = colorindex + 1;
end

% Add new ROIs here
while(1)
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
        break
    end
    
    %draw the bounding polygons and label them
    hold on;
    currcolor    = order(1+mod(colorindex,size(order,1)),:);
    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    
    new_rois{ plane_to_append, cur_roi_idx } = [xv, yv];
    colorindex = colorindex+1;
    cur_roi_idx = cur_roi_idx + 1;
end

end

