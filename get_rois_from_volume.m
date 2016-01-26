function [rois_per_plane] = get_rois_from_volume( cdata )

ac = get_analysis_constants;
order = ac.order;

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

IMAGE_ROWS = 4;
IMAGE_COLS = 4;
PLANES = IMAGE_ROWS * IMAGE_COLS;

f = figure('units','normalized','outerposition',[0 0 1 1]);
all_axes = [];

for i=1:PLANES
    all_axes(i) = subaxis(IMAGE_ROWS, IMAGE_COLS, i, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
    
    % Throw away the first few volumes due to settling time.
    imagesc(mean(squeeze(cdata(:,:,i,3:end)),3));
    colormap gray;
    caxis([0 3000]);
    axis image;
    axis off;
end

clear rois_per_plane;
for i=1:PLANES
    
    nroi = 1;
    colorindex = 0;
    
    axes(all_axes(i));
    ti = title('Active Plane');
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
    
        rois_per_plane{ i, nroi } = [xv, yv];
        nroi = nroi + 1;
        colorindex = colorindex+1;
    end   
    title('');
end

end

