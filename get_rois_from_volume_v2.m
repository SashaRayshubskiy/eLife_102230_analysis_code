function [rois_per_plane] = get_rois_from_volume_v2( cdata )

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

IMAGE_ROWS = 4;
IMAGE_COLS = 4;
PLANES = IMAGE_ROWS * IMAGE_COLS;

if 0
for i=1:PLANES
    all_axes(i) = subaxis(IMAGE_ROWS, IMAGE_COLS, i, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
    
    % Throw away the first few volumes due to settling time.
    imagesc(mean(squeeze(cdata(:,:,i,3:end)),3));
    colormap gray;
    caxis([0 3000]);
    axis image;
    axis off;
end
end

order    = [ rgb('Blue'); rgb('Green'); rgb('Red'); rgb('Black'); rgb('Purple'); rgb('Brown'); rgb('Indigo'); rgb('DarkRed') ];
clear rois_per_plane;
for i=1:PLANES
    
    nroi = 1;
    colorindex = 0;
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    subaxis(1,1,1, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
    
    % Throw away the first few volumes due to settling time.
    imagesc(mean(squeeze(cdata(:,:,i,3:end)),3));
    colormap gray;
    caxis([0 3000]);
    axis image;
    axis off;
    
    ti = title(['Plane: ' num2str(i)]);
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
    close(f);
end

end

