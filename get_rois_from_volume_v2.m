function [rois] = get_rois_from_volume_v2( asid, cdata, analysis_path )

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

PLANES = size(cdata,3);
IMAGE_ROWS = floor(sqrt(PLANES));
IMAGE_COLS = IMAGE_ROWS;

ac = get_analysis_constants;
order = ac.order;

rois_path = [analysis_path '/asid_' num2str(asid) '_rois.mat'];
if(exist( rois_path, 'file') )
    d = load(rois_path);
    rois = d.rois;
else

clear rois;
for i=1:PLANES
    
    nroi = 1;
    colorindex = 0;
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    subaxis(1,1,1, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
    
    % Throw away the first few volumes due to settling time.
    imagesc(mean(squeeze(cdata(:,:,i,3:end)),3));
    colormap gray;
    caxis([0 1500]);
    axis image;
IMAGE_ROWS = 4;
IMAGE_COLS = 4;
PLANES = IMAGE_ROWS * IMAGE_COLS;

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
    
        rois{ i, nroi } = [xv, yv];
        nroi = nroi + 1;
        colorindex = colorindex+1;
    end
    close(f);
end

save(rois_path, 'rois');
end

end

