function [ rois ] = get_rois_from_volume( PLANE, cdata, rois_path )

ac = get_analysis_constants;
order = ac.order;

if(exist( rois_path, 'file') )
    d = load(rois_path);
    rois = d.rois;
else
    clear rois;
    
    nroi = 1;
    colorindex = 0;
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    
    % Throw away the first few volumes due to settling time.
    imagesc(mean(squeeze(cdata(:,:,PLANE,3:end)),3));
    colormap gray;
    caxis([0 1000]);
    axis image;
    axis off;
    
    ti = title(['Plane: ' num2str(PLANE)]);
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
        
        rois{ nroi } = [xv, yv];
        nroi = nroi + 1;
        colorindex = colorindex+1;
    end
    close(f);
    
    save(rois_path, 'rois');
end

end

