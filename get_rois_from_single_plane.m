function [rois] = get_rois_from_single_plane( asid, cdata, analysis_path )

ac = get_analysis_constants;
order = ac.order;

rois_path = [analysis_path '/asid_single_' num2str(asid) '_rois.mat'];
if(exist( rois_path, 'file') )
    d = load(rois_path);
    rois = d.rois;
else
    
    clear rois;
    
    nroi = 1;
    colorindex = 0;
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    
    % Throw away the first few volumes due to settling time.
    imagesc(mean(squeeze(cdata(:,:,3:end)),3));
    colormap gray;
    caxis([0 1000]);
    axis image;
    axis off;
    
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
        
        rois{ 1, nroi } = [xv, yv];
        nroi = nroi + 1;
        colorindex = colorindex+1;
    end
    close(f);

    save(rois_path, 'rois');
end

end

