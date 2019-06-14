function [ roi ] = get_PFLm_rois( imaging_data, analysis_path )

% imaging_data = { trials, x, y, planes, time }
% glom is a cell array of { planes, roi }, there will always be 8 glomeruli

CAXIS_MAX = 250;

roi_mat_file = [analysis_path '/rois.mat'];
if( exist(roi_mat_file,'file') > 0 )
    d = load(roi_mat_file);
    roi = d.roi;
else
    
    f = figure;    
    num_planes = size( imaging_data, 4 );
    for p = 1:num_planes
        subplot(3,4,p);

        cur_ref_img = squeeze(mean(mean( squeeze( imaging_data(:,:,:,p,:)), 4 ),1));
        imshow( cur_ref_img, [], 'InitialMagnification', 'fit' );
        caxis( [0 CAXIS_MAX] );
        title(['Plane: ' num2str(p)]);
    end
    
    saveas(f, [analysis_path '/all_active_planes.fig']);
    saveas(f, [analysis_path '/all_active_planes.png']);
    
    disp('Got here 1');
%     waitforbuttonpress;
        
    % { left PFL.LAL, right PFL.LAL, 8 FB glomeruli starting left to right }
    NUM_roi = 12;    
    roi = cell(NUM_roi,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 190504_VT7338_lexA_LexAOp_GCaMP7f_01
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     roi{1,1} = [2:3];
    %     roi{2,1} = [6:8];
    %     roi{3,1} = [11:12];
    %     roi{4,1} = [9:11];
    %     roi{5,1} = [7:10];
    %     roi{6,1} = [6:8];
    %     roi{7,1} = [6:8];
    %     roi{8,1} = [7:10];
    %     roi{9,1} = [9:11];
    %     roi{10,1} = [11:12];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 190506_VT7338_lexA_LexAOp_GCaMP7f_02
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    roi{1,1} = [1:7];     % Left PFL.LAL axon 
    roi{2,1} = [1:8];    % Left PFL-FB - PFL.LAL neurite
    roi{3,1} = [6:9];     % Right PFL.LAL axon 
    roi{4,1} = [1:8];     % Right PFL-FB - PFL.LAL neurite
    
    roi{5,1} = [11:12];
    roi{6,1} = [11:12];
    roi{7,1} = [11:12];
    roi{8,1} = [11:12];
    roi{9,1} = [11:12];
    roi{10,1} = [11:12];
    roi{11,1} = [11:12];
    roi{12,1} = [11:12];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    order_avg    = [ rgb('Blue'); rgb('Green'); rgb('Red'); rgb('Black'); rgb('Purple'); rgb('Brown'); rgb('Indigo'); rgb('DarkRed') ];
    order_single = [ rgb('LightBlue'); rgb('LightGreen'); rgb('LightSalmon'); rgb('LightGray'); rgb('Violet'); rgb('Bisque'); rgb('Plum'); rgb('LightPink') ];
    nroi = 1;
    npts = 1;
    colorindex = 0;

    my_figs = cell(1,NUM_roi);    
    
    for g = 1:NUM_roi
        my_figs{g} = figure;    
        cur_planes = roi{g,1};
        
        % subplot(5,2,g);
        hold on;
        cur_ref_img = squeeze(mean(squeeze(mean( squeeze(max(imaging_data(:,:,:,cur_planes,:), [], 4)), 1)), 3));
        imshow( cur_ref_img, [], 'InitialMagnification', 'fit' );
        caxis( [0 CAXIS_MAX] );
        title(['ROIs: ' num2str(g) ' planes [ ' num2str( cur_planes ) ']']);        
    end
        
    for g = 1:NUM_roi
        % subplot(4,2,g);
        f = figure(my_figs{g});
        [xv, yv] = (getline(gca, 'closed'));
        
        currcolor_avg = order_avg(1+mod(colorindex,size(order_avg,1)),:);
        plot(xv, yv, 'Linewidth', 1,'Color',currcolor_avg);
        text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor_avg,'FontSize',12);
        colorindex = colorindex+1;
        roi{g,2} = [xv, yv];
        saveas(f, [analysis_path '/rois_' num2str(g) '.fig']);
        saveas(f, [analysis_path '/rois_' num2str(g) '.png']);
    end
    
    save( roi_mat_file, 'roi' );
end

end

