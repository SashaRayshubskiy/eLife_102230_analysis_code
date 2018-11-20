function [glom] = get_left_half_PN_glomeruli_10202018( EPG_data, analysis_path )

% EPG_data = { trials, x, y, planes, time }
% glom is a cell array of { planes, roi }, there will always be 8 glomeruli

CAXIS_MAX = 100;

glom_mat_file = [analysis_path '/glomeruli.mat'];
if( exist(glom_mat_file,'file') > 0 )
    d = load(glom_mat_file);
    glom = d.glom;
else
    
    f = figure;    
    num_planes = size( EPG_data, 4 );
    for p = 1:num_planes
        subplot(3,2,p);

        cur_ref_img = squeeze(mean(mean( squeeze(EPG_data(:,:,:,p,:)), 4 ),1));
        imshow( cur_ref_img, [], 'InitialMagnification', 'fit' );
        caxis( [0 CAXIS_MAX] );
        title(['Plane: ' num2str(p)]);
    end
    
    saveas(f, [analysis_path '/all_active_planes.fig']);
    saveas(f, [analysis_path '/all_active_planes.png']);
    
    
    NUM_GLOM = 8;
    glom = cell(NUM_GLOM,2);
    
    % 181019_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_13
%     glom{1,1} = 5; glom{2,1} = 4; glom{3,1} = [1:2];  glom{4,1} = [1:2];
%     glom{5,1} = [1:2]; glom{6,1} = [1:2]; glom{7,1} = [1:5];  glom{8,1} = [3:5];
    
    % 181019_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_12
%     glom{1,1} = 6; glom{2,1} = [4:5]; glom{3,1} = [1:3];  glom{4,1} = [1:3];
%     glom{5,1} = [2:4]; glom{6,1} = [3:4]; glom{7,1} = [4:6];  glom{8,1} = [4:6];
    
    % 181022_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_14
%     glom{1,1} = 6; glom{2,1} = [5]; glom{3,1} = [1:2];  glom{4,1} = [1:2];
%     glom{5,1} = [1:3]; glom{6,1} = [1:3]; glom{7,1} = [4:5];  glom{8,1} = [4:5];

    % 181022_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_15 - EB
     glom{1,1} = [1:3]; glom{2,1} = [1:3]; glom{3,1} = [1:3];  glom{4,1} = [4];
     glom{5,1} = [5:6]; glom{6,1} = [5:6]; glom{7,1} = [5:6];  glom{8,1} = [4];
    
    f = figure;    
    
    for g = 1:NUM_GLOM
        cur_planes = glom{g,1};
        
        subplot(4,2,g);
        hold on;
        cur_ref_img = squeeze(mean(squeeze(mean( squeeze(max(EPG_data(:,:,:,cur_planes,:), [], 4)), 1)), 3));
        imshow( cur_ref_img, [], 'InitialMagnification', 'fit' );
        caxis( [0 CAXIS_MAX] );
        title(['Glomerulus: ' num2str(g) ' planes [ ' num2str( cur_planes ) ']']);        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % Specify plane for glomerulus 1
    % 
    % This can be fly specific.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%     glom{1,1} = 6;
%     glom{2,1} = [4:5];
%     for g = 3:8
%         glom{g,1} = [1:3];
%     end
    
    % Fly 181015_Lex_6f_60D05_Gal4_P2X2_PEN1_recomb_09
%    glom{1,1} = 6; glom{2,1} = [4:5]; glom{3,1} = 3;  glom{4,1} = [1:3]; 
%    glom{5,1} = [1:3]; glom{6,1} = [1:3]; glom{7,1} = [3:5];  glom{8,1} = [3:5];     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    hold on;
    
    order_avg    = [ rgb('Blue'); rgb('Green'); rgb('Red'); rgb('Black'); rgb('Purple'); rgb('Brown'); rgb('Indigo'); rgb('DarkRed') ];
    order_single = [ rgb('LightBlue'); rgb('LightGreen'); rgb('LightSalmon'); rgb('LightGray'); rgb('Violet'); rgb('Bisque'); rgb('Plum'); rgb('LightPink') ];
    nroi = 1;
    npts = 1;
    colorindex = 0;
        
    for g = 1:NUM_GLOM
        subplot(4,2,g);
        [xv, yv] = (getline(gca, 'closed'));
        
        currcolor_avg = order_avg(1+mod(colorindex,size(order_single,1)),:);
        plot(xv, yv, 'Linewidth', 1,'Color',currcolor_avg);
        text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor_avg,'FontSize',12);
        colorindex = colorindex+1;
        glom{g,2} = [xv, yv];
    end
    
    save( glom_mat_file, 'glom' );
    
    saveas(f, [analysis_path '/glomeruli_rois.fig']);
    saveas(f, [analysis_path '/glomeruli_rois.png']);
end

end

