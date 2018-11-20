function [glom_F_per_trial] = get_PB_F_per_trial( glom, F_data )

% This function handles planes to glomeruli accumulation

% F_data 
% { trials, x, y, planes, time }

% glom is a cell array of { planes, roi }, there will always be 8 glomeruli

num_gloms = length(glom);
num_trials = size(F_data,1);
nframes    = size(F_data,5);
ysize = size( F_data, 2 );
xsize = size( F_data, 3 );

[x, y] = meshgrid(1:xsize, 1:ysize);

glom_F_per_trial = zeros( num_trials, num_gloms, nframes );

for g = 1 : num_gloms
    
    cur_planes = glom{ g, 1 }; 
    cur_roi    = glom{ g, 2 };        
    
    xv = cur_roi(:,1);
    yv = cur_roi(:,2);
            
    inpoly = inpolygon(x,y,xv,yv);
    
    cur_F_data = squeeze(mean( F_data(:,:,:,cur_planes,:), 4 ));
    
    for tr = 1 : num_trials 
                
        cur_F_tc_in_roi = squeeze(sum(sum(squeeze(cur_F_data(tr,:,:,:)).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));  
    
        glom_F_per_trial(tr,g, :) = cur_F_tc_in_roi;
    end    
end

cc = 32;

end

