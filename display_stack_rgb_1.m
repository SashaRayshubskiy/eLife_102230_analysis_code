function display_stack_rgb_1( stack, stack_path )

vid = VideoWriter([stack_path '_rgb_stack.avi']);
vid.FrameRate = 1;
vid.Quality = 100;
open(vid);

f = figure('units','normalized','outerposition',[0 0 1 1]);
for s = 3:size(stack,4)-2
%for s = 300
    
    im1 = squeeze(mean(stack( :, :, 2, s-2:s+2 ), 4 )); 
    im2 = squeeze(mean(stack( :, :, 1, s-2:s+2 ), 4 )); 
    
%     imagesc(im1)
%     hold on;
%     colormap gray(256)
%     overlayImg = imagesc(im2);
%         
%     aF = 0.5*ones(size(im2,1), size(im2,2));
%     
%     set(overlayImg, 'AlphaData', aF );
    
    %imagesc(imfuse(im1, im2));
    
    imagesc(im1+im2);
    colormap gray;
    colorbar;
    caxis([0 10000])
    axis image;
    title( ['frame - ' num2str( s ) ] );

    writeVideo(vid, getframe(f));
    pause(0.1);    
end

close(vid);

end

