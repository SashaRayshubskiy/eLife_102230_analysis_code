function display_stack_one_chan( stack, stack_path )

vid = VideoWriter([stack_path '.avi']);
vid.FrameRate = 1;
vid.Quality = 100;
open(vid);

f = figure('units','normalized','outerposition',[0 0 1 1]);
for s = 3:size(stack,3)-2
    
    imagesc(squeeze(flipud(stack( :, :, s ))));
    axis image;
    colormap gray;
    caxis([0 7000]);
    title( [' frame - ' num2str( s ) ] );

    writeVideo(vid, getframe(f));
    pause(0.1);    

end

close(vid);

end

