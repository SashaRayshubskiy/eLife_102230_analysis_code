function display_stack( channel, stack, stack_path )

vid = VideoWriter([stack_path '_chan_' num2str(channel) '.avi']);
vid.FrameRate = 1;
vid.Quality = 100;
open(vid);

f = figure('units','normalized','outerposition',[0 0 1 1]);
for s = 3:size(stack,4)-2
%for s = 300
    
    imagesc(squeeze(mean(stack( :, :, channel, s-2:s+2 ), 4 )));
    axis image;
    colormap gray;
    caxis([0 5500]);
    title( ['chan(' num2str( channel ) '): frame - ' num2str( s ) ] );

    writeVideo(vid, getframe(f));
end

close(vid);

end

