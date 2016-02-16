function display_stack( stack, stack_path )

vid = VideoWriter([stack_path '.avi']);
vid.FrameRate = 1;
vid.Quality = 100;
open(vid);

f = figure('units','normalized','outerposition',[0 0 1 1]);
for s = 3:size(stack,3)-2
    imagesc(squeeze(mean(stack(:,:,s-2:s+2),3)));
    colormap gray;
    caxis([0 5500]);

    writeVideo(vid, getframe(f));
end

close(vid);

end

