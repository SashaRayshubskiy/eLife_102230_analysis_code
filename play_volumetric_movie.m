function play_volumetric_movie( cdata, VPS, figsave_prefix )

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

IMAGE_ROWS = 4;
IMAGE_COLS = 4;
PLANES = IMAGE_ROWS * IMAGE_COLS;


%vidObj = VideoWriter([figsave_prefix '.avi']);
%vidObj.Quality = 100;
%vidObj.FrameRate = VPS;
%open(vidObj);

tmp_output_dir = '/tmp/video_output_dir/';
if(~exist(tmp_output_dir, 'dir'))
    mkdir(tmp_output_dir);
end

f = figure('units','normalized','outerposition',[0 0 1 1]);
for v=1:size(cdata,4)
    for i=1:PLANES
        subaxis(IMAGE_ROWS, IMAGE_COLS, i, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
        
        % Throw away the first few volumes due to settling time.
        imagesc(squeeze(cdata(:,:,i,v)));
        colormap gray;
        caxis([0 2000]);
        axis image;
        axis off;
    end
    
    subaxis(IMAGE_ROWS, IMAGE_COLS, 1, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
    title(['Time: ' num2str(v*1.0/VPS)]);
    
    %saveas(f, [tmp_output_dir '/tmp_img_' sprintf('%05d', v) '.png']);
    pause(0.05);    
end

%system(['ffmpeg -r 1 -f image2 -i ' tmp_output_dir '/tmp_img_%05d.png  ' figsave_prefix '.avi']);

% rmdir(tmp_output_dir);

end

