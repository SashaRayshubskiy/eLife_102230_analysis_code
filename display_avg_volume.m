function display_avg_volume( cdata, figsave_prefix )

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

IMAGE_ROWS = 4;
IMAGE_COLS = 4;
PLANES = IMAGE_ROWS * IMAGE_COLS;

f = figure;
for i=1:PLANES
    subaxis(IMAGE_ROWS, IMAGE_COLS, i, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
    
    % Throw away the first few volumes due to settling time.
    imagesc(mean(squeeze(cdata(:,:,i,3:end)),3));
    colormap gray;
    caxis([0 2000]);
    axis image;
    axis off;
end

saveas(f, [figsave_prefix '_mean_slices_in_volume.fig']);
saveas(f, [figsave_prefix '_mean_slices_in_volume.png']);

end

