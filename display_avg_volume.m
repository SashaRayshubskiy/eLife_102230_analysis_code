function display_avg_volume( down_params, cdata, figsave_prefix )

SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

x_size = size(cdata,1);
y_size = size(cdata,2);
PLANES = size(cdata,3);
IMAGE_ROWS = floor(sqrt(PLANES));
IMAGE_COLS = IMAGE_ROWS;

f = figure('Color', 'none');
for i=1:PLANES
    subaxis(IMAGE_ROWS, IMAGE_COLS, i, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
    
    % Throw away the first few volumes due to settling time.
    imagesc([1:y_size*down_params(2)], [1:x_size*down_params(1)], mean(squeeze(cdata(:,:,i,3:end)),3));
    colormap gray;
    caxis([0 500]);
    axis image;
    axis off;
end

saveas(f, [figsave_prefix '_mean_slices_in_volume.fig']);
saveas(f, [figsave_prefix '_mean_slices_in_volume.png']);

end

