function [ mask ] = get_dead_pixel_mask( ref_img )

[xsize, ysize] = size(ref_img);

THRESHOLD = 100;

mask = zeros(xsize, ysize);

mask( find(ref_img > THRESHOLD ) ) = 1;


end

