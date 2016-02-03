function [ img_filt ] = filter_dithered_image( img )

[ x_size, y_size ] = size(img);
 
img_filt = zeros(x_size, y_size);

PIX_BLOCK = 5;

total_pixel_cnt_in_block = PIX_BLOCK * PIX_BLOCK;

for i = 3:x_size-2
    for j = 3:y_size-2

       img_filt(i,j) =  get_filt_dithered_pixel_value_v2( img(i-2:i+2, j-2:j+2), total_pixel_cnt_in_block );

    end
end

end

function [ new_pixel_val ] = get_filt_dithered_pixel_value_v2( pb, pb_cnt )
% Assumes a 5x5 pixel block

center_pixel = pb(3,3);

if(center_pixel>0)
    new_pixel_val = center_pixel;
else
    new_pixel_val = 0;
end
end

function [ new_pixel_val ] = get_filt_dithered_pixel_value( pb, pb_cnt )
% Assumes a 5x5 pixel block

center_pixel = pb(3,3);

num_pos = sum(sum(find(pb > 0)));
per_pos_pixels = 1.0*num_pos ./ pb_cnt;

% If the majority of pixels in block are positive, then this negative
% pixel is likely noise and should be an average of it's neighbors. And
% vice versa.
if( per_pos_pixels > 0.5 )
    if( center_pixel > 0 )
        new_pixel_val = center_pixel;
    else
        new_pixel_val = mean(mean(pb));
    end
else
    if( center_pixel > 0 )
        new_pixel_val = mean(mean(pb));
    else
        new_pixel_val = center_pixel;
    end
end

end