function [ mov_down ] = downsample_with_mask( mov, mask, dx, dy )

[ xsize, ysize, nframes ] = size( mov );

d_sizex = xsize/dx;
d_sizey = ysize/dy;

mov_down = zeros( d_sizex, d_sizey, nframes, 'double' );

pixels_in_block = dx*dy;

for r_b = 1:d_sizex
    for c_b = 1:d_sizey
        used_pixels = 0.0;
        cur_intens = zeros(nframes,1);
        for rr_b = 1:dx
            for cc_b = 1:dy
                r = (r_b-1)*dx + rr_b;
                c = (c_b-1)*dy + cc_b;
                
                cur_pixel_mask = mask(r,c);
                
                if( cur_pixel_mask == 1 )
                    cur_intens = cur_intens + squeeze(mov(r,c,:));
                    used_pixels = used_pixels + 1.0;
                end
            end
        end
        
        if( (1.0*used_pixels)/(1.0*pixels_in_block) > .75 )
            mov_down(r_b,c_b,:) = cur_intens./used_pixels;
        end
    end
end

end

