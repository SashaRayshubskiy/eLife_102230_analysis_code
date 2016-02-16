function [ img_exp ] = expand_img( img, dx, dy )

[xsize, ysize] = size(img);

img_exp = zeros(xsize*dx, ysize*dy);

for r = 1:xsize
    for c = 1:ysize
        cur_val = img(r,c);
        for rr = 1:dx
            for cc = 1:dy
                cur_rr = (r-1)*dx + rr;
                cur_cc = (c-1)*dy + cc;
                
                img_exp(cur_rr, cur_cc) = cur_val;
            end
        end
    end
end

end

