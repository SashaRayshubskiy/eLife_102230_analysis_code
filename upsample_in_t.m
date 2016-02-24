function [ mov_up ] = upsample_in_t( mov, n )

[xsize, ysize, nframes] = size(mov);

new_nframes =(nframes-1) * n;

mov_up = zeros(xsize, ysize, new_nframes);

org_frames = [1:nframes];
new_frames = [1:1/n:nframes];

for x = 1:xsize
    for y = 1:ysize
        ss = spline( org_frames, squeeze(mov(x,y,:)), new_frames(1:end-1) );
        mov_up(x,y,:) = ss;
    end
end

end

