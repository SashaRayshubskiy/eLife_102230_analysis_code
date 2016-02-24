function [ out ] = fft_filter_3D( data, cutoff, samplerate )

[xsize, ysize, nframes] = size(data);

out = zeros( xsize, ysize, nframes );

for x = 1:xsize
    for y = 1:ysize
        out(x,y,:) = fft_filter( squeeze(data(x,y,:)), cutoff, samplerate );
    end
end

end

