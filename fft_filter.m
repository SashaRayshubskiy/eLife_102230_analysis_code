function [ out ] = fft_filter( data, cutoff, samplerate )

    a = fft( data );
    t = [0:1/samplerate:(size(data,1)-1)/samplerate];
    
    tcourse = [0:1/max(t):(1/mean(diff(t)))];
   
    bb = min(find(tcourse > cutoff));
       
    a(bb:end-bb,:) = 0;
        
    out = real(ifft(a));

    %figure; 
    %hold on;
    %plot(t, data);
    %plot(t, out, 'r--');
end