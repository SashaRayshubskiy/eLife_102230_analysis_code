function filtsig = hanningsmooth_circ(signal,win)
% Filter signal with a hanning window of width win
% Cuts off ends of convolved signal to align to original signal


half_win = win/2;

signal_padded = [ signal(end-half_win:end); signal; signal(1:half_win) ];

smooth_win = hanning(win)/sum(hanning(win));
filtsig = conv(signal_padded,smooth_win);

sigrange = [win:win+length(signal)-1];
filtsig = filtsig(sigrange);

assert( length(signal) == length( filtsig) );

%figure; plot(signal); hold on; plot(filtsig,'r');
