function filtsig = hanningsmooth(signal,win);
% Filter signal with a hanning window of width win
% Cuts off ends of convolved signal to align to original signal


smooth_win = hanning(win)/sum(hanning(win));
filtsig = conv(signal,smooth_win);

sigrange = [win/2:win/2+length(signal)-1];
filtsig = filtsig(sigrange);

%figure; plot(signal); hold on; plot(filtsig,'r');
