function [ cur_psth ] = calculate_psth_bak(voltage, VOLTAGE_SR, SPIKE_THRESHOLD, BIN_SIZE);

psth_dt = BIN_SIZE / (1.0*VOLTAGE_SR);

%%%
% Get spikes from a voltage trace (current clamp)
%%%
cur_volt = voltage;
VmFilt = medfilt1(cur_volt, 0.08 * VOLTAGE_SR, 'truncate');
delta_Vm = cur_volt - VmFilt;

d = diff(delta_Vm);

volts_thresh = d;
volts_thresh(find(d < SPIKE_THRESHOLD)) = 0;

%[~, locs] = findpeaks(volts_thresh, 'MinPeakProminence',0.02, 'Annotate','extents');
[~, locs] = findpeaks(volts_thresh, 'Annotate','extents');

spikes = zeros(1, length(delta_Vm));
spikes(locs+1) = 1;
%%%%%%%%%%%%%%

figure;
hold on;
plot(d);
plot(d(locs),  'gx')



%%%
% Get spikes from a voltage trace (current clamp)
%%%
cur_psth = sum(reshape(spikes, [BIN_SIZE, size(voltage,1)/BIN_SIZE]),1) / psth_dt;
cur_psth = hanningsmooth(cur_psth, 50);
%%%%%%%%%%%%%%

end


        