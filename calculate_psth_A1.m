function [ cur_psth ] = calculate_psth_A1( volts_t, vel_t, voltage, VOLTAGE_SR, SPIKE_THRESHOLD, BIN_SIZE )

USE_INST_FIRING_RATE = 0;
USE_EXP_FILTER = 0;
USE_HAMMING_FILTER = 1;

disp(['Got here: 2']);
%%%
% Get spikes from a voltage trace (current clamp)
%%%
VmFilt = medfilt1(voltage, 0.08 * VOLTAGE_SR, 'truncate');
delta_Vm = voltage - VmFilt;

%d = diff(delta_Vm);
%delta_Vm = d;

disp(['Got here: 1']);

volts_thresh = delta_Vm;
volts_thresh( find(delta_Vm < SPIKE_THRESHOLD) ) = 0;

if 0
    STEP_SIZE = 30 * 10000;
    START = 1;
    END   = STEP_SIZE;
    
    for ii = 1:10
        ff = figure;
        hold on;
        % plot(volts_t(1,2:end), d(1:30000));
        plot(volts_t(START:END), delta_Vm(START:END));
        % plot(volts_t(START+1:END+1), d(START:END));
        %plot(volts_t(START+1:END+1), volts_thresh(START+1:END+1));
        
        waitforbuttonpress;
        % close(ff);
        
        START = END+1;
        END   = END+STEP_SIZE;
    end
end

% [~, locs] = findpeaks(volts_thresh, 'MinPeakProminence', 0.02, 'Annotate','extents');
[~, locs] = findpeaks(volts_thresh, 'MinPeakDistance', 0.005 * VOLTAGE_SR );
%[~, locs] = findpeaks( volts_thresh );

disp(['Got here: 3']);

spikes = zeros(1, length(delta_Vm));
spikes(locs+1) = 1;

%%%%%%%%%%%%%%

if ( USE_INST_FIRING_RATE == 1 )
    % Calculate instantaneous firing rate
    cur_psth = zeros(size(spikes)); % initialize the array
    f = find(spikes); % spike times
    for c = 1:length(f)-1 % for all the spikes
        cur_psth(f(c):f(c+1)) = VOLTAGE_SR/(f(c+1)-f(c));
    end
    
    % Bin the firing rate to be on the same sampling rate as yaw
    cur_psth = squeeze( mean( reshape(cur_psth, [BIN_SIZE, length(voltage)/BIN_SIZE]),1) );    
else
    %%%
    % Get spikes from a voltage trace (current clamp)
    %%%    
    psth_dt = BIN_SIZE / (1.0*VOLTAGE_SR);
    disp(['Got here: 4']);
    sum_of_bins = sum(reshape(spikes, [BIN_SIZE,  length(voltage)/BIN_SIZE]),1);
    
    cur_psth =  sum_of_bins./ psth_dt;
    
    if ( USE_EXP_FILTER == 1 )        
        cur_psth = smoothts(cur_psth, 'e', 10);
    elseif ( USE_HAMMING_FILTER == 1 )
        cur_psth = hanningsmooth(cur_psth, 40);
    end
    disp(['Got here: 5']);
    %%%%%%%%%%%%%%
end

if 0
    STEP_SIZE = 30 * 10000;
    START = 1;
    END   = STEP_SIZE;

    STEP_SIZE_v = 30 * 100;
    START_v = 1;
    END_v   = STEP_SIZE_v;
    
    for ii = 1:10
        
        ff = figure;
        subplot(2,1,1);
        
        hold on;
        %     plot(volts_t{1}(1,:), voltage );
        %     plot(volts_t{1}(1,:), voltage .* spikes, 'x', 'color', 'g');
        plot(volts_t(START:END), voltage(START:END) );
        plot(volts_t(START:END), voltage(START:END) .* spikes(START:END)', 'x', 'color', 'g');
        
        subplot(2,1,2);
        %    plot(vel_t{1}( 1, : ), cur_psth);
        plot( vel_t( START_v : END_v ), cur_psth( START_v : END_v ) );
        
        waitforbuttonpress;

        START = END+1;
        END   = END+STEP_SIZE;

        START_v = END_v + 1;
        END_v   = END_v + STEP_SIZE_v;
        % close(ff);
    end
    %break;
end
end