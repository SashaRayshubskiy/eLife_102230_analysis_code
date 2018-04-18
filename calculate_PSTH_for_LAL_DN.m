function [ volt_PSTH ] = calculate_PSTH_for_LAL_DN( t_volts, t_vel_all, volts_all, DN_TYPE )

trial_type_cnt = length(volts_all);

volt_PSTH = cell(trial_type_cnt,1);

settings = sensor_settings;
ephys_SR = settings.sampRate;
ball_SR = settings.sensorPollFreq;

psth_dt_samples = ephys_SR/ball_SR;
psth_dt = psth_dt_samples / (1.0*ephys_SR);

SPIKE_THRESHOLD_LAL_DN = 0.5; % Good for A2
if( strcmp(DN_TYPE, 'A2') == 1 ) 
    SPIKE_THRESHOLD_LAL_DN = 0.5; % Good for A2
elseif( strcmp(DN_TYPE, 'A1') == 1 ) 
    SPIKE_THRESHOLD_LAL_DN = 0.4; % Good for A1
end

for tt=1:trial_type_cnt
    
    volt_PSTH{tt} = zeros(size(volts_all{ tt },1), size(volts_all{ tt },2)/psth_dt_samples); 
    
    for tr=1:size(volts_all{ tt },1)
        
        %%%
        % Get spikes from a voltage trace (current clamp)
        %%%                
        cur_volt = squeeze(volts_all{ tt }( tr, : ));
        
        cur_psth = calculate_psth( t_volts, t_vel_all, cur_volt, ephys_SR, SPIKE_THRESHOLD_LAL_DN, psth_dt_samples );
        
        volt_PSTH{tt}(tr,:) = cur_psth;
    end
    
end
end