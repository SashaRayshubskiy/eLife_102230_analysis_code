function [ pico_stim_data ] = reformat_pico_stim_data_berg1( bdata_raw )

settings = sensor_settings;

trial_type_cnt = size(bdata_raw, 2);
samples_per_trial        = size(bdata_raw{1}, 2);

for i=1:trial_type_cnt
    trial_per_type_cnt = size(bdata_raw{i},1);
           
    pico_stim_data{ i } = zeros(trial_per_type_cnt, samples_per_trial, 'double' );
    
    for j = 1:trial_per_type_cnt
        pico_stim_data{i}(j, :) = squeeze(bdata_raw{i}(j, :, settings.BERG1_PICO_MONITOR )) ;    
    end
end

end

