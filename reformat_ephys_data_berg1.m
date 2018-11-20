function [ ephys_time, ephys_data ] = reformat_ephys_data_berg1( bdata_time_in, bdata_raw )

ephys_time = bdata_time_in;

trial_type_cnt = size(bdata_raw, 2);

for i=1:trial_type_cnt
    trial_per_type_cnt = size(bdata_raw{i},1);
           
    ephys_data{ i } = zeros(trial_per_type_cnt, length(ephys_time), 'double' );
    
    for j = 1:trial_per_type_cnt    
        [ current, voltage ] = get_scaled_voltage_and_current_A_berg1( squeeze(bdata_raw{i}(j, :,:)) );
        ephys_data{i}(j, :) = voltage;    
    end
end

end

