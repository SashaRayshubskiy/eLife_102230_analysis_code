function [ out_base, out_mean, out_95ci ] = get_95_ci_for_spont_vs_cx( FR_data, yaw_data )

% Assumes that the data is already shifted

MAX_FR = 160;
BIN_FACTOR = 20;
BIN_COUNT = MAX_FR / BIN_FACTOR;

yaw_bin_on_FR = cell( 1, BIN_COUNT );

for i = 1:BIN_COUNT
    yaw_bin_on_FR{i} = [];
end

for i = 1:length( FR_data )
    cur_FR  = FR_data( i );
    cur_yaw = yaw_data( i );
    
    FR_bin = floor(cur_FR/BIN_FACTOR) + 1;
    if(FR_bin > BIN_COUNT)
        FR_bin = BIN_COUNT;
    end
    
    yaw_bin_on_FR{ FR_bin }(end+1) = cur_yaw;
end

OCCUPANCY_LIMIT = 2;
for i = 1:BIN_COUNT
    
    out_base( i ) = (i-1) * BIN_FACTOR;
    
    if( length( yaw_bin_on_FR{ i } ) > OCCUPANCY_LIMIT )
        out_mean( i ) = mean(  yaw_bin_on_FR{ i } );
        out_95ci( i ) = 2 * std( yaw_bin_on_FR{ i }, 1 );
    else
        out_mean( i ) = NaN;
        out_95ci( i ) = NaN;
    end
end


