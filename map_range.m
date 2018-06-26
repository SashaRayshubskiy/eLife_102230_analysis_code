function [ to_val ] = map_range( from_range, to_range, from_val )

    if( from_val < from_range(1) )
        to_val = to_range(1);
    elseif( from_val > from_range(2) )
        to_val = to_range(2);
    else
        to_val = ceil(interp1( from_range, to_range, from_val ));
    end
    
end

