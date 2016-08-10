function [ seq ] = convert_bmeta_to_sequential_data( bmeta )

cur_tt_1_idx = 1;
cur_tt_2_idx = 1;

seq = [];

while ((cur_tt_1_idx <= length(bmeta{1})) | (cur_tt_2_idx <= length(bmeta{2})))
    
    if (cur_tt_1_idx <= length(bmeta{1}) )        
        cur_tid_1 = bmeta{1}( cur_tt_1_idx, 2 );
    else
        cur_tid_1 = 10000;
    end
    
    if (cur_tt_2_idx <= length(bmeta{2}) )        
        cur_tid_2 = bmeta{2}( cur_tt_2_idx, 2 );
    else
        cur_tid_2 = 10000;
    end
    
    if(  cur_tid_1 < cur_tid_2 )
        seq( end+1, : ) = [ 1 cur_tt_1_idx ];
        cur_tt_1_idx = cur_tt_1_idx + 1;
    else
        seq( end+1, : ) = [ 2 cur_tt_2_idx ];
        cur_tt_2_idx = cur_tt_2_idx + 1;        
    end

end

end

