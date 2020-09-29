function [fwd_mms] = convert_fwd_to_mms( fwd_au )

A_FWD = 0.00636;
B_FWD = 0.0000;

fwd_mms = (fwd_au - B_FWD) ./ A_FWD;
end

