function [yaw_deg] = convert_yaw_to_degrees( yaw_au )

A_YAW = 0.000395;
B_YAW = 0.0034;

yaw_deg = (yaw_au - B_YAW) ./ A_YAW;
end

