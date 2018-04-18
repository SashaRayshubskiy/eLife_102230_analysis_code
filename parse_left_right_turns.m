function [ left_yaw_ids, right_yaw_ids ] = parse_left_right_turns( t_down, yaw_all_down )

figure;
plot(t_down, yaw_all_down);
xlim([0 t_down(end)]);
end

