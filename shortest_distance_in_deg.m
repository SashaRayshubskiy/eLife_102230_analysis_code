function [ delta_deg ] = shortest_distance_in_deg( prev_deg, cur_deg )

delta_deg = rad2deg(angdiff(deg2rad( prev_deg ), deg2rad( cur_deg )));

end

