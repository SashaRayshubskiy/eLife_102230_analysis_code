function  bump_motion_ids = get_weigh_avg_bump_pos( cur_bump_gloms )

% assumes input { glomerulus, time }

n_gloms  = size( cur_bump_gloms, 1 );
n_frames = size( cur_bump_gloms, 2 );

bump_motion_ids = zeros(1, n_frames);

LARGE_FACTOR = 10000;
x = [1:n_gloms] + LARGE_FACTOR;

for f = 1:n_frames
    cur_gloms = squeeze(cur_bump_gloms(:,f));
    sum_cur_gloms = sum(cur_gloms);
    xx = x .* cur_gloms';
    xx_sum = sum(xx);
   
    cur_bump_pos = (xx_sum / sum_cur_gloms) - LARGE_FACTOR;
    
    bump_motion_ids(f) = cur_bump_pos;
    
end

end

