function [ stack ] = generate_stack_movie( stack_path )

%stacks = open_tif_fast_simple( stack_path, 256, 512, 800, 1, 40, 40*800 )
stacks = squeeze(open_tif_fast(stack_path));

stack = squeeze(mean(stacks,4));

end


