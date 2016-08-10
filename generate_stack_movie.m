function [ stack ] = generate_stack_movie( stack_path )

%stacks = open_tif_fast_simple( stack_path, 256, 512, 700, 1, 22, 22*700*1 );
stacks = squeeze(open_tif_stack(stack_path));

stack = squeeze(mean(stacks,5));

end


