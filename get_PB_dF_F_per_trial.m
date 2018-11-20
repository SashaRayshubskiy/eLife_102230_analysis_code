function glom_dF_F_per_trial = get_PB_dF_F_per_trial( F_per_trial )

% F_per_trial: { trials, glomeruli, time }

glom_dF_F_per_trial = zeros( size(F_per_trial,1), size(F_per_trial,2), size(F_per_trial,3) );



for g = 1:size(F_per_trial,2)    
    cur_glom = squeeze( F_per_trial( :, g, : ) );
    % cur_baseline = mean( cur_glom(:) ) - 2*std( cur_glom(:), 1 );

    % Use the bottom 5% as baseline
    cur_glom_sorted = sort( cur_glom(:) );
    %five_percent_index = ceil( length(cur_glom(:)) * 0.05 );
    five_percent_index = ceil( length(cur_glom(:)) * 0.5 );
    cur_baseline = mean(cur_glom_sorted( 1:five_percent_index ));

    % cur_baseline = mean(cur_glom(:));
        
    for tr = 1:size( F_per_trial )

        cur_F = squeeze(F_per_trial( tr, g, : ));

        % cur_baseline = mean(cur_F);

        glom_dF_F_per_trial( tr, g, :) = (cur_F - cur_baseline) ./ cur_baseline;    
    end
end
end

