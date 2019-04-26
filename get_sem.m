function [ sem ] = get_sem( data, dim )

    sem = std( data, 1, dim ) ./ sqrt( size( data, dim ));
    
end

