function [ sem ] = get_sem( data, dim )

    sem = std( data, dim ) ./ sqrt( size( data, dim ));
    
end

