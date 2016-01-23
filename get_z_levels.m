function [ z_level ] = get_z_levels( tdata, XDIM, PLANES, VOLUMES )

x = [1:XDIM];
z_level = zeros(PLANES, VOLUMES);
for v=1:VOLUMES
    for p=1:PLANES
        y = squeeze(mean(squeeze(tdata(:,:,p,v)),2));
        f = fit(x',y,'gauss2');
        z_level( p, v ) = f.b1;
    end
end

end

