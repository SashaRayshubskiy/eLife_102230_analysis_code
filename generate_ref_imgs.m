function [ ref_imgs ] = generate_ref_imgs( cdata_raw )

xsize = size(cdata_raw{1}, 2);
ysize = size(cdata_raw{1}, 3);
PLANES = size(cdata_raw{1}, 4);
VOLS = size(cdata_raw{1}, 5);
IMAGES = VOLS;

ref_imgs = zeros( PLANES, xsize, ysize );
for p = 1:PLANES
    ref_imgs(p,:,:) = mean(squeeze(cdata_raw{ 1 }(1,:,:,p,1:IMAGES)),3);
end
end

