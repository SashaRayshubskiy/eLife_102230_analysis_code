function [ df_f_in_roi ] = get_df_f_in_roi( trial_cdata, base_begin, base_end, roi )

x_size = size(trial_cdata, 1);
y_size = size(trial_cdata, 2);
nframes = size(trial_cdata, 3);

trial_cdata(~isfinite(trial_cdata)) = 0.0;
baseline_img = squeeze(mean(trial_cdata(:,:,base_begin:base_end),3));
baseline = repmat(baseline_img, [1 1 nframes]);

df_f = (trial_cdata - baseline) ./ baseline;
df_f(~isfinite(df_f)) = 0.0;

[x, y] = meshgrid(1:y_size, 1:x_size);
inpoly = inpolygon( x, y, roi(:,1), roi(:,2) );
df_f_in_roi = squeeze(sum(sum(double(df_f).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));

end

