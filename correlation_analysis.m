%% Get response timeocourse in ROI
REFERENCE_PLANE = 14;

fname = [analysis_path '\sid_' num2str(sid) '_comp_analysis_ref_plane_' num2str(REFERENCE_PLANE) '_clicky_df_f_' num2str(file_writer_cnt)];

BASELINE_START = 0.1;
BASELINE_END = 3.0;
STIM_START = 3.0;
STIM_END = 3.5;
VPS = 6.44;
[roi_points, intens] = clicky_df_f_custom_baseline_3_trial_types(AVG_DATA, REFERENCE_PLANE, VPS, BASELINE_START, BASELINE_END, STIM_START, STIM_END, fname);

file_writer_cnt = file_writer_cnt + 1;

%% Get correlation image in volume
SPACING = 0.01;
PADDING = 0;
MARGIN = 0.05;

BEGIN_CORR = 0.0;
END_CORR = 6.0;

begin_corr_idx = 1; %floor(BEGIN_CORR*FR);
end_corr_idx = floor(END_CORR*VPS);

DATA = double(DATA);

for tt = 1:3
    f2 = figure;
    for i=1:PLANES
        subaxis(SUBAXIS_ROW,SUBAXIS_COL, i, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
        rho = corr(squeeze(intens(1,begin_corr_idx:end_corr_idx))', reshape(squeeze(AVG_DATA{tt}(:,:,i,begin_corr_idx:end_corr_idx)), [size(AVG_DATA{tt},1)*size(AVG_DATA{tt},2) size(AVG_DATA{tt}(:,:,:,begin_corr_idx:end_corr_idx),4) ])' );
        corr_img = reshape(rho', [size(AVG_DATA{tt},1),  size(AVG_DATA{tt},2)]);
        imagesc( corr_img );
        %caxis([-0.2 0.2]);
        axis image;
        axis off;
        colormap jet;
        %colorbar;
    end


saveas(f2, [analysis_path '\sid_' num2str(sid) '_corr_in_volume_' num2str(file_writer_cnt) '_' num2str(tt) '.fig']);
saveas(f2, [analysis_path '\sid_' num2str(sid) '_corr_in_volume_' num2str(file_writer_cnt) '_' num2str(tt) '.png']);
close(f2);
end