%% Calibrate pixels to distance

filepath = '/data/drive0/sasha/160405_calibration_with_0.01mm_ruler/';
filename = 'calibration_512x128p_zoom_25_00001.tif';
load_path = [filepath filename];

raw_data = squeeze(open_tif_fast_single_plane( load_path ));

figure;
imagesc(squeeze(mean(raw_data,3)));