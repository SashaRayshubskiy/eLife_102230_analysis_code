function plot_all_behaviour_with_gcamp_in_roi_in_2D( btrial_meta, cdata_raw, bdata_vel, bdata_vel_time, frame_start_offsets_per_plane, PLANE, VPS )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ac = get_analysis_constants;
settings = sensor_settings;
colorindex = 0;
order    = ac.order;

bdata_seq = convert_bmeta_to_sequential_data(btrial_meta);

INTER_TRIAL_T = 5.0;

nframes = size( cdata_raw{ 1 }, 5 );
x_size = size( cdata_raw{ 1 }, 2 );
y_size = size( cdata_raw{ 1 }, 3 );

t = (([0:nframes-1]))./VPS + frame_start_offsets_per_plane(PLANE);
[x, y] = meshgrid(1:y_size, 1:x_size);

% Get 2 ROIs 
f11 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
imagesc(squeeze(mean(squeeze(cdata_raw{ 1 }( 40, :, :, PLANE, : )),3)));
colormap gray;
axis image;

roiId = 1;
inpoly = [];
while(1)
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
        break
    end
    inpoly(roiId,:,:) = inpolygon(x,y,xv,yv);
    roiId = roiId + 1;
    
    %draw the bounding polygons and label them
    currcolor    = order(1+mod(colorindex,size(order,1)),:);
    hold on;
    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    colorindex = colorindex + 1;
end

f22 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
cur_start_time = 0.0;
for i=1:length(bdata_seq)
    
    cur_tt = bdata_seq(i,1);
    cur_ord = bdata_seq(i,2);
    
    cur_bvel = squeeze(bdata_vel{ cur_tt }( cur_ord, ac.VEL_YAW, : ));
    
    cur_ct = t + cur_start_time;
    cur_cd = squeeze(cdata_raw{ cur_tt }( cur_ord, :, :, PLANE, : ));
    
    cur_cd_trace1 = squeeze(sum(sum(double(cur_cd).*repmat(squeeze(inpoly(1,:,:)), [1, 1, nframes]))))/sum(squeeze(inpoly(1,:)));
    cur_cd_trace2 = squeeze(sum(sum(double(cur_cd).*repmat(squeeze(inpoly(2,:,:)), [1, 1, nframes]))))/sum(squeeze(inpoly(2,:)));
    
    n = floor(length(cur_bvel) / length(cur_cd_trace1));
    round_len = n*length(cur_cd_trace1);
    
    cur_bvel_reformat = reshape(cur_bvel(1:round_len), [n, length(cur_cd_trace1)]);
    %cur_bvel_std = std(cur_bvel_reformat,1);
    cur_bvel_std = std(cur_bvel);
    
    hold on;
    scatter(cur_bvel_std, squeeze(mean(cur_cd_trace1)));    
    
    pause(0.1);
    
end

saveidx = 1;

saveas(f11, [analysis_path '/df_f_vs_yaw_vel_std_roi_' num2str(saveidx) '.fig']);
saveas(f11, [analysis_path '/df_f_vs_yaw_vel_std_roi_' num2str(saveidx) '.png']);

saveas(f22, [analysis_path '/df_f_vs_yaw_vel_std_tc_' num2str(saveidx) '.fig']);
saveas(f22, [analysis_path '/df_f_vs_yaw_vel_std_tc_' num2str(saveidx) '.png']);

end

