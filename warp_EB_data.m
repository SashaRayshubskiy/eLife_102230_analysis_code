function warp_EB_data( stim_events, VPS, analysis_path )

settings = sensor_settings;

PEAK_EB_BUMP_FRAME = 20;
REFERENCE_EB_TC = 14;

ball_FPS  = settings.sensorPollFreq;
ephys_FPS = settings.sampRate;

EB_tw_data = stim_events.EB;
yaw_tw_data = stim_events.yaw;
ephys_tw_data = stim_events.ephys;

num_stims = length(EB_tw_data);

f = figure;
subplot(3,1,1);
for i = 1:num_stims
    hold on;
    plot( EB_tw_data{i} );
    %title(['Stim id: ' num2str(i) ]);
    %waitforbuttonpress;
end

title('Full stim events');
%plot( mean(EB_tw_data), 'LineWidth', 4.0 );

subplot(3,1,2);
EB_frames_to_warp = cell( 1, num_stims );
filtered_EB_data = cell( 1, num_stims );

for i = 1:num_stims
    hold on;
    filtered_EB_data{i} = medfilt1( EB_tw_data{i}, 10);
    
    plot( filtered_EB_data{i} );
    
    EB_frames_to_warp{i} = filtered_EB_data{i}(PEAK_EB_BUMP_FRAME:end);
end

%plot( mean(filtered_EB_data), 'LineWidth', 4.0 );

subplot(3,1,3);

reference_tc = EB_frames_to_warp{ REFERENCE_EB_TC };

EB_warped   = cell(1, num_stims);
yaw_warped   = cell(1, num_stims);
ephys_warped = cell(1, num_stims);

for i = 1:num_stims
    hold on;
    
    tic;
    [dummy, ix, iy ] = dtw( reference_tc, EB_frames_to_warp{i} );
    toc;
    disp(['Finished dtw for i: ' num2str(i)]);
    
    EB_warped{i}    = EB_frames_to_warp{i}( iy );
    
    % Convert EB frames to yaw frames by using time
    
    YAW_START_FRAME = floor((PEAK_EB_BUMP_FRAME/VPS) * ball_FPS);
    yaw_iy = floor((iy/VPS) * ball_FPS);
    cur_yaw_snip = yaw_tw_data{i}( YAW_START_FRAME:end );
    
    yaw_iy_final = [];
    for j = 1:length(yaw_iy)
        cur_yaw_iy = yaw_iy(j);
        if( cur_yaw_iy < length(cur_yaw_snip))
            yaw_iy_final(j) = cur_yaw_iy;
        end
    end
    
    yaw_warped{i}   = cur_yaw_snip(yaw_iy_final);
  
    EPHYS_START_FRAME = floor((PEAK_EB_BUMP_FRAME/VPS) * ephys_FPS);
    ephys_iy = floor((iy/VPS) * ephys_FPS);
    cur_ephys_snip = ephys_tw_data{i}( EPHYS_START_FRAME:end );
    
    ephys_iy_final = [];
    for j = 1:length(ephys_iy)
        cur_ephys_iy = ephys_iy(j);
        if( cur_ephys_iy < length(cur_ephys_snip))
            ephys_iy_final(j) = cur_yaw_iy;
        end
    end
        
    ephys_warped{i}   = cur_ephys_snip( ephys_iy_final );
    
    plot( EB_warped{i} );
end

title('Peak of the bump move to baseline');
xlabel('EB Frame');

saveas(f, [analysis_path '/time_warped_EB_traces_upsample.fig']);
saveas(f, [analysis_path '/time_warped_EB_traces_upsample.png']);

% Compute averages of yaw and ephys after the EB based time warp
f = figure;

subplot(3,1,1);
max_len = -1;
for i = 1:num_stims
    cur_len = length( EB_warped{i} );
    if( cur_len > max_len )
        max_len = cur_len;
    end
end

EB_warped_array = zeros( num_stims, max_len );

hold on;
for i = 1:num_stims
    EB_warped_array(i,1:length(EB_warped{i})) = EB_warped{i};    
    plot( EB_warped_array(i,:) );
end

plot( mean(EB_warped_array), 'LineWidth', 3.0 );

subplot(3,1,2);
max_len = -1;
for i = 1:num_stims
    cur_len = length( yaw_warped{i} );
    if( cur_len > max_len )
        max_len = cur_len;
    end
end

yaw_warped_array = zeros( num_stims, max_len );

hold on;
for i = 1:num_stims
    yaw_warped_array(i,1:length(yaw_warped{i})) = yaw_warped{i};    
    plot( yaw_warped_array(i,:) );
end

plot( mean(yaw_warped_array), 'LineWidth', 3.0 );

subplot(3,1,3);
max_len = -1;
for i = 1:num_stims
    cur_len = length( yaw_warped{i} );
    if( cur_len > max_len )
        max_len = cur_len;
    end
end

ephys_warped_array = zeros( num_stims, max_len );

hold on;
for i = 1:num_stims
    ephys_warped_array(i,1:length(ephys_warped{i})) = ephys_warped{i};    
    plot( ephys_warped_array(i,:) );
end

plot( mean( ephys_warped_array ), 'LineWidth', 3.0 );

saveas(f, [analysis_path '/avg_warped_yaw_ephys_after_EB_time_warp.fig']);
saveas(f, [analysis_path '/avg_warped_yaw_ephys_after_EB_time_warp.png']);

end