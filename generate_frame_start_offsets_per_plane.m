function [ frame_start_offsets_per_plane ] = generate_frame_start_offsets_per_plane( planes, b_rawdata, b_time )

ac = get_analysis_constants;
one_trial_bdata = squeeze(b_rawdata{ ac.BOTH }(1,:,:));
frame_clock     = squeeze(one_trial_bdata(:,5));

%figure;
%hold on;
%plot(b_time, frame_clock);

FSTATE_HIGH = 1;
FSTATE_LOW = 2; 
cur_state = FSTATE_LOW;
FRAME_STATE_CHANGE_THRESHOLD = 2.0;

frame_begins_t = [];
frame_ends_t = [];

for t=1:length(frame_clock)-1

    cur_frame_signal = frame_clock(t);
    next_frame_signal = frame_clock(t+1);
    
    if( cur_state == FSTATE_LOW )
        % disp(['cfs: ' num2str(cur_frame_signal) ' nfs: ' num2str(next_frame_signal)]);
        if((cur_frame_signal - next_frame_signal) < -1.0*FRAME_STATE_CHANGE_THRESHOLD )
            frame_begins_t(end+1) = b_time(t+1);
            cur_state = FSTATE_HIGH;
        end
    elseif( cur_state == FSTATE_HIGH )        
        if((cur_frame_signal - next_frame_signal) > FRAME_STATE_CHANGE_THRESHOLD )
            frame_ends_t(end+1) = b_time(t);
            cur_state = FSTATE_LOW;
        end
    else
        disp(['ERROR: State: ' num2str(cur_state) ' is not recognized.']);
    end
end

%plot(frame_begins_t, 5.0, 'xr');
%plot(frame_ends_t, 5.0, 'xb');

FPLANES = planes;

VOLUMES = floor(length(frame_begins_t) / FPLANES);

frame_begins_in_vol = reshape( frame_begins_t(1:VOLUMES*FPLANES), [ FPLANES, VOLUMES ] );

frame_start_offsets_per_plane_per_vol = zeros(FPLANES, VOLUMES);

for v = 1:size(frame_begins_in_vol,2)
    
    first_plane_time = frame_begins_in_vol(1,v);
    
    for p = 1:FPLANES
        frame_start_offsets_per_plane_per_vol(p,v) = frame_begins_in_vol(p,v) - first_plane_time;
    end
end
    
frame_start_offsets_per_plane = zeros(1, FPLANES);
for p = 1:FPLANES
    frame_start_offsets_per_plane(p) = squeeze(mean(frame_start_offsets_per_plane_per_vol(p,:),2));
end

if 0
f = figure;
hold on;

for p = 1:size(frame_begins_in_vol,1)
    avg_offset(p) = squeeze(mean(frame_start_offsets_per_plane(p,:),2));
    err_in_offset(p) = squeeze(std(frame_start_offsets_per_plane(p,:)));       
end

errorbar( [1:size(frame_begins_in_vol,1)], avg_offset, err_in_offset );
xlabel('Planes');
ylabel('Frame start offset (s)');

saveas(f, [analysis_path '/frame_start_offsets_per_plane.fig']);
saveas(f, [analysis_path '/frame_start_offsets_per_plane.png']);
end

end

