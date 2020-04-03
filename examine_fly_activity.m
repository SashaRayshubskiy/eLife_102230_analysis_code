function examine_fly_activity( basedir, directories )

figure;

DIR_NUM = length( directories );

for d = 1:DIR_NUM
      
    cur_dir = directories{d}{1};
    cur_sid = directories{d}{2};
    
    datapath = [ basedir '/' cur_dir ];
    
    disp( [ 'Processing: ' datapath ] );
    
    [ cur_data ] = load_EB_data_behaviour_only( datapath, cur_sid );
    
    ball_time = cur_data.bdata_vel_time;
    ball_data = cur_data.bdata_vel;

    cur_d = ball_data{ 1 };
    
    total_speed = [];
    for tr = 1:size(cur_d,1)
        
        cur_fwd  = squeeze( cur_d( tr, 1, : ) );
        cur_side = squeeze( cur_d( tr, 2, : ) );
        cur_yaw  = squeeze( cur_d( tr, 3, : ) );

        cur_speed = abs(cur_fwd) + abs(cur_side) + abs(cur_yaw);
        
        total_speed = vertcat( total_speed, cur_speed );
        
        % disp(['Trial: ' num2str(tr)]);
    end
    
    subplot(DIR_NUM, 1, d)
    histogram( total_speed );
    xlabel('Speed ');
    disp(['Avg speed: ' num2str( mean(total_speed) ) ]);
    
end

end

