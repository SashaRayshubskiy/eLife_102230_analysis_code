function extract_key_variables( basedir, directories )

for d = 1:length( directories )
      
    cur_dir = directories{d}{1};
    cur_sid = directories{d}{2};
    
    datapath = [ basedir '/' cur_dir ];
    
    disp( [ 'Processing: ' datapath ] );
    
    [ cur_data ] = load_EB_data( datapath, cur_sid );
    
    analysis_path = [datapath  '/analysis'];
    
    if(~exist(analysis_path, 'dir'))
        mkdir(analysis_path);
    end
    
    % Get max projection 
    START_PLANE = 5; % For 16 planes
    EPG_data_1       = squeeze( cur_data.cdata_raw{1}( :, :, :, 1, START_PLANE:end, : ));
    % Get max for each plane
    % EPG_data = squeeze(max(EPG_data_1, [], 4));

    glom = get_EB_glomeruli_16_planes( EPG_data_1, analysis_path );

    % Use trials_to_include from this point forward
    glom_EPG_F_per_trial = get_PB_F_per_trial( glom, EPG_data_1 );
    
    ephys_time = cur_data.ephys_time;
    ephys_data = cur_data.ephys_data;
    
    ball_time = cur_data.bdata_vel_time;
    ball_data = cur_data.bdata_vel;
        
    stim_data = cur_data.pico_stim_data;
    VPS       = cur_data.VPS;
    
    cur_var_path = [analysis_path '/raw_glom_ball_ephys.mat'];
    save( cur_var_path, 'glom_EPG_F_per_trial', 'ephys_time', 'ephys_data', 'ball_time', 'ball_data', 'stim_data', 'VPS' );
    disp( [ 'Wrote: ' cur_var_path ] );
    
    clear cur_data;
    clear EPG_data_1; 
    clear EPG_data;
end

end

