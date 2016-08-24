function [ btraces_per_condition, avg_df_f_per_condition_per_plane ] = collect_two_behavioral_condition_and_df_f_per_cond( condition_trials, cdata_raw, bdata_vel, VPS, trial_exclusion_list, btrial_meta)

ac = get_analysis_constants;
settings = sensor_settings;

PLANES = size( cdata_raw{1}, 4 );

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

base_begin = 1;
baseline_time = 2.0;
base_end = floor(baseline_time*VPS);

total_time = prestim + stim + poststim;

first_stim_t = prestim;
last_stim_t = stim + prestim;

btraces_per_condition = cell(2,size( condition_trials, 1 ));

x_size = size(cdata_raw{ 1 }, 2);
y_size = size(cdata_raw{ 1 }, 3);
nframes = size(cdata_raw{ 1 }, 5);
avg_df_f_per_condition_per_plane = zeros( 3, 2, PLANES, x_size, y_size, nframes );

for trial_type = 1:size( condition_trials, 1 )
        
    for cond_ord = 1:2
        
        cur_condition_trial_ords = condition_trials{ trial_type, cond_ord };
        
        avg_raw_per_plane = zeros( PLANES, x_size, y_size, nframes );

        used_trial_idx = 1;
        for trial_ord_idx = 1:length( cur_condition_trial_ords)
        
            cur_trial_ord = cur_condition_trial_ords( trial_ord_idx );
            
            is_excluded = is_excluded_trial( cur_trial_ord, trial_exclusion_list{trial_type}, btrial_meta{trial_type} );
            if( is_excluded == 1 )
                continue;
            end
            
            disp(['About to process trial_type: ' num2str(trial_type) ' trial: ' num2str(cur_trial_ord)]);
        
            % Collect behavioral data
            btraces_per_condition{ cond_ord, trial_type }( used_trial_idx, ac.VEL_FWD, : ) = squeeze( bdata_vel{ trial_type }( cur_trial_ord, ac.VEL_FWD, : ) );
            btraces_per_condition{ cond_ord, trial_type }( used_trial_idx, ac.VEL_YAW, : ) = squeeze( bdata_vel{ trial_type }( cur_trial_ord, ac.VEL_YAW, : ) );
            
            % Collect calcium data
            trial_cdata = squeeze( cdata_raw{ trial_type }( cur_trial_ord, :,:,:,:,:) );            
            
            for p = 1:PLANES                           
                cdata_in_plane = squeeze(trial_cdata(:,:,p,:));
                
                avg_raw_per_plane(p,:,:,:) =  squeeze(avg_raw_per_plane(p,:,:,:)) + cdata_in_plane;

                if 0
                    baseline_img = squeeze(mean(cdata_in_plane(:,:,base_begin:base_end),3));
                    baseline = repmat(baseline_img, [1 1 nframes]);
                    
                    cur_df_f = (cdata_in_plane - baseline) ./ baseline;
                    
                    % get a df/f (x,y,v)
                    avg_df_f_per_condition_per_plane(cond_ord, p, :,:,:) = squeeze(avg_df_f_per_condition_per_plane(cond_ord, p, :,:,:)) + cur_df_f;
                end
            end
            
            used_trial_idx = used_trial_idx + 1;    
        end
    
        for p = 1:PLANES
            avg_raw_per_plane(p,:,:,:) = squeeze(avg_raw_per_plane(p,:,:,:)) ./ length( cur_condition_trial_ords );
            
            baseline_img = squeeze(mean(avg_raw_per_plane(p,:,:,base_begin:base_end),4));
            baseline = repmat(baseline_img, [1 1 nframes]);
            
            avg_df_f_per_condition_per_plane(trial_type, cond_ord, p, :,:,:) = (squeeze(avg_raw_per_plane(p,:,:,:)) - baseline) ./ baseline;
        end

        pause(0.01); % This is so I can cancel if the program runs for too long.
    end
end
end

