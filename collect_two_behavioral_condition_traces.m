function [ btraces_per_condition, ctraces_in_roi_per_condition ] = collect_two_behavioral_condition_traces( condition_trials, cdata_raw, bdata_vel, rois )

ac = get_analysis_constants;
settings = sensor_settings;

PLANES = size(cdata_raw{1},4);

prestim = settings.pre_stim;
stim    = settings.stim;
poststim    = settings.post_stim;

base_begin = 1;
base_end = floor(prestim*VPS);

total_time = prestim + stim + poststim;

first_stim_t = prestim;
last_stim_t = stim + prestim;

[ysize, xsize] = size(squeeze(cdata_raw{1}(1,:,:,1,1,1)));
[x, y] = meshgrid(1:xsize, 1:ysize);

nframes = size( squeeze(cdata_raw{1}(1,:,:,1,1,:)), 3 );

for trial_type = 1:size( condition_trials, 1 )
        
    for cond_ord = 1:2
        
        cur_condition_trial_ords = condition_trials{ cond_ord, trial_type };
        
        for trial_ord_idx = 1:size( cur_condition_trial_ords, 1 )
        
            cur_trial_ord = cur_condition_trial_ords( trial_ord_idx );
            
            trial_cdata = squeeze( cdata_raw{ trial_type }( cur_trial_ord, :,:,:,:,:) );            
        
            % Collect behavioral data
            btraces_per_condition( cond_ord, trial_type, trial_ord_idx, ac.VEL_FWD, : ) = squeeze( bdata_vel{ trial_type }( cur_trial_ord, ac.VEL_FWD, : ) );
            btraces_per_condition( cond_ord, trial_type, trial_ord_idx, ac.VEL_YAW, : ) = squeeze( bdata_vel{ trial_type }( cur_trial_ord, ac.VEL_YAW, : ) );
            
            % Collect calcium data
            for p=1:PLANES
                            
                cur_roi_idx = 1;
                cur_roi = rois{ p, cur_roi_idx };
                
                cdata_in_plane = squeeze(trial_cdata(:,:,p,:));
                
                while( ~isempty( cur_roi ) )
                    xv = cur_roi(:,1);
                    yv = cur_roi(:,2);
                    
                    inpoly = inpolygon(x,y,xv,yv);

                    tmp = squeeze(sum(sum(cdata_in_plane .* repmat(inpoly, [1, 1, nframes])))) / sum(inpoly(:));
                    baseline = repmat(mean(tmp(base_begin:base_end)), [1 1 size(tmp,2)]);
                    itrace = (tmp-baseline) ./ baseline;                
                    
                    ctraces_in_roi_per_condition( cond_ord, trial_type, trial_ord_idx, cur_roi_idx, : ) = itrace;
                    
                    cur_roi_idx = cur_roi_idx + 1;
                    if( cur_roi_idx > size(rois,2))
                        break;
                    end
                    
                    cur_roi = rois{p, cur_roi_idx};
                end                
            end
        end
        
        pause(0.01); % This is so I can cancel if the program runs for too long.
    end
end
end

