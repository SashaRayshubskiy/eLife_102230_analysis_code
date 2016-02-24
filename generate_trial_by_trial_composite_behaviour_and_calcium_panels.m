function generate_trial_by_trial_composite_behaviour_and_calcium_panels( asid, sid, cdata_raw, bdata_vel, btrial_meta, bdata_vel_time, VPS, analysis_path, rois )

global slash;

% tbt = trial-by-trial
tbt_analysis_path = [analysis_path '/asid_' num2str(asid) '_trial_by_trial/'];

if(~exist(tbt_analysis_path, 'dir'))
    mkdir( tbt_analysis_path );
end

aconstants = get_analysis_constants;
trial_type_cnt = length(cdata_raw);
settings = sensor_settings;

for trial_type = 1:trial_type_cnt
    
    cur_trial_cnt = size( cdata_raw{ trial_type }, 1 );
    cur_trial_type_str = aconstants.task_str{trial_type};
    
    for trial = 1:cur_trial_cnt
        %cur_trial_cdata = squeeze(cdata_raw{ trial_type }(trial,:,:,:,:,:));
        cur_trial_bdata = squeeze(bdata_vel{ trial_type }(trial,:,:));
        cur_trial_id = squeeze(btrial_meta{ trial_type }(trial, 2));
        
        cur_tbt_filename_prefix = [ tbt_analysis_path '/roi_avg_volume_sid_' num2str(sid) '_' cur_trial_type_str '_tid_' num2str(cur_trial_id)];
        generate_volume_avg_with_rois( squeeze(cdata_raw{ trial_type }(trial,:,:,:,:,:)), rois, cur_tbt_filename_prefix );
        
        cur_tbt_filename_prefix = [ tbt_analysis_path '/time_courses_in_volume_sid_' num2str(sid) '_' cur_trial_type_str '_tid_' num2str(cur_trial_id)];        
        generate_volume_time_courses_with_rois( squeeze(cdata_raw{ trial_type }(trial,:,:,:,:,:)), cur_trial_bdata, bdata_vel_time, rois, VPS, cur_tbt_filename_prefix );
    end
end

end

