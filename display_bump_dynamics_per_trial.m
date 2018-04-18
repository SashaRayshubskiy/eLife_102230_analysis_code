function display_bump_dynamics_per_trial( cdata_raw, bdata_vel_time, bdata_vel, VPS, analysis_path )

SHOW_AVG_ONLY = 1;

% DEBUG flag section
DEBUG_DISPLAY_BUMP_DYNAMICS_FOR_ALL_TRIALS = 0;
DEBUG_SHOW_PVA = 0;
DEBUG_SHOW_TC_PER_WEDGE = 0;

first_stim_t = 3.0;
last_stim_t = 4.0;

COMPARE_INITIAL_TURNS = 0;
COMPARE_COUNTER_TURNS = 1;

if( COMPARE_COUNTER_TURNS == 1 )
    % Generate counter turn trials
    ct_list = generate_counter_turn_trials( bdata_vel_time, bdata_vel, analysis_path );
elseif( COMPARE_INITIAL_TURNS == 1 )
    ct_list = generate_initial_turn_trials( bdata_vel_time, bdata_vel, analysis_path );
end

% Examine the overall bump image and plot the wedges
cur_trial_img = squeeze(cdata_raw{1}( :,:,:,:,:));

cur_trial_img_data_p = permute( cur_trial_img, [1 4 3 2 5]);

bump_2D_per_trial = mean(squeeze(mean(cur_trial_img_data_p, 4)),4);

dX = 8;

bump_2D_down_per_trial = squeeze(mean(reshape( bump_2D_per_trial, [size(bump_2D_per_trial, 1), size(bump_2D_per_trial, 2), dX, size(bump_2D_per_trial,3)/dX] ),3));

bump_2D_all = squeeze(mean(bump_2D_down_per_trial));

f = figure;
hold on;
imagesc( bump_2D_all );

% These parameters will be specific per fly
diameter_in = 4;
diameter_out = 22;

% Params for 6f - 03
xCenter = 16;
yCenter = 14;

% Params for 6f - 02
%xCenter = 16;
%yCenter = 14;

% Params for 6f - 01
% xCenter = 18;
% yCenter = 13;

% % Params for op6s - 01
% xCenter = 16;
% yCenter = 13;

% Params for op6s - 02
%xCenter = 18;
%yCenter = 13;

[rois, wedge_unit_vectors ] = get_EPG_wedge_rois(diameter_in, diameter_out, xCenter, yCenter);

display_wedge_rois( rois );

WEDGE_COUNT = 16;

% Get 'dark' space ROI for df/f image
dark_roi_param_file = [analysis_path '/dark_roi_params.mat'];
if( exist( dark_roi_param_file, 'file' ) == 2 )
    ppp = load(dark_roi_param_file);
    dark_space_roi = ppp.dark_space_roi;
        
    xv = dark_space_roi(:,1);
    yv = dark_space_roi(:,2);
    plot(xv, yv, 'Linewidth', 1,'Color','b');
    text(mean(xv),mean(yv),num2str(1),'Color','b','FontSize',12);    
else    
    dark_space_roi = get_roi();
    save( dark_roi_param_file, 'dark_space_roi' );    
end

saveas(f, [analysis_path '/EPG_bump_rois.fig']);
saveas(f, [analysis_path '/EPG_bump_rois.png']);

nframes = size(cdata_raw{1}, 5 );

t = (([0:nframes-1]))./VPS;

[x_size, y_size] = size( bump_2D_all );
[x, y] = meshgrid(1:y_size, 1:x_size);

TRIAL_TYPE_CNT = 2;

PVA_angle_ts = cell(1,TRIAL_TYPE_CNT);
PVA_mag_ts   = cell(1,TRIAL_TYPE_CNT);

f1 = figure;
f2 = figure;

pause(1);

if( DEBUG_SHOW_PVA == 1)
    TEST_TRIAL_CNT = 100;
    START_TEST_TRIAL = 1;
    END_TEST_TRIAL = START_TEST_TRIAL + TEST_TRIAL_CNT;
    
%    f3 = figure('units','normalized','outerposition',[0 0 1 1]);
    f3 = figure();
    vid = VideoWriter([analysis_path '/bump_PVA_tracking_debug.avi']);
    vid.FrameRate = 1;
    vid.Quality = 100;
    open(vid);
end

ac = get_analysis_constants;

for tt = 1:TRIAL_TYPE_CNT

    if( tt == ac.BOTH )
        cur_color_avg = rgb('DimGray');
        cur_color_single = rgb('DarkGray');
    elseif( tt == ac.LEFT )
        cur_color_avg = rgb('FireBrick');
        cur_color_single = rgb('LightSalmon');
    elseif( tt == ac.RIGHT )
        cur_color_avg = rgb('SeaGreen');
        cur_color_single = rgb('PaleGreen');
    end
    
    start_t = 1;
    
    trial_cnt = size(cdata_raw{tt}, 1);
    
    PVA_angle_ts{tt} = zeros( trial_cnt, 2, nframes );
    PVA_mag_ts{tt} = zeros( trial_cnt, 1, nframes );
    
    for tr = 1:trial_cnt       

        cur_trial_img_data = squeeze(cdata_raw{tt}( tr,:,:,:,:));
        
        % Crush out the shortest axis.
        cur_trial_img_data_p = squeeze(mean(permute( cur_trial_img_data, [3 2 1 4]),3));
        
        % Bin down in the longest axis
        bump_2D_down_per_trial = squeeze(mean(reshape( cur_trial_img_data_p, [size(cur_trial_img_data_p, 1), dX, size(cur_trial_img_data_p, 2)/dX, size(cur_trial_img_data_p,3)] ),2));        
        
        % Get 'dark' space time course for proper df/f baseline
        inpoly_dark = inpolygon(x,y,dark_space_roi(:,1),dark_space_roi(:,2));
        dark_tc = squeeze(sum(sum(double(bump_2D_down_per_trial).*repmat(inpoly_dark, [1, 1, nframes]))))/sum(inpoly_dark(:));
        dark_baseline = repmat( squeeze(mean(dark_tc)), [length(dark_tc) 1] );
        
        % If baseline is negative, use default
        dark_baseline(find(dark_baseline < 0 )) = 100;
        
        % Calculate time courses in wedge 
        tc_per_wedge = zeros(WEDGE_COUNT, nframes);
        for w = 1:WEDGE_COUNT
            xv = rois{ w }{ 1 };
            yv = rois{ w }{ 2 };
            inpoly = inpolygon(x,y,xv,yv);
            
            wedge_tc = squeeze(sum(sum(double(bump_2D_down_per_trial).*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
            
            wedge_tc_df_f = (wedge_tc - dark_baseline) ./ dark_baseline;
            wedge_tc_df_f(find(wedge_tc_df_f < 0)) = 0.1;
            tc_per_wedge(w,:) = wedge_tc_df_f;            
        end
                
        % Calculate PVA for each time step
        wedge_df_f_sum = sum(tc_per_wedge);
        wedge_id_tc = [];
        PVA_ts = [];
        for ts = 1:size(tc_per_wedge,2)
            cur_ts_df_f = squeeze(tc_per_wedge(:,ts));
            
            cur_wedge_df_f_sum = squeeze(wedge_df_f_sum(ts));
            
            PVA_ts_x_0 = squeeze(wedge_unit_vectors(:,1));
            PVA_ts_x_1 = PVA_ts_x_0 .* cur_ts_df_f;
            PVA_ts_x_sum = sum(PVA_ts_x_1,1);
            PVA_ts_x = PVA_ts_x_sum / cur_wedge_df_f_sum;

            PVA_ts_y_0 = squeeze(wedge_unit_vectors(:,2));
            PVA_ts_y_1 = PVA_ts_y_0 .* cur_ts_df_f;
            PVA_ts_y_sum = sum(PVA_ts_y_1,1);
            PVA_ts_y = PVA_ts_y_sum / cur_wedge_df_f_sum;
            
            PVA_ts(1) = PVA_ts_x;
            PVA_ts(2) = PVA_ts_y;                       
            
            % PVA_ts = sum(wedge_unit_vectors .* repmat(cur_ts_df_f, [1 2])) ./ wedge_df_f_sum(ts);
            wedge_id = get_wedge_for_angle( PVA_ts );            

            if 0
                figure(f1);
                imagesc(squeeze(bump_2D_down_per_trial(:,:,ts)));
                caxis([0 2000]);
                colormap redbluecmap;
                colorbar;
                hold on;
                
                display_wedge_rois( rois );
                
                cur_ang = PVA_ts;
                
                % Map magnitude from 0-50 to 0-5
                %cur_mag_scaled = cur_mag * 10.0;
                cur_mag_scaled = 10.0;
                
                x1 = xCenter;
                x2 = xCenter + cur_mag_scaled * cur_ang( 1 );
                y1 = yCenter;
                y2 = yCenter + cur_mag_scaled * cur_ang( 2 );
                
                plot( [x1 x2], [y1 y2], 'color', 'b', 'LineWidth', 2.0 );
                waitforbuttonpress;
            end
            %figure(f);
            %hold on;
            %plot([xCenter xCenter+10*PVA_ts(1)],[yCenter yCenter+10*PVA_ts(2)], 'k', 'LineWidth', 1.0);
            %waitforbuttonpress;
            
            wedge_id_tc(ts) = wedge_id;
            
            PVA_angle_ts{tt}(tr,:,ts) = PVA_ts;
            
            
            % Consider taking an average of several wedges around this
            % wedge
            PVA_mag_ts{tt}(tr,ts) = squeeze(tc_per_wedge(wedge_id,ts));
        end        
        
        if(DEBUG_SHOW_TC_PER_WEDGE == 1)
            f4 = figure;
            imagesc(tc_per_wedge);
            hold on;
            plot( wedge_id_tc, 'x', 'color', 'r');
            waitforbuttonpress;
            
        end            
        
        if((tt == 1) && ( DEBUG_SHOW_PVA == 1))
            figure(f3);
            subplot(1,1,1);
            
            for ts = 1:size(bump_2D_down_per_trial, 3)
                imagesc(squeeze(bump_2D_down_per_trial(:,:,ts)));
                %caxis([0 8]);
                caxis([0 2000]);
                colormap redbluecmap;
                colorbar;
                hold on;
                
                display_wedge_rois( rois );
                
                cur_mag = PVA_mag_ts{tt}(tr,ts);
                cur_ang = squeeze(PVA_angle_ts{tt}(tr,:,ts));
                
                % Map magnitude from 0-50 to 0-5
                %cur_mag_scaled = cur_mag * 10.0;
                cur_mag_scaled = 10.0;
                
                x1 = xCenter;
                x2 = xCenter + cur_mag_scaled * cur_ang( 1 );
                y1 = yCenter;
                y2 = yCenter + cur_mag_scaled * cur_ang( 2 );
                
                plot( [x1 x2], [y1 y2], 'color', 'b', 'LineWidth', 2.0 );
                title(['Tr_type: ' num2str(tt) ' trial: ' num2str(tr) ' ts: ' num2str(ts)]);
                
                writeVideo(vid, getframe(f3));
                pause(0.5);
                % waitforbuttonpress;
                hold off;                
            end
            
            if( tr == END_TEST_TRIAL )
                close(vid);
            end
        end
        
        % Display bump dynamics for all trials
        if( DEBUG_DISPLAY_BUMP_DYNAMICS_FOR_ALL_TRIALS == 1 )
            figure(f1);
            subplot(2,1,tt)
            hold on;
            end_t   = start_t + nframes;
            imagesc( [start_t:end_t], [], tc_per_wedge );
            caxis([0 20]);
            
            start_t = start_t + size(tc_per_wedge, 2);
            
            if(tt == 1)
                title('Left trials');
            else
                title('Right trials');
            end
        end
    end                  
    
    % Plot t vs. bump mag    for counter-turn vs. no-counter-turn {left-red, right-green}
    % Plot t vs. bump angle  for counter-turn vs. no-counter-turn
    
    PVA_mag_per_cond = cell(1,2);
    PVA_ang_per_cond = cell(1,2);
    
    prestim_t = find( (t < first_stim_t) & (t > (first_stim_t-0.15)));
    if(length(prestim_t) == 0)
        disp('ERROR: prestim_t is empty');
    end
    
    REF_REG = 180.0;
    
    for cond = 1:2
        
        cur_cond_trials = ct_list{ tt, cond };
        
        PVA_mag_per_cond{ cond } = [];
        PVA_ang_per_cond{ cond } = [];
        
        for tr = 1:length(cur_cond_trials)
            cur_cond_tr = cur_cond_trials(tr);
            
            cur_PVA_mag = squeeze(PVA_mag_ts{tt}(cur_cond_tr,:));
            cur_PVA_ang = squeeze(PVA_angle_ts{tt}(cur_cond_tr, :, : ));
            
            PVA_mag_per_cond{ cond }(tr,:) = cur_PVA_mag;
            
            theta_ts = rad2deg(atan2( cur_PVA_ang(2,:) , cur_PVA_ang(1,:) ));
            
            idx = find(theta_ts < 0);
            theta_ts(idx) = theta_ts(idx) + 360;
            
            % Get the offset to rotate the rest of the phases in the
            % trialx
            pre_stim_theta_ts = theta_ts(prestim_t(end));
            deg_offset = REF_REG - pre_stim_theta_ts;
            
            % Rotate the dt before stim onset bump vector to a common reference angle (zero).
            
            theta_ts_rotated = theta_ts + deg_offset;
            
            idx = find( theta_ts_rotated < 0 );
            theta_ts_rotated(idx) = theta_ts_rotated(idx) + 360;

            idx = find( theta_ts_rotated > 360 );
            theta_ts_rotated(idx) = theta_ts_rotated(idx) - 360;
            
            PVA_ang_per_cond{ cond }(tr,:) = theta_ts_rotated;
        end
    end
    
    figure(f2);
    for cond = 1:2
        if(cond == 1)
            subplot_1 = 1;
            subplot_2 = 3;
        else
            subplot_1 = 2;
            subplot_2 = 4;
        end
        
        trial_cnt = size( PVA_ang_per_cond{ cond }, 1 );
        if( SHOW_AVG_ONLY == 0 )
            for tr = 1:trial_cnt
                subplot(2,2,subplot_1);
                hold on;
                cur_PVA_mag = squeeze( PVA_mag_per_cond{ cond }(tr,:) );
                plot(t, cur_PVA_mag, 'color', cur_color_single, 'LineWidth', 1.0 );
                
                
                subplot(2,2,subplot_2);
                hold on;
                cur_PVA_ang = squeeze( PVA_ang_per_cond{ cond }(tr,:) );
                plot(t, cur_PVA_ang, 'color', cur_color_single, 'LineWidth', 1.0 );
            end
        end
        
        % Plot averages for PVA Mag
        subplot(2,2,subplot_1);
        hold on;
        avg_PVA_mag = mean( PVA_mag_per_cond{ cond } );
        sem_PVA_mag = get_sem( PVA_mag_per_cond{ cond }, 1);
        
        fh = fill( [t, fliplr(t)], ...
            [(avg_PVA_mag+sem_PVA_mag) fliplr((avg_PVA_mag-sem_PVA_mag))], ...
            cur_color_single);
        set(fh, 'EdgeColor', 'None');
        
        plot(t, avg_PVA_mag, 'color', cur_color_avg, 'LineWidth', 2.0 );
        
        xlim([0, t(end)]);
        ylim([0 50]);
        
        yy = ylim;
        y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');
        
        if( COMPARE_COUNTER_TURNS == 1)
            if(cond == 1)
                title('Counter-turn present, bump dF/F mag');
            else
                title('Counter-turn absent, bump dF/F mag');
            end
        elseif( COMPARE_INITIAL_TURNS == 1)
            if(cond == 1)
                title('Initial-turn present, bump dF/F mag');
            else
                title('Initial-turn absent, bump dF/F mag');
            end
        end
        
        % Plot averages for PVA Mag
        subplot(2,2,subplot_2);
        hold on;
        avg_PVA_ang = mean( PVA_ang_per_cond{ cond } );
        sem_PVA_ang = get_sem( PVA_ang_per_cond{ cond }, 1);
        
        fh = fill( [t, fliplr(t)], ...
            [(avg_PVA_ang + sem_PVA_ang) fliplr((avg_PVA_ang - sem_PVA_ang))], ...
            cur_color_single);
        set(fh, 'EdgeColor', 'None');
        
        plot(t, avg_PVA_ang, 'color', cur_color_avg, 'LineWidth', 2.0 );
        
        xlim([0, t(end)]);
        ylim([140 210]);
        
        yy = ylim;
        y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
        hh = fill([ first_stim_t first_stim_t last_stim_t last_stim_t ],[y_min y_max y_max y_min ], rgb('Wheat'));
        set(gca,'children',circshift(get(gca,'children'),-1));
        set(hh, 'EdgeColor', 'None');
        
        title('bump dF/F phase');
    end
end

%saveas(f1, [analysis_path '/bump_dynamics_all_data.fig']);
%saveas(f1, [analysis_path '/bump_dynamics_all_data.png']);

analysis_str = '';
if( COMPARE_COUNTER_TURNS == 1)
    analysis_str = 'counter_turn_vs_no_counter_turn';
elseif( COMPARE_INITIAL_TURNS == 1)
    analysis_str = 'initial_turn_vs_no_initial_turn';
end


if(SHOW_AVG_ONLY == 1)
    saveas(f2, [analysis_path '/bump_dynamics_' analysis_str '_avg.fig']);
    saveas(f2, [analysis_path '/bump_dynamics_' analysis_str '_avg.png']);
else
    saveas(f2, [analysis_path '/bump_dynamics_' analysis_str '.fig']);
    saveas(f2, [analysis_path '/bump_dynamics_' analysis_str '.png']);
end
end
