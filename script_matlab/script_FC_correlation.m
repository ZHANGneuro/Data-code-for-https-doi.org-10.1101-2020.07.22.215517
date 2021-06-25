

global_path = '/Users/boo/Documents/degree_PhD/data_fmri'; % global_path =  '/Users/boo/Desktop/fmri_script/';
cluster_path = '/Users/boo/Documents/degree_PhD/data_fmri/analysis_contrast_direction/ROIs/';
cluster_dir = dir(fullfile(cluster_path, '*_final'));
cluster_list = {cluster_dir.name}';

ima_4d_b = '/Users/boo/Documents/degree_PhD/data_fmri/analysis_contrast_direction/raw_timeseries_b/';
ima_4d_dir_b = dir(fullfile(ima_4d_b, '*.nii'));
ima_4d_pool_b =  strcat(ima_4d_b, {ima_4d_dir_b.name}');

ima_4d_l = '/Users/boo/Documents/degree_PhD/data_fmri/analysis_contrast_direction/raw_timeseries_l/';
ima_4d_dir_l = dir(fullfile(ima_4d_l, '*.nii'));
ima_4d_pool_l =  strcat(ima_4d_l, {ima_4d_dir_b.name}');

ima_4d_r = '/Users/boo/Documents/degree_PhD/data_fmri/analysis_contrast_direction/raw_timeseries_r/';
ima_4d_dir_r = dir(fullfile(ima_4d_r, '*.nii'));
ima_4d_pool_r =  strcat(ima_4d_r, {ima_4d_dir_b.name}');

ima_4d_pool = {ima_4d_pool_b, ima_4d_pool_l, ima_4d_pool_r};

target_folder_b = '/Users/boo/Documents/degree_PhD/data_fmri/analysis_contrast_direction/raw_fc_r_map_b/';
target_folder_l = '/Users/boo/Documents/degree_PhD/data_fmri/analysis_contrast_direction/raw_fc_r_map_l/';
target_folder_r = '/Users/boo/Documents/degree_PhD/data_fmri/analysis_contrast_direction/raw_fc_r_map_r/';
output_folder = {target_folder_b, target_folder_l, target_folder_r};

for ith_roi = 1:length(cluster_list)
    
    roi_name = cluster_list{ith_roi};
    
    roi_timeseries_path = [cluster_path, roi_name];
    roi_timeseries_dir = dir(fullfile(roi_timeseries_path, '*.txt'));
    roi_timeseries_pool =  strcat(roi_timeseries_path, '/',{roi_timeseries_dir.name}');
    
    roi_ts_b = roi_timeseries_pool(find(contains(roi_timeseries_pool, '/b_')));
    roi_ts_l = roi_timeseries_pool(find(contains(roi_timeseries_pool, '/l_')));
    roi_ts_r = roi_timeseries_pool(find(contains(roi_timeseries_pool, '/r_')));
    
    roi_ts = {roi_ts_b, roi_ts_l, roi_ts_r};
    for ith_info = 1:3
        
        curr_info = roi_ts{ith_info};
        curr_4d_ima = ima_4d_pool{ith_info};
        curr_output_path = output_folder{ith_info};
        
        for ith_roi_timeseries= 1:length(curr_info)
            
            path_timeseries =curr_info{ith_roi_timeseries};
            path_roi = fopen(path_timeseries,'r');
            roi_series = textscan(path_roi,'%s');
            roi_series = str2double(roi_series{1});
            fclose(path_roi);
            
            ima_4d_path = curr_4d_ima{ith_roi_timeseries};
            ima_4d_struct = MRIread(ima_4d_path);
            ima_4d = ima_4d_struct.vol;
            
            sub = extractBetween(path_timeseries,  '_sub', '_s');
            session = extractBetween(path_timeseries, ['sub', sub{1},'_s'], '.');
            sub = str2num(sub{1});
            session = str2num(session{1});
            
            mean_ima4d = mean(ima_4d,4);
            [d1, d2, d3] = ind2sub(size(mean_ima4d),find(mean_ima4d~=0));
            ith_nonzero = find(mean_ima4d~=0);
            
            temp_1d_ima = zeros(length(d1), 1);
            
            empty_ima = ima_4d(:,:,:,1);
            empty_ima(:) = 0;
            
            parfor (ith_voxel = 1:length(d1), 6)
                
                voxelseries_4d = ima_4d(d1(ith_voxel), d2(ith_voxel), d3(ith_voxel), :);
                voxelseries_4d = voxelseries_4d(:);
                
                RHO = corr([roi_series voxelseries_4d]);
                temp_1d_ima(ith_voxel, 1) = RHO(2);
                
                [num2str(ith_roi), '-' ,num2str(ith_roi_timeseries), '-', num2str(ith_voxel)]
            end
            
            empty_ima(ith_nonzero) = 0.5 * (log((1+temp_1d_ima)./(1-temp_1d_ima)));
            
            ima_4d_struct.vol=empty_ima;
            MRIwrite(ima_4d_struct,[curr_output_path, roi_name ,'_sub_', num2str(sub), '_s',  num2str(session), '.nii']);
            
        end
    end
end











