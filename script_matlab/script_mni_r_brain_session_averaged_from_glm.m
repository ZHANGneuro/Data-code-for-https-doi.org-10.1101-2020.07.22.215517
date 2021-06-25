

root_path =  '/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_direction_20200217/';
% info_pool = {'beta_back_2s_t', 'beta_left_2s_t', 'beta_right_2s_t', 'beta_back_2s_tn', 'beta_left_2s_tn', 'beta_right_2s_tn'};
info_pool = {'beta_back_4s', 'beta_left_4s', 'beta_right_4s'};

for ith_info = 1:length(info_pool)
    
    period_key = info_pool{ith_info};
    
    for sub = 8:26
        
        target_folder_path = [root_path, period_key,'/'];
        brain_dir = dir(fullfile(target_folder_path, ['mni_sub',num2str(sub),'*']));
        brain_dir = strcat(target_folder_path,  {brain_dir.name}');
        
        s1_struct = MRIread(brain_dir{1});
        s2_struct = MRIread(brain_dir{2});
        s3_struct = MRIread(brain_dir{3});
        s4_struct = MRIread(brain_dir{4});
        s1_ima = s1_struct.vol;
        s2_ima = s2_struct.vol;
        s3_ima = s3_struct.vol;
        s4_ima = s4_struct.vol;
        
        ave_r_brain = (s1_ima + s2_ima + s3_ima +s4_ima)/4;
        temp_struc  = s1_struct;
        temp_struc.vol = ave_r_brain;
        MRIwrite(temp_struc,[root_path, period_key, '/mean_mni_sub', num2str(sub),'.nii']);
    end
end

