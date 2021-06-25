



%% targeting period
ego_back_path = '/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210207_detail_for_wf_an_2s_for_t/beta_back_2s/';
ego_left_path = '/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210207_detail_for_wf_an_2s_for_t/beta_left_2s/';
ego_right_path = '/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210207_detail_for_wf_an_2s_for_t/beta_right_2s/';

ego_back_dir = dir(fullfile(ego_back_path, 'mean*.nii'));
ego_left_dir = dir(fullfile(ego_left_path, 'mean*.nii'));
ego_right_dir = dir(fullfile(ego_right_path, 'mean*.nii'));

ego_back_pool =  strcat(ego_back_path, {ego_back_dir.name}');
ego_left_pool =  strcat(ego_left_path, {ego_left_dir.name}');
ego_right_pool =  strcat(ego_right_path, {ego_right_dir.name}');

for ith_sub = 8:26
    
    ego_back_sub_path = ego_back_pool{find(cellfun(@(x) contains(x, ['sub' num2str(ith_sub)]), ego_back_pool))};
    ego_left_sub_path = ego_left_pool{find(cellfun(@(x) contains(x, ['sub' num2str(ith_sub)]), ego_left_pool))};
    ego_right_sub_path = ego_right_pool{find(cellfun(@(x) contains(x, ['sub' num2str(ith_sub)]), ego_right_pool))};
    
    ima_struct_shell = MRIread(ego_back_sub_path);
    ego_back_struct = MRIread(ego_back_sub_path);
    ego_left_struct = MRIread(ego_left_sub_path);
    ego_right_struct = MRIread(ego_right_sub_path);
    
    ego_back_ima = ego_back_struct.vol;
    ego_left_ima = ego_left_struct.vol;
    ego_right_ima = ego_right_struct.vol;
    
    % back - both
    final_ima = ego_back_ima - (ego_left_ima + ego_right_ima)/2;
    ima_struct_shell.vol = final_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210207_detail_for_wf_an_2s_for_t/beta_back_both_2s/mean_mni_sub', num2str(ith_sub), '.nii' ]);
    
    % both - back
    final_ima = (ego_left_ima + ego_right_ima)/2 - ego_back_ima;
    ima_struct_shell.vol = final_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210207_detail_for_wf_an_2s_for_t/beta_both_back_2s/mean_mni_sub', num2str(ith_sub), '.nii']);
    
    % left - back
    final_ima = ego_left_ima - ego_back_ima;
    ima_struct_shell.vol = final_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210207_detail_for_wf_an_2s_for_t/beta_LB_2s/mean_mni_sub', num2str(ith_sub), '.nii']);
    
    % right - back
    final_ima = ego_right_ima - ego_back_ima;
    ima_struct_shell.vol = final_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210207_detail_for_wf_an_2s_for_t/beta_RB_2s/mean_mni_sub', num2str(ith_sub), '.nii']);
end







%% targeting period
ego_back_t_path = '/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_back_2s_t/';
ego_left_t_path = '/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_left_2s_t/';
ego_right_t_path = '/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_right_2s_t/';

ego_back_tn_path = '/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_back_2s_tn/';
ego_left_tn_path = '/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_left_2s_tn/';
ego_right_tn_path = '/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_right_2s_tn/';

ego_back_t_dir = dir(fullfile(ego_back_t_path, 'mean*.nii'));
ego_left_t_dir = dir(fullfile(ego_left_t_path, 'mean*.nii'));
ego_right_t_dir = dir(fullfile(ego_right_t_path, 'mean*.nii'));

ego_back_tn_dir = dir(fullfile(ego_back_tn_path, 'mean*.nii'));
ego_left_tn_dir = dir(fullfile(ego_left_tn_path, 'mean*.nii'));
ego_right_tn_dir = dir(fullfile(ego_right_tn_path, 'mean*.nii'));

ego_back_t_pool =  strcat(ego_back_t_path, {ego_back_t_dir.name}');
ego_left_t_pool =  strcat(ego_left_t_path, {ego_left_t_dir.name}');
ego_right_t_pool =  strcat(ego_right_t_path, {ego_right_t_dir.name}');

ego_back_tn_pool =  strcat(ego_back_tn_path, {ego_back_tn_dir.name}');
ego_left_tn_pool =  strcat(ego_left_tn_path, {ego_left_tn_dir.name}');
ego_right_tn_pool =  strcat(ego_right_tn_path, {ego_right_tn_dir.name}');

for ith_sub = 8:26
    
    ego_back_t_sub_path = ego_back_t_pool{find(cellfun(@(x) contains(x, ['sub' num2str(ith_sub)]), ego_back_t_pool))};
    ego_left_t_sub_path = ego_left_t_pool{find(cellfun(@(x) contains(x, ['sub' num2str(ith_sub)]), ego_left_t_pool))};
    ego_right_t_sub_path = ego_right_t_pool{find(cellfun(@(x) contains(x, ['sub' num2str(ith_sub)]), ego_right_t_pool))};
    
    ego_back_tn_sub_path = ego_back_tn_pool{find(cellfun(@(x) contains(x, ['sub' num2str(ith_sub)]), ego_back_tn_pool))};
    ego_left_tn_sub_path = ego_left_tn_pool{find(cellfun(@(x) contains(x, ['sub' num2str(ith_sub)]), ego_left_tn_pool))};
    ego_right_tn_sub_path = ego_right_tn_pool{find(cellfun(@(x) contains(x, ['sub' num2str(ith_sub)]), ego_right_tn_pool))};
    
    ima_struct_shell = MRIread(ego_back_t_sub_path);
    
    ego_back_t_struct = MRIread(ego_back_t_sub_path);
    ego_left_t_struct = MRIread(ego_left_t_sub_path);
    ego_right_t_struct = MRIread(ego_right_t_sub_path);
    
    ego_back_tn_struct = MRIread(ego_back_tn_sub_path);
    ego_left_tn_struct = MRIread(ego_left_tn_sub_path);
    ego_right_tn_struct = MRIread(ego_right_tn_sub_path);
    
    ego_back_t_ima = ego_back_t_struct.vol;
    ego_left_t_ima = ego_left_t_struct.vol;
    ego_right_t_ima = ego_right_t_struct.vol;
    
    ego_back_tn_ima = ego_back_tn_struct.vol;
    ego_left_tn_ima = ego_left_tn_struct.vol;
    ego_right_tn_ima = ego_right_tn_struct.vol;
    
    % back - both
    final_t_ima = ego_back_t_ima - (ego_left_t_ima + ego_right_t_ima)/2;
    ima_struct_shell.vol = final_t_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_back_both_2s_t/mean_mni_sub', num2str(ith_sub), '.nii' ]);
    
    % both - back
    final_t_ima = (ego_left_t_ima + ego_right_t_ima)/2 - ego_back_t_ima;
    ima_struct_shell.vol = final_t_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_both_back_2s_t/mean_mni_sub', num2str(ith_sub), '.nii']);
    
    % left - right
    final_t_ima = ego_left_t_ima - ego_right_t_ima;
    ima_struct_shell.vol = final_t_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_LR_2s_t/mean_mni_sub', num2str(ith_sub), '.nii']);
    
    % right - left
    final_t_ima = ego_right_t_ima - ego_left_t_ima;
    ima_struct_shell.vol = final_t_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_RL_2s_t/mean_mni_sub', num2str(ith_sub), '.nii']);
    
    
    % back - both tn
    final_tn_ima = ego_back_tn_ima - (ego_left_tn_ima + ego_right_tn_ima)/2;
    ima_struct_shell.vol = final_tn_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_back_both_2s_tn/mean_mni_sub', num2str(ith_sub), '.nii' ]);
    
    % both - back tn
    final_tn_ima = (ego_left_tn_ima + ego_right_tn_ima)/2 - ego_back_tn_ima;
    ima_struct_shell.vol = final_tn_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_both_back_2s_tn/mean_mni_sub', num2str(ith_sub), '.nii']);
    
    % left - right tn
    final_tn_ima = ego_left_tn_ima - ego_right_tn_ima;
    ima_struct_shell.vol = final_tn_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_LR_2s_tn/mean_mni_sub', num2str(ith_sub), '.nii']);
    
    % right - left tn
    final_tn_ima = ego_right_tn_ima - ego_left_tn_ima;
    ima_struct_shell.vol = final_tn_ima;
    MRIwrite(ima_struct_shell,  ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210206/beta_RL_2s_tn/mean_mni_sub', num2str(ith_sub), '.nii']);
end

