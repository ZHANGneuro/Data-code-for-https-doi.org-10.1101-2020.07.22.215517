





%% transform to 4d

root_path = ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_drection_20210207_detail_for_wf_an_2s_for_t/'];
each_cluster_dir = dir(fullfile(root_path, 'beta_RB_2s*'));
cluster_pool =  strcat(root_path, {each_cluster_dir.name}');

% root_path = ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_direction/'];
% each_cluster_dir = dir(fullfile(root_path, 'raw_ego*'));
% cluster_pool =  strcat(root_path, {each_cluster_dir.name}');

for ith_cluster = 1:length(cluster_pool)
    
    target_folder_path = cluster_pool{ith_cluster};
    
    brain_dir = dir(fullfile([target_folder_path, '/', 'mean_mni*']));
%     brain_dir = dir(fullfile([target_folder_path, '/', 'fsl152*']));
    brain_dir = strcat(target_folder_path, '/',  {brain_dir.name}');
    
    key_info = extractBetween(brain_dir{1}, 'analysis_contrast_drection_20210207_detail_for_wf_an_2s_for_t/', '/mean_mni');
%     key_info = extractBetween(brain_dir{1}, 'raw_ego_', '/fsl152');
    key_info = key_info{1};
    
    %% exclude voxel outside of brain
    mni152_strct = MRIread('/Users/bo/Documents/brainmask/fsl_standard/MNI152_T1_2mm_brain.nii.gz');
    mni152_ima=mni152_strct.vol;
    ith_pool = find(mni152_ima~=0);
    for n = 1:length(brain_dir)
        empty_brain = mni152_ima;
        empty_brain(:,:,:)=0;
        sub_num = cell2mat(extractBetween(brain_dir{n}, 'mean_mni_sub', '.nii'));
        sub_brain = MRIread(brain_dir{n});
        sub_ima = sub_brain.vol;
        empty_brain(ith_pool) = sub_ima(ith_pool);
        sub_brain.vol = empty_brain;
        MRIwrite(sub_brain, [target_folder_path, '/fsl152_sub', sub_num,'.nii']);
    end
    
    %% transform to 4d
    brain_dir = dir(fullfile(target_folder_path, 'fsl152*'));
    brain_dir = strcat(target_folder_path,  '/',{brain_dir.name}');
    spm('defaults','fmri');
    spm_jobman('initcfg');
    matlabbatch =[];
    matlabbatch{1}.spm.util.cat.vols = brain_dir;
    matlabbatch{1}.spm.util.cat.name = '4d.nii';
    matlabbatch{1}.spm.util.cat.dtype = 4;
    spm_jobman('run',matlabbatch);
    
    
    out_path = '/Users/bo/Desktop/';
    path_string = [target_folder_path, '/4d.nii'];
    imgs_strct = MRIread(path_string);
    [~,voxel_uncorrp] = ttest(imgs_strct.vol, 0, 'Dim', 4,'tail', 'right');
    voxel_uncorrp = 1- voxel_uncorrp;
    voxel_uncorrp(find(isnan(voxel_uncorrp)))=0;
    sub_brain.vol = voxel_uncorrp;
    MRIwrite(sub_brain, [out_path, '/unc_',key_info  , '.nii']);
    
end




