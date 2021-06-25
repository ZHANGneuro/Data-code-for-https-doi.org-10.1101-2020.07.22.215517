

%% transform to 4d
key_wordd='folder_fc_back_parietal';
path = ['/Users/bo/Documents/data_yuji_lab/data_fmri/analysis_contrast_direction/' ,key_wordd,  '/'];

mni152_strct = MRIread('/Users/bo/Documents/data_yuji_lab/data_fmri/brainmask/fsl_standard/MNI152_T1_2mm_brain_mask.nii.gz');
mni152_ima=mni152_strct.vol;

% ith_pool = find(mni152_ima~=0);
% brain_dir = dir(fullfile( [ path , 'fsl*'] ));
% brain_dir = strcat(path,  {brain_dir.name}');
% for n = 1:length(brain_dir)
%     empty_brain = mni152_ima;
%     empty_brain(:,:,:)=0;
%     sub_num = cell2mat(extractBetween(brain_dir{n}, 'fsl152_sub', '.nii'));
%     sub_brain = MRIread(brain_dir{n});
%     sub_ima = sub_brain.vol;
%     empty_brain(ith_pool) = sub_ima (ith_pool);
%     sub_brain.vol = empty_brain;
%     MRIwrite(sub_brain, [path, 'fsl152_sub', num2str(sub_num), '.nii']);
% end

%% transform to 4d
brain_dir = dir(fullfile( [ path , 'mean*'] ));
brain_dir = strcat(path,  {brain_dir.name}');
spm('defaults','fmri');
spm_jobman('initcfg');
matlabbatch =[];
matlabbatch{1}.spm.util.cat.vols = brain_dir;
matlabbatch{1}.spm.util.cat.name = ['4d.nii'];
matlabbatch{1}.spm.util.cat.dtype = 4;
spm_jobman('run',matlabbatch);

%% export p and t value
imgs_strct = MRIread([path, '4d.nii']);
% out_path = path;
out_path= '/Users/boo/Desktop';
[~,voxel_uncorrp,~, STATS] = ttest(imgs_strct.vol, 0, 'Dim', 4,'tail', 'right');
voxel_uncorrp = 1- voxel_uncorrp;
voxel_uncorrp(find(isnan(voxel_uncorrp)))=0;

imgs_output = MRIread([path, 'mean_mni_sub8_parietal.nii']);
imgs_output.vol = voxel_uncorrp;
MRIwrite(imgs_output, [out_path, '/unc_' ,key_wordd, '.nii']);

% imgs_output.vol = STATS.tstat;
% MRIwrite(imgs_output, [out_path, '/tvalue_' ,key_wordd, '.nii']);






% %%
% path = '/Users/boo/Documents/degree_PhD/data_fmri/analysis_contrast_period/hnd_statistic_ima/';
% brain_dir = dir(fullfile( [path, 'tvalue_hnd_*'] ));
% brain_dir = strcat(path,  {brain_dir.name}');
% 
% for ith = 1:length(brain_dir)
%     ima_path = brain_dir{ith};
%     keyword = extractBetween(ima_path, 'tvalue_hnd_data_', '.nii');
%     keyword = keyword{1};
%     imgs_strct = MRIread(ima_path);
%     img = imgs_strct.vol;
%     img(find(img<2.55 & img>-2.55))=NaN;
%     imgs_strct.vol = img;
%     MRIwrite(imgs_strct, [path, 'tvalue_cropped_2.5_' ,keyword, '.nii']);
% end

