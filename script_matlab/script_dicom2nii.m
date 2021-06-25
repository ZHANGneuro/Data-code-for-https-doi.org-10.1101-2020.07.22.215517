



%%

root_path = ['/Users/bo/Documents/data_liujia_lab/task_greeble/pilot_beh_mri_results/MRI_raw_bold_4s_all_session_spm/'];
each_cluster_dir = dir(fullfile(root_path, 'LJZ*'));
cluster_pool =  strcat(root_path, {each_cluster_dir.name}');
sub_pool = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10'};
for ith_sub = 1:10
    cur_sub = sub_pool{ith_sub};
    for ith_sess = 1:8
        cur_file_list = cluster_pool(find(contains(cluster_pool, ['LJZB', cur_sub]) & contains(cluster_pool, ['sms_bold_run', num2str(ith_sess)])  ));
        for ith_file = 1:length(cur_file_list)
            cur_file_path = cur_file_list{ith_file};
            get_tr = extractBetween(cur_file_path, ['run', num2str(ith_sess), '-'], '.dcm');
            get_tr = str2num(get_tr{1});
            output_path = [root_path, 'sub', num2str(ith_sub), '_sess', num2str(ith_sess), '_TR', num2str(get_tr)];
            dicm2nii(cur_file_path, output_path, 1);
        end
    end
end


root_path = ['/Users/bo/Documents/data_liujia_lab/task_greeble/pilot_beh_mri_results/MRI_raw_bold_4s_all_session_spm/'];
each_cluster_dir = dir(fullfile(root_path, 'sub*'));
cluster_pool =  strcat(root_path, {each_cluster_dir.name}');
for ith_sub = 1:10
    for ith_sess = 1:8
        cur_file_list = cluster_pool(find(contains(cluster_pool, ['sub', num2str(ith_sub)]) & contains(cluster_pool, ['sess', num2str(ith_sess)])  ));
        for ith_file = 1:length(cur_file_list)
            cur_file_path = cur_file_list{ith_file};
            get_tr = strsplit(cur_file_path, 'TR');
            ith_tr = get_tr{2};
            source_file = [cur_file_path, '/sms_bold_run', num2str(ith_sess),'.nii.gz'];
            des_file = [root_path, '/bold_sub', num2str(ith_sub), '_s', num2str(ith_sess), '_tr', ith_tr, '.nii.gz'];
            copyfile(source_file, des_file);
        end
    end
end


root_path = ['/Users/bo/Documents/data_liujia_lab/task_greeble/pilot_beh_mri_results/MRI_raw_bold_4s_all_session_spm/'];
each_cluster_dir = dir(fullfile(root_path, 'bold*'));
cluster_pool =  strcat(root_path, {each_cluster_dir.name}');
for ith_sub = 1:10
        cur_file_list = cluster_pool(find(contains(cluster_pool, ['sub', num2str(ith_sub), '_'])));
        export_table = cell(length(cur_file_list), 3);
        sess_list = zeros(length(cur_file_list),1);
        tr_list = zeros(length(cur_file_list),1);
        for ith_row = 1:length(cur_file_list)
            cur_file_path = cur_file_list{ith_row};
            get_sess = extractBetween(cur_file_path, ['sub', num2str(ith_sub),'_s'], '_tr');
            get_tr = extractBetween(cur_file_path, 'tr', '.nii');
            export_table{ith_row, 1} = cur_file_path;
            export_table{ith_row, 2} = get_sess{1};
            export_table{ith_row, 3} = get_tr{1};
            sess_list(ith_row, 1) = str2num(get_sess{1});
            tr_list(ith_row, 1) = str2num(get_tr{1});
        end
        [B, I] = sortrows([sess_list, tr_list], 'ascend');
        export_table = export_table(I, :);

        %% transform to 4d
        spm('defaults','fmri');
        spm_jobman('initcfg');
        matlabbatch =[];
        matlabbatch{1}.spm.util.cat.vols = export_table(:,1);
        matlabbatch{1}.spm.util.cat.name = ['bold_4d_all_session_sub' , num2str(ith_sub), '.nii'];
        matlabbatch{1}.spm.util.cat.dtype = 4;
        spm_jobman('run',matlabbatch);
end


root_path = ['/Users/bo/Documents/data_liujia_lab/task_greeble/pilot_beh_mri_results/MRI_raw_bold_4s_all_session_spm/'];
each_cluster_dir = dir(fullfile(root_path, 'bold*'));
cluster_pool =  strcat(root_path, {each_cluster_dir.name}');
for ith_file = 1:length(cluster_pool)
    
    cur_file_path = cluster_pool{ith_file};
    get_file_first_part = strsplit(cur_file_path, '_tr');
    file_first_part = get_file_first_part{1};
    
    get_tr = strsplit(get_file_first_part{2}, '.');
    value_tr = get_tr{1};
    
    if length(value_tr)==1
        export_tr = ['00', num2str(value_tr)];
        des_path = [file_first_part, '_tr', export_tr, '.nii'];
        movefile(cur_file_path, des_path);
    end
    if length(value_tr)==2
        export_tr = ['0', num2str(value_tr)];
        des_path = [file_first_part, '_tr', export_tr, '.nii'];
        movefile(cur_file_path, des_path);
    end
end

%%
% root_path = ['/Users/bo/Documents/data_liujia_lab/task_greeble/pilot_beh_mri_results/'];
% each_cluster_dir = dir(fullfile(root_path, 'LJZ*'));
% cluster_pool =  strcat(root_path, {each_cluster_dir.name}');
% for ith_cluster = 1:length(cluster_pool)
%     curr_dir = cluster_pool{ith_cluster};
%     sub_dir = dir(fullfile(curr_dir, '*sms_bold_run*'));
%     sub_pool =  strcat(sub_dir(1).folder, '/',{sub_dir.name}');
%     
%     export_cell = cell(8, 2);
%     for ith_sess = 1:length(sub_pool)
%         cur_path = sub_pool{ith_sess};
%         curr_ima_strct = MRIread([cur_path, '/sms_bold_run', num2str(ith_sess),'.nii.gz']);
%         export_cell{ith_sess, 1} = curr_ima_strct;
%         export_cell{ith_sess, 2} = [cur_path, '/sms_bold_run', num2str(ith_sess),'.nii.gz'];
%     end
%     
%     export_struct = export_cell{1};
%     export_ima = export_struct.vol;
%     for ith_ima = 2:length(export_cell)
%         export_ima = cat(4, export_ima, export_cell{ith_ima}.vol);
%     end
%     export_struct.vol = export_ima;
%     MRIwrite(export_struct, [curr_dir, '/combined_4d_bold_sub', num2str(ith_cluster), '.nii.gz']);
% end




% 
% root_path = ['/Users/bo/Documents/data_liujia_lab/task_greeble/pilot_beh_mri_results/'];
% each_cluster_dir = dir(fullfile(root_path, 'LJZ*'));
% cluster_pool =  strcat(root_path, {each_cluster_dir.name}');
% for ith_cluster = 1:length(cluster_pool)
%     curr_dir = cluster_pool{ith_cluster};
%     sub_dir = dir(fullfile(curr_dir, '*gre_field_mapping_sms_mag-4.92*'));
%     sub_pool =  strcat(sub_dir(1).folder, '/',{sub_dir.name}');
%     for ith_sess = 1:length(sub_pool)
%         cur_path = sub_pool{ith_sess};
%         dicm2nii(cur_path, cur_path, 1)
%     end
% end


% root_path = ['/Users/bo/Documents/data_liujia_lab/task_greeble/pilot_beh_mri_results/'];
% each_cluster_dir = dir(fullfile(root_path, 'LJZ*'));
% cluster_pool =  strcat(root_path, {each_cluster_dir.name}');
% for ith_cluster = 1:length(cluster_pool)
%     curr_dir = cluster_pool{ith_cluster};
%     sub_dir = dir(fullfile(curr_dir, '*gre_field_mapping_sms_pha-7.38*'));
%     sub_pool =  strcat(sub_dir(1).folder, '/',{sub_dir.name}');
%     cur_path = sub_pool{1};
%     dicm2nii(cur_path, cur_path, 1)
% end

