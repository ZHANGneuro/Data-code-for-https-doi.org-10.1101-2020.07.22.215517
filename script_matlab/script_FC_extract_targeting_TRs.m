% function split_bandpass_residual_into_period

global_path = '/Users/boo/Documents/degree_PhD/data_fmri/'; % global_path =  '/gpfs/share/home/1401110602/fmri_script/';
ima_4d_root = '/Users/boo/Documents/degree_PhD/data_fmri/residual_fc_5mm_derivative_bandpass/';
ima_4d_dir = dir(fullfile(ima_4d_root, '*.nii'));
ima_4d_pool =  strcat(ima_4d_root, {ima_4d_dir.name}');

for ith_residual= 1:length(ima_4d_pool)
    
    ima_4d_path = ima_4d_pool{ith_residual};
    ima_4d_struct = MRIread(ima_4d_path);
    ima_4d = ima_4d_struct.vol;
    
    temp = extractBetween(ima_4d_path, 'bandpass_', 'nii');
    temp = temp{1};
    sub = extractBetween(temp,  'sub', '_s');
    session = extractBetween(temp, '_s', '.');
    sub = str2num(sub{1});
    session = str2num(session{1});
    
    % load behavior data, load into matrix, make it seconds
    [table_rawdata, table_time] = extract_behavioral_table_fmri(global_path, sub);
    col_count_raw_time = length(table_time(1,:));
    raw_time_by_session = cell(1,col_count_raw_time);
    raw_data_by_session = cell(1,col_count_raw_time);
    for each_col = 1:col_count_raw_time
        indexes = find(cellfun(@(x,y)  ~strcmp(x, 'NA') & strcmp(y, num2str(session)), ...
            table_time{1,12}, table_time{1,13}));
        raw_time_by_session{1, each_col} = table_time{1, each_col}(indexes);
        raw_data_by_session{1, each_col} = table_rawdata{1, each_col}(indexes);
    end
    num_row = length(raw_time_by_session{1,1});
    num_col = length(raw_time_by_session(1,:));
    txt_matrix=nan(num_row, num_col);
    for ith_txtcol = 1: num_col
        for ith_txtrow = 1: num_row
            txt_matrix(ith_txtrow, ith_txtcol) = str2num(raw_time_by_session{1, ith_txtcol}{ith_txtrow});
        end
    end
    txt_matrix = txt_matrix(:,[2, 4:11]);
    txt_matrix = txt_matrix* 0.001-6;
    TR_ith_matrix = round(txt_matrix/2);
    TR_ith_matrix( find(TR_ith_matrix(:,1)<0), 1) = 0;
    TR_ith_matrix = TR_ith_matrix + 1;
    TR_ith_matrix = TR_ith_matrix + 2;
    TR_ith_matrix(:,6:11) = TR_ith_matrix(:,4:9);
    TR_ith_matrix(:,4) = TR_ith_matrix(:,3)+1;
    TR_ith_matrix(:,5) = TR_ith_matrix(:,4)+1;
    
    
    TR_series_back = [];
    TR_series_left = [];
    TR_series_right = [];
    
    indexes_back = find(cellfun(@(x,y)  strcmp(x, 'back'), raw_data_by_session{1,10}));
    indexes_left = find(cellfun(@(x,y)  strcmp(x, 'left'), raw_data_by_session{1,10}));
    indexes_right = find(cellfun(@(x,y)  strcmp(x, 'right'), raw_data_by_session{1,10}));
    
    TR_ith_matrix_b = TR_ith_matrix(indexes_back, [9 10]);
    TR_ith_matrix_l = TR_ith_matrix(indexes_left, [9 10]);
    TR_ith_matrix_r = TR_ith_matrix(indexes_right, [9 10]);
    
    TR_ith_matrix_b = TR_ith_matrix_b(:);
    TR_ith_matrix_l = TR_ith_matrix_l(:);
    TR_ith_matrix_r = TR_ith_matrix_r(:);
    
    empty_ima_b = ima_4d(:,:,:, TR_ith_matrix_b);
    empty_ima_l = ima_4d(:,:,:, TR_ith_matrix_l);
    empty_ima_r = ima_4d(:,:,:, TR_ith_matrix_r);
    
    ima_4d_struct.vol=empty_ima_b;
    MRIwrite(ima_4d_struct,['/Users/boo/Documents/degree_PhD/data_fmri/analysis_contrast_direction/raw_timeseries_b/', 'fc_bandpass_sub_', num2str(sub), '_s',  num2str(session), '.nii']);
    
    ima_4d_struct.vol=empty_ima_l;
    MRIwrite(ima_4d_struct,['/Users/boo/Documents/degree_PhD/data_fmri/analysis_contrast_direction/raw_timeseries_l/', 'fc_bandpass_sub_', num2str(sub), '_s',  num2str(session), '.nii']);
    
    ima_4d_struct.vol=empty_ima_r;
    MRIwrite(ima_4d_struct,['/Users/boo/Documents/degree_PhD/data_fmri/analysis_contrast_direction/raw_timeseries_r/', 'fc_bandpass_sub_', num2str(sub), '_s',  num2str(session), '.nii']);
    
    ith_residual
    
end


