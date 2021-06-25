

sub_list = 8:26;
final_matrix = cell(length(sub_list), 4);
counter = 1;
for ith_sub = sub_list
    file_path = ['/Users/bo/Documents/data_yuji_lab/data_MEG/eyedate_fmri/el_', num2str(ith_sub),'.edf'];
    if exist(file_path,'file')
        edf = Edf2Mat(file_path);
        Result_orig = edf.RawEdf.FEVENT;
        sizeR = size(Result_orig);
        Result = cell(sizeR(2),9);
        for i  = 1:sizeR(2)
            tempme = Result_orig(i).message;
            tempme(find(isspace(tempme))) = [];
            Result{i,1} = i;
            Result{i,2} = tempme;
            Result{i,3} = Result_orig(i).codestring;
            Result{i,4} = Result_orig(i).gstx;
            Result{i,5} = Result_orig(i).gsty;
            Result{i,6} = Result_orig(i).genx;
            Result{i,7} = Result_orig(i).geny;
            Result{i,8} = Result_orig(i).sttime;
            Result{i,9} = Result_orig(i).entime;
        end
        Result(find(cellfun(@(x) isempty(x), Result(:,2))),2) = {'nan'};
        
        samples_struct = edf.RawEdf.FSAMPLE;
        time = samples_struct.time';
        time_start = min(time);
        time_end = max(time);
        
        sample_x = samples_struct.gx';
        sample_y = samples_struct.gy';
        if sample_x(1,1)<-10000
            the_sample = [sample_x(:,2) sample_y(:,2) time];
        elseif sample_x(1,2)<-10000
            the_sample = [sample_x(:,1) sample_y(:,1) time];
        end
        
        index_start_clip = find(cellfun(@(x) x==time_start, Result(:,8)));
        Result = Result(index_start_clip(end):end, :);
        
        for ith_ss = 1:4
            index_start = find(cellfun(@(x) strcmp(x, ['session:', num2str(ith_ss),'trial:0screen:ITI']), Result(:,2)));
            index_end = find(cellfun(@(x) strcmp(x, ['session:', num2str(ith_ss),'trial:40screen:trial_end']), Result(:,2)));
            
            if index_end
                temp_table = Result(index_start:index_end, :);
            else
                temp_table = Result(index_start:length(Result(:,1)), :);
            end
            
            final_matrix{counter, 1} = ith_sub;
            final_matrix{counter, 2} = ith_ss;
            final_matrix{counter, 3} = temp_table;
            final_matrix{counter, 4} = the_sample;
            counter = counter + 1;
            [num2str(ith_sub), '_',  num2str(ith_ss)]
        end
    end
end
save('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/mri_eyeData_raw_for_position.mat', 'final_matrix', '-v7.3');



% 
%%
final_matrix = load('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/mri_eyeData_raw_for_position.mat');
final_matrix_copy = final_matrix.final_matrix;
global_path = '/Users/bo/Documents/data_yuji_lab/data_fmri/';

for ith_row = 1:length(final_matrix_copy(:,1))
    curr_result_matrix = final_matrix_copy{ith_row, 3};
    curr_gaze_matrix = final_matrix_copy{ith_row, 4};
    cur_sub = final_matrix_copy{ith_row, 1};
    cur_sess = final_matrix_copy{ith_row, 2};
    
    %
    path_txt = ['/Users/bo/Documents/data_yuji_lab/data_MEG/eyedate_fmri/sub_',num2str(cur_sub), '_formal_rawdata.txt'];
    rawFile = fopen(path_txt,'rt');
    table_rawdata = textscan(rawFile,'%s %s %s %s %s %s %s %s %s %s %s %s %s');
    fclose(rawFile);
    [table_rawdata, table_time] = extract_behavioral_table_fmri (global_path, cur_sub);
    table_rawdata_by_session = cell(1,17);
    table_rawdata_by_time = cell(1,13);
    temp_ith = find( cellfun(@(x) strcmp(x, num2str(cur_sess)), table_time{1,13}));
    for each_col = 1:length(table_rawdata(1,:))
        table_rawdata_by_session{1, each_col} = table_rawdata{1, each_col}(temp_ith);
    end
    for each_col = 1 : length(table_time(1,:))
        table_rawdata_by_time{1, each_col} = table_time{1, each_col}(temp_ith);
    end
    need_choice_list = table_rawdata_by_session{1, 10};
    
    %
    inserted_matrix = cell(40,3);
    for ith_trial = 1:40
        trial_start_index = find(cellfun(@(x) strcmp(x, ['session:', num2str(cur_sess),'trial:', num2str(ith_trial-1), 'screen:show_target']), curr_result_matrix(:,2)));
        trial_end_index = find(cellfun(@(x) strcmp(x, ['session:', num2str(cur_sess),'trial:', num2str(ith_trial-1), 'screen:show_response_cue']), curr_result_matrix(:,2))); % 4s
%         trial_end_index = find(cellfun(@(x) strcmp(x, ['session:', num2str(cur_sess),'trial:', num2str(ith_trial-1), 'screen:target_noise']), curr_result_matrix(:,2))); % 2s
        
        iti_start_index = find(cellfun(@(x) strcmp(x, ['session:', num2str(cur_sess),'trial:', num2str(ith_trial-1), 'screen:ITI']), curr_result_matrix(:,2)));
        iti_end_index = find(cellfun(@(x) strcmp(x, ['session:', num2str(cur_sess),'trial:', num2str(ith_trial-1), 'screen:walking']), curr_result_matrix(:,2)));

        if trial_start_index
            trial_matrix_t = curr_result_matrix(trial_start_index:trial_end_index, :);
            inserted_matrix{ith_trial, 1} = trial_matrix_t;
            inserted_matrix{ith_trial, 2} = need_choice_list{ith_trial};
            inserted_matrix{ith_trial, 3} = [trial_matrix_t{1,8} trial_matrix_t{end,8}];
            
            trial_matrix_iti = curr_result_matrix(iti_start_index:iti_end_index, :);
            inserted_matrix{ith_trial, 4} = [trial_matrix_iti{1,8} trial_matrix_iti{end,8}];
        end
        [num2str(ith_row), '_',  num2str(ith_trial)]
    end
    final_matrix_copy{ith_row, 5} = inserted_matrix;
end

list_back = cell(1,6);
list_left = cell(1,6);
list_right = cell(1,6);
ith_row_back = 1;
ith_row_left = 1;
ith_row_right = 1;
for ith_row = 1:length(final_matrix_copy)
    
    curr_xy_pool = final_matrix_copy{ith_row, 4};
    curr_matrix = final_matrix_copy{ith_row, 5};

    back_index = find(cellfun(@(x) strcmp(x, 'back'), curr_matrix(:,2)));
    left_index = find(cellfun(@(x) strcmp(x, 'left'), curr_matrix(:,2)));
    right_index = find(cellfun(@(x) strcmp(x, 'right'), curr_matrix(:,2)));
    
    for ith = 1:length(back_index)
        get_the_index = back_index(ith);
        tp_start = curr_matrix{get_the_index,3}(1);
        tp_end = curr_matrix{get_the_index,3}(2);
        ith_start = find(curr_xy_pool(:,3)==tp_start);
        ith_end = find(curr_xy_pool(:,3)==tp_end);
        ITI_tp_start = curr_matrix{get_the_index,4}(1);
        ITI_tp_end = curr_matrix{get_the_index,4}(2);
        ith_ITI_start = find(curr_xy_pool(:,3)==ITI_tp_start);
        ith_ITI_end = find(curr_xy_pool(:,3)==ITI_tp_end);
        t_xy = curr_xy_pool(ith_start:ith_end, 1:2);
        t_xy = double(t_xy);
        t_xy(:,1) = t_xy(:,1)-1024/2;
        t_xy(:,2) = t_xy(:,2)-768/2;
        t_xy(find(t_xy(:,1)>1024/2 | t_xy(:,1)<=-1024/2 | t_xy(:,2)<=-768/2 | t_xy(:,2)>768/2),:)=[]; % carefull to add or not here
        t_xy(:,1) = t_xy(:,1) * (419/1024);
        t_xy(:,2) = t_xy(:,2) * (315/768);
        t_xy = radtodeg(atan(t_xy/751));

        bl_xy = curr_xy_pool(ith_ITI_start:ith_ITI_end, 1:2);
        bl_xy = double(bl_xy);
        bl_xy(find(bl_xy(:,1)>1024 | bl_xy(:,1)<=0 | bl_xy(:,2)<=0 | bl_xy(:,2)>768),:)=[];
        bl_xy(:,1) = bl_xy(:,1)-1024/2;
        bl_xy(:,2) = bl_xy(:,2)-768/2;
        bl_xy = mean(bl_xy, 1);
%         t_xy(:,1) = t_xy(:,1) - bl_xy(1);
%         t_xy(:,2) = t_xy(:,2) - bl_xy(2);
        list_back{ith_row_back, 1} = final_matrix_copy{ith_row, 1};
        list_back{ith_row_back, 2} = final_matrix_copy{ith_row, 2};
        list_back{ith_row_back, 3} = ith;
        list_back{ith_row_back, 4} = t_xy;
        list_back{ith_row_back, 5} = mean(t_xy(:,1));
        list_back{ith_row_back, 6} = mean(t_xy(:,2));
        ith_row_back = ith_row_back + 1;
    end
    for ith = 1:length(left_index)
        get_the_index = left_index(ith);
        tp_start = curr_matrix{get_the_index,3}(1);
        tp_end = curr_matrix{get_the_index,3}(2);
        ith_start = find(curr_xy_pool(:,3)==tp_start);
        ith_end = find(curr_xy_pool(:,3)==tp_end);
        ITI_tp_start = curr_matrix{get_the_index,4}(1);
        ITI_tp_end = curr_matrix{get_the_index,4}(2);
        ith_ITI_start = find(curr_xy_pool(:,3)==ITI_tp_start);
        ith_ITI_end = find(curr_xy_pool(:,3)==ITI_tp_end);
        t_xy = curr_xy_pool(ith_start:ith_end, 1:2);
        t_xy = double(t_xy);
        t_xy(:,1) = t_xy(:,1)-1024/2;
        t_xy(:,2) = t_xy(:,2)-768/2;
        t_xy(find(t_xy(:,1)>1024/2 | t_xy(:,1)<=-1024/2 | t_xy(:,2)<=-768/2 | t_xy(:,2)>768/2),:)=[]; % carefull to add or not here
        t_xy(:,1) = t_xy(:,1) * (419/1024);
        t_xy(:,2) = t_xy(:,2) * (315/768);
        t_xy = radtodeg(atan(t_xy/751));
        
        bl_xy = curr_xy_pool(ith_ITI_start:ith_ITI_end, 1:2);
        bl_xy = double(bl_xy);
        bl_xy(find(bl_xy(:,1)>1024 | bl_xy(:,1)<=0 | bl_xy(:,2)<=0 | bl_xy(:,2)>768),:)=[];
        bl_xy(:,1) = bl_xy(:,1)-1024/2;
        bl_xy(:,2) = bl_xy(:,2)-768/2;
        bl_xy = mean(bl_xy, 1);
%         t_xy(:,1) = t_xy(:,1) - bl_xy(1);
%         t_xy(:,2) = t_xy(:,2) - bl_xy(2);
        list_left{ith_row_left, 1} = final_matrix_copy{ith_row, 1};
        list_left{ith_row_left, 2} = final_matrix_copy{ith_row, 2};
        list_left{ith_row_left, 3} = ith;
        list_left{ith_row_left, 4} = t_xy;
        list_left{ith_row_left, 5} = mean(t_xy(:,1));
        list_left{ith_row_left, 6} = mean(t_xy(:,2));
        ith_row_left = ith_row_left + 1;
    end
    for ith = 1:length(right_index)
        get_the_index = right_index(ith);
        tp_start = curr_matrix{get_the_index,3}(1);
        tp_end = curr_matrix{get_the_index,3}(2);
        ith_start = find(curr_xy_pool(:,3)==tp_start);
        ith_end = find(curr_xy_pool(:,3)==tp_end);
        ITI_tp_start = curr_matrix{get_the_index,4}(1);
        ITI_tp_end = curr_matrix{get_the_index,4}(2);
        ith_ITI_start = find(curr_xy_pool(:,3)==ITI_tp_start);
        ith_ITI_end = find(curr_xy_pool(:,3)==ITI_tp_end);
        t_xy = curr_xy_pool(ith_start:ith_end, 1:2);
        t_xy = double(t_xy);
        t_xy(:,1) = t_xy(:,1)-1024/2;
        t_xy(:,2) = t_xy(:,2)-768/2;
        t_xy(find(t_xy(:,1)>1024/2 | t_xy(:,1)<=-1024/2 | t_xy(:,2)<=-768/2 | t_xy(:,2)>768/2),:)=[]; % carefull to add or not here
        t_xy(:,1) = t_xy(:,1) * (419/1024);
        t_xy(:,2) = t_xy(:,2) * (315/768);
        t_xy = radtodeg(atan(t_xy/751));
        
        bl_xy = curr_xy_pool(ith_ITI_start:ith_ITI_end, 1:2);
        bl_xy = double(bl_xy);
        bl_xy(find(bl_xy(:,1)>1024 | bl_xy(:,1)<=0 | bl_xy(:,2)<=0 | bl_xy(:,2)>768),:)=[];
        bl_xy(:,1) = bl_xy(:,1)-1024/2;
        bl_xy(:,2) = bl_xy(:,2)-768/2;
        bl_xy = mean(bl_xy, 1);
%         t_xy(:,1) = t_xy(:,1) - bl_xy(1);
%         t_xy(:,2) = t_xy(:,2) - bl_xy(2);
        list_right{ith_row_right, 1} = final_matrix_copy{ith_row, 1};
        list_right{ith_row_right, 2} = final_matrix_copy{ith_row, 2};
        list_right{ith_row_right, 3} = ith;
        list_right{ith_row_right, 4} = t_xy;
        list_right{ith_row_right, 5} = mean(t_xy(:,1));
        list_right{ith_row_right, 6} = mean(t_xy(:,2));
        ith_row_right = ith_row_right + 1;
    end
end
list_back(find(cellfun(@(x) isnan(x), list_back(:,5))),:)=[];
list_left(find(cellfun(@(x) isnan(x), list_left(:,5))),:)=[];
list_right(find(cellfun(@(x) isnan(x), list_right(:,5))),:)=[];

list_back = cell2mat(list_back(:, [1,2,3,5,6]));
list_left = cell2mat(list_left(:, [1,2,3,5,6]));
list_right = cell2mat(list_right(:, [1,2,3,5,6]));
% % compute visual angle
% hight_limit_half = radtodeg(atan((315/2)/751));
% width_limit_half = radtodeg(atan((419/2)/751));



%% plot heatmap for position
% col 4 5 w h
data_pool = {list_back, list_left, list_right};
label_pool = {'back', 'left','right'};
for ith_matrix = 1:length(data_pool)
    temp_matrix = data_pool{ith_matrix};
    fid = fopen(['/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/gaze_', label_pool{ith_matrix},'_bar_anova_4s_mri.txt'],'w');
    for n = 1:length(temp_matrix(:,1))
        fprintf(fid, '%3.4f\t %3.4f\t %3.4f\t %3.4f\t %3.4f\n', temp_matrix(n,:));
    end
    fclose(fid);
end














