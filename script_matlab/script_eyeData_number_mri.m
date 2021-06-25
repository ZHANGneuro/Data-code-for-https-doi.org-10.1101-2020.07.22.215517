

final_matrix = cell(length(8:26), 4);
for ith_sub = 8:26
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
        index_s1_start = find(cellfun(@(x) strcmp(x, 'session:1trial:0screen:ITI'), Result(:,2)));
        index_s2_start = find(cellfun(@(x) strcmp(x, 'session:2trial:0screen:ITI'), Result(:,2)));
        index_s3_start = find(cellfun(@(x) strcmp(x, 'session:3trial:0screen:ITI'), Result(:,2)));
        index_s4_start = find(cellfun(@(x) strcmp(x, 'session:4trial:0screen:ITI'), Result(:,2)));
        
        final_matrix{ith_sub, 1} = Result(index_s1_start:index_s2_start,:);
        final_matrix{ith_sub, 2} = Result(index_s2_start:index_s3_start,:);
        final_matrix{ith_sub, 3} = Result(index_s3_start:index_s4_start,:);
        final_matrix{ith_sub, 4} = Result(index_s4_start:end,:);
        [num2str(ith_sub)]
    end
end




%% extract trials
final_matrix2 = cell(19*4*40,4);
final_matrix = load('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/raw_eyedata_sub_by_sess_mri.mat');
final_matrix = final_matrix.final_matrix;

counter = 1;
for ith_sub = 8:26
    for ith_sess = 1:4
        curr_matrix = final_matrix{ith_sub, ith_sess};
        if ~isempty(curr_matrix)
            for ith_trial_index = 1:40
                index_start = find(cellfun(@(x) strcmp(x, ['session:',num2str(ith_sess),'trial:', num2str(ith_trial_index-1), 'screen:ITI']), curr_matrix(:,2)));
                index_end = find(cellfun(@(x) strcmp(x, ['session:',num2str(ith_sess),'trial:', num2str(ith_trial_index), 'screen:ITI']), curr_matrix(:,2)));
                
                if ~isempty(index_start) | ~isempty(index_end)
                    info_trial = curr_matrix{index_start, 2};
                    
                    final_matrix2{counter, 1} = ith_sub;
                    final_matrix2{counter, 2} = ith_sess;
                    final_matrix2{counter, 3} = info_trial;

                    if ith_trial_index == 40
                        final_matrix2{counter, 4} = curr_matrix(index_start:end, :);
                    else
                        final_matrix2{counter, 4} = curr_matrix(index_start:index_end, :);
                    end
                    counter = counter + 1;
                    ith_sub
                    ith_sess
                    ith_trial_index
                    
                end
            end
        end
    end
end






%% extract periods
final_matrix_trials = load('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/eyedata_formated_by_trial_mri.mat');
final_matrix_trials = final_matrix_trials.final_matrix2;
for ith_row = 1:length(final_matrix_trials(:,1))
    
%     trialid = strsplit(final_matrix_trials{ith_row,3}, 'TRIALID');
%     trialid = str2num(trialid{2})+1;
%     final_matrix_trials{ith_row,5} = trialid;
    
    curr_matrix = final_matrix_trials{ith_row, 4};
    
    %% check
    index_t_start = find(cellfun(@(x) contains(x, 'show_target'), curr_matrix(:,2)));
    index_t_end = find(cellfun(@(x) contains(x, 'show_response_cue'), curr_matrix(:,2)));
%     index_t_end = find(cellfun(@(x) contains(x, 'target_noise'), curr_matrix(:,2)));
    
    index_hnd_checker = find(cellfun(@(x) contains(x, 'show_idle_head_image'), curr_matrix(:,2)));
    
    final_matrix_trials{ith_row, 5} = curr_matrix(index_t_start:index_t_end, :);
    if index_hnd_checker
        final_matrix_trials{ith_row, 6} = 'hnd';
    end
end
% combine left, right, back info into matrix
for ith_row = 1:length(final_matrix_trials(:,1))
    curr_sub = final_matrix_trials{ith_row,1};
    curr_sess = final_matrix_trials{ith_row,2};
    curr_trial = final_matrix_trials{ith_row,3};
    curr_trial = str2num(cell2mat(extractBetween(curr_trial, 'trial:', 'screen:')))+1;
    
    [table_rawdata, table_time] = extract_behavioral_table_fmri ('/Users/bo/Documents/data_yuji_lab/data_fmri', curr_sub);
    table_rawdata_by_session = cell(1,17);
    temp_ith = find( cellfun(@(x,y) strcmp(x, num2str(curr_sess)), table_time{1,13}));
    for each_col = 1:length(table_rawdata(1,:))
        table_rawdata_by_session{1, each_col} = table_rawdata{1, each_col}(temp_ith);
    end

    final_matrix_trials{ith_row, 7} = table_rawdata_by_session{1,10}{curr_trial};
end
final_matrix_trials(find(cellfun(@(x) strcmp(x, 'hnd'), final_matrix_trials(:,6))),:)=[];





%% separate b l r c
% 11 no.fix.t 12 no.sac.t 13 no.fix.noise 14 no.sac.noise
the_matrix = load('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/eyedata_merged_with_beh_mri_4s.mat');
the_matrix = the_matrix.final_matrix_trials;

matrix_back = the_matrix(find(cellfun(@(x) strcmp(x, 'back'), the_matrix(:,7))),:);
matrix_left = the_matrix(find(cellfun(@(x) strcmp(x, 'left'), the_matrix(:,7))),:);
matrix_right = the_matrix(find(cellfun(@(x) strcmp(x, 'right'), the_matrix(:,7))),:);
final_matrix= {matrix_back, matrix_left, matrix_right};
for ith_con = 1:3
    cur_con = final_matrix{ith_con};
    for ith_row = 1:length(cur_con)
        tp_check_unit = cur_con{ith_row, 5};
        final_matrix{ith_con}{ith_row, 8} = length(find(cellfun(@(x) strcmp(x, 'ENDFIX'), tp_check_unit(:,3))));
        final_matrix{ith_con}{ith_row, 9} = length(find(cellfun(@(x) strcmp(x, 'ENDSACC'), tp_check_unit(:,3))));
    end
end


%% conpute mean for each sub and sess
the_matrix = load('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/eyedata_mean_for_pool_mri_4s.mat');
the_matrix = the_matrix.final_matrix;
num_sub = length(8:26);
for ith_con = 1:3
    
    cur_matrix = the_matrix{1,ith_con};
    output_matrix = nan(num_sub,2);
    
    t_fix_matrix = nan(num_sub, 4);
    t_sac_matrix = nan(num_sub, 4);
    for ith_sub = 8:26
        for ith_sess= 1:4
            t_fix_matrix(ith_sub, ith_sess) = nanmean(cell2mat(cur_matrix(find(cellfun(@(x,y) x==ith_sub & y==ith_sess, cur_matrix(:,1),cur_matrix(:,2))),8)));
            t_sac_matrix(ith_sub, ith_sess) = nanmean(cell2mat(cur_matrix(find(cellfun(@(x,y) x==ith_sub & y==ith_sess, cur_matrix(:,1),cur_matrix(:,2))),9)));
        end
    end
    t_fix_matrix_sess_mean = nanmean(t_fix_matrix,2);
    t_sac_matrix_sess_mean = nanmean(t_sac_matrix,2);
    output_matrix(:,1) = t_fix_matrix_sess_mean(8:end);
    output_matrix(:,2) = t_sac_matrix_sess_mean(8:end);
    
    if ith_con ==1
        fid = fopen('/Users/bo/Desktop/mri_eyedata_back_4s.txt','w');
        for n = 1:length(output_matrix(:,1))
            fprintf(fid,'%1.4f %1.4f\n', output_matrix(n,:));
        end
        fclose(fid);
    elseif ith_con ==2
        fid = fopen('/Users/bo/Desktop/mri_eyedata_left_4s.txt','w');
        for n = 1:length(output_matrix(:,1))
            fprintf(fid,'%1.4f %1.4f\n', output_matrix(n,:));
        end
        fclose(fid);
    elseif ith_con ==3
        fid = fopen('/Users/bo/Desktop/mri_eyedata_right_4s.txt','w');
        for n = 1:length(output_matrix(:,1))
            fprintf(fid,'%1.4f %1.4f\n', output_matrix(n,:));
        end
        fclose(fid);
    end
end








