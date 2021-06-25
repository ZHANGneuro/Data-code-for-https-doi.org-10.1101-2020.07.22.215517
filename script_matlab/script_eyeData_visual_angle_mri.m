

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
%     index_t_end = find(cellfun(@(x) contains(x, 'show_response_cue'), curr_matrix(:,2)));
    index_t_end = find(cellfun(@(x) contains(x, 'target_noise'), curr_matrix(:,2)));
    
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
the_matrix = load('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/eyedata_merged_with_beh_mri.mat');
the_matrix = the_matrix.final_matrix_trials;

for ith = 1:length(the_matrix(:,1))
    tp_check_unit = the_matrix{ith, 5};
    index_fix = find(cellfun(@(x) strcmp(x, 'ENDFIX'), tp_check_unit(:,3)));
    index_sca = find(cellfun(@(x) strcmp(x, 'ENDSACC'), tp_check_unit(:,3)));
    
    the_matrix{ith, 8} = length(index_fix);
    the_matrix{ith, 9} = length(index_sca);
    the_matrix{ith, 10} = tp_check_unit(index_fix, [6,7]);
    the_matrix{ith, 11} = tp_check_unit(index_sca, [6,7]);
end


output_table = cell(1, 3); % left right back
ego_pool = {'left', 'right', 'back'};
counter = 1;
for ith_sub = 8:26
    for ith_sess = 1:4
        for ith_ego = 1:3
            temp_matrix = the_matrix(find(cellfun(@(x, y, z) x==ith_sub & y==ith_sess & strcmp(z, ego_pool{ith_ego}), the_matrix(:,1), the_matrix(:,2), the_matrix(:,7))), :);
            
            output_table{counter, 1} = ith_sub;
            output_table{counter, 2} = ith_sess;
            output_table{counter, 3} = ego_pool{ith_ego};
            
            len_rows = length(temp_matrix(:,1));
            fixation_coor = [];
            saccade_axis = [];
            for ith = 1:len_rows
                if ~isempty(temp_matrix{ith, 10})
                    fixation_coor = [fixation_coor; temp_matrix{ith, 10}];
                    saccade_axis = [saccade_axis; temp_matrix{ith, 11}];
                end
            end
            output_table{counter, 4} = mean(cell2mat(fixation_coor), 1);
            output_table{counter, 5} = mean(cell2mat(saccade_axis),1);
            output_table{counter, 6} = cell2mat(fixation_coor);
            output_table{counter, 7} = cell2mat(saccade_axis);
            counter = counter + 1;
        end
    end
end
output_table(find(cellfun(@(x) isempty(x), output_table(:,4))),:)=[];

output_table2 = cell(1, 3); % left right back
ego_pool = {'left', 'right', 'back'};
counter = 1;
for ith_sub = 8:26
    for ith_ego = 1:3
        sub_ego_pool = output_table(find(cellfun(@(x,y) x==ith_sub & strcmp(y, ego_pool{ith_ego}), output_table(:, 1), output_table(:, 3))), :);
        
        len_row = length(sub_ego_pool(:,1));
        temp_fixaton = nan(len_row,2);
        temp_saccade = nan(len_row,2);
        for ith_temp = 1:len_row
            temp_fixaton(ith_temp, :) = sub_ego_pool{ith_temp,4};
            temp_saccade(ith_temp, :) = sub_ego_pool{ith_temp,5};
        end
        mean_fix = mean(temp_fixaton,1);
        mean_sac = mean(temp_saccade,1);
        output_table2{counter, 1} = ith_sub;
        output_table2{counter, 2} = ith_ego;
        output_table2{counter, 3} = mean_fix(1);
        output_table2{counter, 4} = mean_fix(2);
        output_table2{counter, 5} = mean_sac(1);
        output_table2{counter, 6} = mean_sac(2);
        counter = counter + 1;
    end
end

output_table2(find(cellfun(@(x) isnan(x), output_table2(:,3))),:)=[];
output_table2 = cell2mat(output_table2);
output_table2(:,[3,5]) = output_table2(:,[3,5])-1024/2;
output_table2(:,[4,6]) = output_table2(:,[4,6])-768/2;
output_table2(find(output_table2(:,3)>=1024/2 | output_table2(:,5)>=1024/2 | output_table2(:,3)<=-1024/2 | ...
    output_table2(:,5)<=-1024/2| output_table2(:,4)<=-768/2 | output_table2(:,6)<=-768/2| output_table2(:,4)>=768/2 | ...
    output_table2(:,6)>=768/2),:)=[]; % carefull to add or not here
output_table2(:,[3,5]) = output_table2(:,[3,5]) * (419/1024);
output_table2(:,[4,6]) = output_table2(:,[4,6]) * (315/768);
output_table2(:,3:6) = radtodeg(atan(output_table2(:,3:6)/751));

% 1, 2, 3 left right back


fid = fopen('/Users/bo/Desktop/mri_fix_sac_visual_angle_back_2s.txt','w');
for n = 1:length(output_table2(:,1))
    fprintf(fid,'%1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\n', output_table2(n,:));
end
fclose(fid);






