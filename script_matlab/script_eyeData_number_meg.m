

final_matrix = cell(length(2:13), 4);
for ith_sub = 4:13
    for ith_sess = 1:4
        file_path = ['/Users/bo/Documents/data_yuji_lab/data_MEG/data_behav/el_', num2str(ith_sub), '_s', num2str(ith_sess),'.edf'];
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
            
            final_matrix{ith_sub, ith_sess} = Result;
            [num2str(ith_sub), '_',  num2str(ith_sess)]
        end
        % sample_x = edf.Samples.gx;
        % sample_y = edf.Samples.gy;
    end
end




%% extract trials
final_matrix2 = cell(13,4);
final_matrix_cheaker = cell(13*4,3);
final_matrix = load('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData/raw_eyedata_sub_by_sess.mat');
final_matrix = final_matrix.final_matrix;

counter = 1;
counter_cheaker = 1;
for ith_sub = 1:13
    for ith_sess = 1:4
        
        curr_matrix = final_matrix{ith_sub, ith_sess};
        
        if ~isempty(curr_matrix)
            
            %% clearning cell format
            curr_matrix(find(cellfun(@(x) isempty(x), curr_matrix(:,2))),2) = {'nan'};
            
            %%
            index_array = find(cellfun(@(x) contains(x, 'TRIALID'), curr_matrix(:,2)));
            final_matrix_cheaker{counter_cheaker, 1} = ith_sub;
            final_matrix_cheaker{counter_cheaker, 2} = ith_sess;
            final_matrix_cheaker{counter_cheaker, 3} = length(index_array);
            counter_cheaker = counter_cheaker + 1;
            
            for ith_index = 1:length(index_array)
                
                info_trial = curr_matrix{index_array(ith_index), 2};
                
                final_matrix2{counter, 1} = ith_sub;
                final_matrix2{counter, 2} = ith_sess;
                final_matrix2{counter, 3} = info_trial;
                
                if ith_index == length(index_array)
                    final_matrix2{counter, 4} = curr_matrix(index_array(ith_index):curr_matrix{end,1}, :);
                else
                    final_matrix2{counter, 4} = curr_matrix(index_array(ith_index):index_array(ith_index+1), :);
                end
                counter = counter + 1;
            end
        end
    end
end






%% extract periods
final_matrix_trials = load('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/eyedata_formated_by_trial_meg.mat');
final_matrix_trials = final_matrix_trials.final_matrix2;
for ith_row = 1:length(final_matrix_trials(:,1))
    
    trialid = strsplit(final_matrix_trials{ith_row,3}, 'TRIALID');
    trialid = str2num(trialid{2})+1;
    final_matrix_trials{ith_row,5} = trialid;
    
    curr_matrix = final_matrix_trials{ith_row, 4};
    
    %% check
    index_t_start = find(cellfun(@(x) contains(x, 'targetingperson'), curr_matrix(:,2)));
    index_t_end = find(cellfun(@(x) contains(x, 'noiseaftertargeting'), curr_matrix(:,2)));
    
    index_tc_start = find(cellfun(@(x) contains(x, 'targetingcontrol'), curr_matrix(:,2)));
    index_tc_end = find(cellfun(@(x) contains(x, 'noiseaftertargeting'), curr_matrix(:,2)));
    
    index_tn_start = find(cellfun(@(x) contains(x, 'noiseaftertargeting'), curr_matrix(:,2)));s
    index_tn_end = find(cellfun(@(x) contains(x, 'cuedisplayed'), curr_matrix(:,2)));
    
    index_hnd_checker = find(cellfun(@(x) contains(x, 'hndtrialheaddisplayed'), curr_matrix(:,2)));
    
    final_matrix_trials{ith_row, 6} = curr_matrix(index_t_start:index_t_end, :); % targeting period
    final_matrix_trials{ith_row, 7} = curr_matrix(index_tc_start:index_tc_end, :); % control trial
    final_matrix_trials{ith_row, 8} = curr_matrix(index_tn_start:index_tn_end, :); % targeting noise period
    if index_hnd_checker
        final_matrix_trials{ith_row, 8} = 'hnd';
    end
end
% combine left, right, back info into matrix
for ith_row = 1:length(final_matrix_trials(:,1))
    curr_sub = final_matrix_trials{ith_row,1};
    curr_sess = final_matrix_trials{ith_row,2};
    cur_trial = final_matrix_trials{ith_row,5};
    path_txt = ['/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/sub_',num2str(curr_sub), '_s', num2str(curr_sess),'_rawdata.txt'];
    rawFile = fopen(path_txt,'rt');
    table_rawdata = textscan(rawFile,'%s %s %s %s %s %s %s %s %s %s %s %s %s');
    fclose(rawFile);
  
    final_matrix_trials{ith_row, 9} = table_rawdata{1,2}{cur_trial};
    final_matrix_trials{ith_row, 10} = table_rawdata{1,11}{cur_trial};
end
final_matrix_trials(find(cellfun(@(x) strcmp(x, 'hnd')|strcmp(x, 'sc'), final_matrix_trials(:,9))),:)=[];



%% separate b l r c
% 11 no.fix.t 12 no.sac.t 13 no.fix.noise 14 no.sac.noise
the_matrix = load('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/eyedata_merged_with_beh_meg_both1&2s.mat');
the_matrix = the_matrix.final_matrix_trials;

matrix_back = the_matrix(find(cellfun(@(x) strcmp(x, 'back'), the_matrix(:,10))),:);
matrix_left = the_matrix(find(cellfun(@(x) strcmp(x, 'left'), the_matrix(:,10))),:);
matrix_right = the_matrix(find(cellfun(@(x) strcmp(x, 'right'), the_matrix(:,10))),:);
matrix_control = the_matrix(find(cellfun(@(x) strcmp(x, 'control'), the_matrix(:,10))),:);
final_matrix= {matrix_back, matrix_left, matrix_right, matrix_control};
for ith_con = 1:4
    cur_con = final_matrix{ith_con};
    for ith_row = 1:length(cur_con)
        if ith_con==4
            tp_check_unit = cur_con{ith_row, 7}; % control trial
        else
            tp_check_unit = cur_con{ith_row, 6}; % targeting period
        end
        tn_check_unit = cur_con{ith_row, 8}; % targeting noise period
        final_matrix{ith_con}{ith_row, 11} = length(find(cellfun(@(x) strcmp(x, 'ENDFIX'), tp_check_unit(:,3)))); % numfix tp
        final_matrix{ith_con}{ith_row, 12} = length(find(cellfun(@(x) strcmp(x, 'ENDSACC'), tp_check_unit(:,3)))); % numsac tp
        final_matrix{ith_con}{ith_row, 13} = length(find(cellfun(@(x) strcmp(x, 'ENDFIX'), tn_check_unit(:,3)))); % numfix tnp
        final_matrix{ith_con}{ith_row, 14} = length(find(cellfun(@(x) strcmp(x, 'ENDSACC'), tn_check_unit(:,3)))); % numsac tnp
    end
end


%% conpute mean for each sub and sess col 12-saccade t, 14-saccade tn
the_matrix = load('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/eyedata_mean_for_pool_meg_both1&2s.mat');
the_matrix = the_matrix.final_matrix;
num_sub = length(4:13);
for ith_con = 1:4
    
    cur_matrix = the_matrix{1,ith_con};
    output_matrix = nan(num_sub,4);
    
    t_fix_matrix = nan(num_sub, 4);
    t_sac_matrix = nan(num_sub, 4);
    tn_fix_matrix = nan(num_sub, 4);
    tn_sac_matrix = nan(num_sub, 4);
    for ith_sub = 4:13
        for ith_sess= 1:4
            t_fix_matrix(ith_sub, ith_sess) = nanmean(cell2mat(cur_matrix(find(cellfun(@(x,y) x==ith_sub & y==ith_sess, cur_matrix(:,1),cur_matrix(:,2))),11)));
            t_sac_matrix(ith_sub, ith_sess) = nanmean(cell2mat(cur_matrix(find(cellfun(@(x,y) x==ith_sub & y==ith_sess, cur_matrix(:,1),cur_matrix(:,2))),12)));
            tn_fix_matrix(ith_sub, ith_sess) = nanmean(cell2mat(cur_matrix(find(cellfun(@(x,y) x==ith_sub & y==ith_sess, cur_matrix(:,1),cur_matrix(:,2))),13)));
            tn_sac_matrix(ith_sub, ith_sess) = nanmean(cell2mat(cur_matrix(find(cellfun(@(x,y) x==ith_sub & y==ith_sess, cur_matrix(:,1),cur_matrix(:,2))),14)));
        end
    end
    t_fix_matrix_sess_mean = nanmean(t_fix_matrix,2);
    t_sac_matrix_sess_mean = nanmean(t_sac_matrix,2);
    tn_fix_matrix_sess_mean = nanmean(tn_fix_matrix,2);
    tn_sac_matrix_sess_mean = nanmean(tn_sac_matrix,2);
    output_matrix(:,1) = t_fix_matrix_sess_mean(4:end); % tp fix
    output_matrix(:,2) = t_sac_matrix_sess_mean(4:end); % tp sac
    output_matrix(:,3) = tn_fix_matrix_sess_mean(4:end); % tnp fix
    output_matrix(:,4) = tn_sac_matrix_sess_mean(4:end); % tnp sac
    
    if ith_con ==1
        fid = fopen('/Users/bo/Desktop/eyedata_back.txt','w');
        for n = 1:length(output_matrix(:,1))
            fprintf(fid,'%1.4f %1.4f %1.4f %1.4f\n', output_matrix(n,:));
        end
        fclose(fid);
    elseif ith_con ==2
        fid = fopen('/Users/bo/Desktop/eyedata_left.txt','w');
        for n = 1:length(output_matrix(:,1))
            fprintf(fid,'%1.4f %1.4f %1.4f %1.4f\n', output_matrix(n,:));
        end
        fclose(fid);
    elseif ith_con ==3
        fid = fopen('/Users/bo/Desktop/eyedata_right.txt','w');
        for n = 1:length(output_matrix(:,1))
            fprintf(fid,'%1.4f %1.4f %1.4f %1.4f\n', output_matrix(n,:));
        end
        fclose(fid);
    elseif ith_con ==4
        fid = fopen('/Users/bo/Desktop/eyedata_control.txt','w');
        for n = 1:length(output_matrix(:,1))
            fprintf(fid,'%1.4f %1.4f %1.4f %1.4f\n', output_matrix(n,:));
        end
        fclose(fid);
    end
end








