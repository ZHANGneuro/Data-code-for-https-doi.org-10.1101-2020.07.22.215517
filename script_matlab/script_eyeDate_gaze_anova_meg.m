

% compute visual angle
hight_limit_half = radtodeg(atan((315/2)/751));
width_limit_half = radtodeg(atan((419/2)/751));


sub_list = 4:13;
final_matrix = cell(length(sub_list), 4);
counter = 1;
for ith_sub = 4:13
    for ith_sess = 1:4
        file_path = ['/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/el_', num2str(ith_sub), '_s', num2str(ith_sess),'.edf'];
        if exist(file_path,'file')
            edf = Edf2Mat(file_path);
            Result_orig = edf.RawEdf.FEVENT;
            
            edf1 = edfmex(file_path);
            
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
            
            sample_x = edf.Samples.gx;
            sample_y = edf.Samples.gy;
            if sample_x(1,1)<-10000
                the_sample = [sample_x(:,2) sample_y(:,2)];
            elseif sample_x(1,2)<-10000
                the_sample = [sample_x(:,1) sample_y(:,1)];
            end
            
            index_start_clip = find(cellfun(@(x) strcmp(x, 'STARTPARSE'), Result(:,3)));
            Result = Result(index_start_clip:end, :);
            baseline_tp = Result{1,8};
            for ith_row = 1:length(Result(:,1))
                Result{ith_row,8} = Result{ith_row,8} - baseline_tp;
            end
            
            final_matrix{counter, 1} = ith_sub;
            final_matrix{counter, 2} = ith_sess;
            final_matrix{counter, 3} = Result;
            final_matrix{counter, 4} = the_sample;
            counter = counter + 1;
            [num2str(ith_sub), '_',  num2str(ith_sess)]
        end
    end
end


%
%% extract trials
final_matrix = load('/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/meg_eyeData_raw_for_position.mat');
final_matrix = final_matrix.final_matrix;
for ith_row = 1:length(final_matrix(:,1))
    curr_matrix = final_matrix{ith_row, 3};
    curr_matrix(find(cellfun(@(x) isempty(x), curr_matrix(:,2))),2) = {'nan'};
    index_array = find(cellfun(@(x) contains(x, 'TRIALID'), curr_matrix(:,2)));
    
    temp_matrix = cell(1,2);
    counter = 1;
    for ith_index = 1:length(index_array)
        temp_matrix{counter, 1} = curr_matrix{index_array(ith_index), 2}; % trial info
        if ith_index == length(index_array)
            temp_matrix{counter, 2} = curr_matrix(index_array(ith_index):length(curr_matrix(:,1)), :);
        else
            temp_matrix{counter, 2} = curr_matrix(index_array(ith_index):index_array(ith_index+1), :);
        end
        counter = counter + 1;
    end
    final_matrix{ith_row, 5} = temp_matrix;
end



% extract periods
% combine left, right, back info into matrix
for ith_row = 1:length(final_matrix(:,1))
    
    cur_sub = final_matrix{ith_row, 1};
    cur_sess = final_matrix{ith_row, 2};
    
    path_txt = ['/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_',num2str(cur_sub), '_s', num2str(cur_sess),'_rawdata.txt'];
    rawFile = fopen(path_txt,'rt');
    table_rawdata = textscan(rawFile,'%s %s %s %s %s %s %s %s %s %s %s %s %s');
    fclose(rawFile);
    
    cur_matrix = final_matrix{ith_row, 5};
    for ith_trial = 1:length(cur_matrix(:,1))
        
        trialid = strsplit(cur_matrix{ith_trial,1}, 'TRIALID');
        trialid = str2num(trialid{2})+1;
        cur_matrix{ith_trial,3} = trialid;
        
        event_list = cur_matrix{ith_trial,2};
        
        index_iti_start = find(cellfun(@(x) contains(x, 'ITI'), event_list(:,2)));
        index_iti_end = find(cellfun(@(x) contains(x, 'walking'), event_list(:,2)));
        
        index_t_start = find(cellfun(@(x) contains(x, 'targetingperson'), event_list(:,2)));
        index_t_end = find(cellfun(@(x) contains(x, 'noiseaftertargeting'), event_list(:,2)));
        
        index_tc_start = find(cellfun(@(x) contains(x, 'targetingcontrol'), event_list(:,2)));
        index_tc_end = find(cellfun(@(x) contains(x, 'noiseaftertargeting'), event_list(:,2)));
        
        index_tn_start = find(cellfun(@(x) contains(x, 'noiseaftertargeting'), event_list(:,2)));
        index_tn_end = find(cellfun(@(x) contains(x, 'cuedisplayed'), event_list(:,2)));
        
        index_hnd_checker = find(cellfun(@(x) contains(x, 'hndtrialheaddisplayed'), event_list(:,2)));
        
        cur_matrix{ith_trial, 4} = event_list(index_iti_start:index_iti_end, :);
        cur_matrix{ith_trial, 5} = event_list(index_t_start:index_t_end, :);
        cur_matrix{ith_trial, 6} = event_list(index_tc_start:index_tc_end, :);
        cur_matrix{ith_trial, 7} = event_list(index_tn_start:index_tn_end, :);
        if index_hnd_checker
            cur_matrix{ith_trial, 7} = 'hnd';
        end
        cur_matrix{ith_trial, 8} = table_rawdata{1,2}{ith_trial};
        cur_matrix{ith_trial, 9} = table_rawdata{1,11}{ith_trial};
    end
    cur_matrix(find(cellfun(@(x) strcmp(x, 'hnd')|strcmp(x, 'sc'), cur_matrix(:,8))),:)=[];
    final_matrix{ith_row, 6} = cur_matrix;
end


%% separate b l r c
for ith_row = 1:length(final_matrix(:,1))
    
    the_matrix = final_matrix{ith_row, 6};
    
    matrix_back = the_matrix(find(cellfun(@(x) strcmp(x, 'back'), the_matrix(:,9))),:);
    matrix_left = the_matrix(find(cellfun(@(x) strcmp(x, 'left'), the_matrix(:,9))),:);
    matrix_right = the_matrix(find(cellfun(@(x) strcmp(x, 'right'), the_matrix(:,9))),:);
    matrix_control = the_matrix(find(cellfun(@(x) strcmp(x, 'control'), the_matrix(:,9))),:);
    
    cell_back = cell(1,4);
    counter = 1;
    for ith = 1:length(matrix_back(:,1))
        cell_back{ith,1} = matrix_back{ith, 5}{1,8}; % t
        cell_back{ith,2} = matrix_back{ith, 5}{end,8};
        cell_back{ith,3} = matrix_back{ith, 7}{1,8}; % tnoise
        cell_back{ith,4} = matrix_back{ith, 7}{end,8};
        cell_back{ith,5} = matrix_back{ith, 4}{1,8}; % iti
        cell_back{ith,6} = matrix_back{ith, 4}{end,8};
        counter = counter + 1;
    end
    
    cell_left = cell(1,4);
    counter = 1;
    for ith = 1:length(matrix_left(:,1))
        cell_left{ith,1} = matrix_left{ith, 5}{1,8}; % t
        cell_left{ith,2} = matrix_left{ith, 5}{end,8};
        cell_left{ith,3} = matrix_left{ith, 7}{1,8}; % tnoise
        cell_left{ith,4} = matrix_left{ith, 7}{end,8};
        cell_left{ith,5} = matrix_left{ith, 4}{1,8}; % iti
        cell_left{ith,6} = matrix_left{ith, 4}{end,8};
        counter = counter + 1;
    end
    
    cell_right = cell(1,4);
    counter = 1;
    for ith = 1:length(matrix_right(:,1))
        cell_right{ith,1} = matrix_right{ith, 5}{1,8}; % t
        cell_right{ith,2} = matrix_right{ith, 5}{end,8};
        cell_right{ith,3} = matrix_right{ith, 7}{1,8}; % tnoise
        cell_right{ith,4} = matrix_right{ith, 7}{end,8};
        cell_right{ith,5} = matrix_right{ith, 4}{1,8}; % iti
        cell_right{ith,6} = matrix_right{ith, 4}{end,8};
        counter = counter + 1;
    end
    
    cell_control = cell(1,4);
    counter = 1;
    for ith = 1:length(matrix_control(:,1))
        cell_control{ith,1} = matrix_control{ith, 6}{1,8}; % t
        cell_control{ith,2} = matrix_control{ith, 6}{end,8};
        cell_control{ith,3} = matrix_control{ith, 7}{1,8}; % tnoise
        cell_control{ith,4} = matrix_control{ith, 7}{end,8};
        cell_control{ith,5} = matrix_control{ith, 4}{1,8}; % iti
        cell_control{ith,6} = matrix_control{ith, 4}{end,8};
        counter = counter + 1;
    end
    
    final_matrix{ith_row, 7} = cell_back;
    final_matrix{ith_row, 8} = cell_left;
    final_matrix{ith_row, 9} = cell_right;
    final_matrix{ith_row, 10} = cell_control;
end

selec_dura = 4; % col ith
screen_matrix_back_t = cell(1, 6);
screen_matrix_left_t = cell(1, 6);
screen_matrix_right_t = cell(1, 6);
ith_row_back = 1;
ith_row_left = 1;
ith_row_right = 1;
for ith_sub = 4:13
    for ith_sess = 1:4
        
        target_index = find(cellfun(@(x,y) x==ith_sub & y==ith_sess, final_matrix(:,1), final_matrix(:,2)));
        if ~isempty(target_index)
            % each 6 col: t tend tn tnend iti itend
            matrix_back = final_matrix{target_index, 7};
            matrix_left = final_matrix{target_index, 8};
            matrix_right = final_matrix{target_index, 9};
            matrix_control = final_matrix{target_index, 10};
            matrix_position = final_matrix{target_index, 4};
            
            %back-t
            bl_position = [];
            for ith = 1:length(matrix_back(:,1))
                bl_position = [bl_position; matrix_position(matrix_back{ith,5}:matrix_back{ith,6}, :)];
            end
            bl_x = nanmean(bl_position(:,1))- 1024/2;
            bl_y = nanmean(bl_position(:,2))- 768/2;
            
            for ith = 1:length(matrix_back(:,1))
                cur_position = matrix_position(matrix_back{ith,1}:matrix_back{ith,selec_dura}, :);
                cur_position(find(isnan(cur_position(:,1))),:)=[];
                cur_position(:,1) = cur_position(:,1)-1024/2;
                cur_position(:,2) = cur_position(:,2)-768/2;
                cur_position(:,1) = cur_position(:,1)-bl_x;
                cur_position(:,2) = cur_position(:,2)-bl_y;
                cur_position(:,1) = cur_position(:,1) * (419/1024);
                cur_position(:,2) = cur_position(:,2) * (315/768);
                cur_position = radtodeg(atan(cur_position/751));
                screen_matrix_back_t{ith_row_back, 1} = ith_sub;
                screen_matrix_back_t{ith_row_back, 2} = ith_sess;
                screen_matrix_back_t{ith_row_back, 3} = ith;
                screen_matrix_back_t{ith_row_back, 4} = cur_position;
                screen_matrix_back_t{ith_row_back, 5} = mean(cur_position(:,1));
                screen_matrix_back_t{ith_row_back, 6} = mean(cur_position(:,2));
                ith_row_back = ith_row_back + 1;
            end
            
            %left-t
            for ith = 1:length(matrix_left(:,1))
                bl_position = [bl_position; matrix_position(matrix_left{ith,5}:matrix_left{ith,6}, :)];
            end
            bl_x = nanmean(bl_position(:,1))- 1024/2;
            bl_y = nanmean(bl_position(:,2))- 768/2;
            for ith = 1:length(matrix_left(:,1))
                cur_position = matrix_position(matrix_left{ith,1}:matrix_left{ith,selec_dura}, :);
                cur_position(find(isnan(cur_position(:,1))),:)=[];
                cur_position(:,1) = cur_position(:,1)-1024/2;
                cur_position(:,2) = cur_position(:,2)-768/2;
                cur_position(:,1) = cur_position(:,1)-bl_x;
                cur_position(:,2) = cur_position(:,2)-bl_y;
                cur_position(:,1) = cur_position(:,1) * (419/1024);
                cur_position(:,2) = cur_position(:,2) * (315/768);
                cur_position = radtodeg(atan(cur_position/751));
                screen_matrix_left_t{ith_row_left, 1} = ith_sub;
                screen_matrix_left_t{ith_row_left, 2} = ith_sess;
                screen_matrix_left_t{ith_row_left, 3} = ith;
                screen_matrix_left_t{ith_row_left, 4} = cur_position;
                screen_matrix_left_t{ith_row_left, 5} = mean(cur_position(:,1));
                screen_matrix_left_t{ith_row_left, 6} = mean(cur_position(:,2));
                ith_row_left = ith_row_left + 1;
            end
            
            %right-t
            bl_position = [];
            for ith = 1:length(matrix_right(:,1))
                bl_position = [bl_position; matrix_position(matrix_right{ith,5}:matrix_right{ith,6}, :)];
            end
            bl_x = nanmean(bl_position(:,1))- 1024/2;
            bl_y = nanmean(bl_position(:,2))- 768/2;
            for ith = 1:length(matrix_right(:,1))
                cur_position = matrix_position(matrix_right{ith,1}:matrix_right{ith,selec_dura}, :);
                cur_position(find(isnan(cur_position(:,1))),:)=[];
                cur_position(:,1) = cur_position(:,1)-1024/2;
                cur_position(:,2) = cur_position(:,2)-768/2;
                cur_position(:,1) = cur_position(:,1)-bl_x;
                cur_position(:,2) = cur_position(:,2)-bl_y;
                cur_position(:,1) = cur_position(:,1) * (419/1024);
                cur_position(:,2) = cur_position(:,2) * (315/768);
                cur_position = radtodeg(atan(cur_position/751));
                screen_matrix_right_t{ith_row_right, 1} = ith_sub;
                screen_matrix_right_t{ith_row_right, 2} = ith_sess;
                screen_matrix_right_t{ith_row_right, 3} = ith;
                screen_matrix_right_t{ith_row_right, 4} = cur_position;
                screen_matrix_right_t{ith_row_right, 5} = mean(cur_position(:,1));
                screen_matrix_right_t{ith_row_right, 6} = mean(cur_position(:,2));
                ith_row_right = ith_row_right + 1;
            end
        end
    end
end

screen_matrix_back_t(find(cellfun(@(x) isnan(x), screen_matrix_back_t(:,5))),:)=[];
screen_matrix_left_t(find(cellfun(@(x) isnan(x), screen_matrix_left_t(:,5))),:)=[];
screen_matrix_right_t(find(cellfun(@(x) isnan(x), screen_matrix_right_t(:,5))),:)=[];

screen_matrix_back_t = cell2mat(screen_matrix_back_t(:, [1,2,3,5,6]));
screen_matrix_left_t = cell2mat(screen_matrix_left_t(:, [1,2,3,5,6]));
screen_matrix_right_t = cell2mat(screen_matrix_right_t(:, [1,2,3,5,6]));

%% plot heatmap for position
data_pool = {screen_matrix_back_t, screen_matrix_left_t, screen_matrix_right_t};
label_pool = {'back', 'left','right'};
for ith_matrix = 1:length(data_pool)
    temp_matrix = data_pool{ith_matrix};
    fid = fopen(['/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/gaze_', label_pool{ith_matrix},'_bar_anova_4s_meg_bl_corrected.txt'],'w');
    for n = 1:length(temp_matrix(:,1))
        fprintf(fid, '%3.4f\t %3.4f\t %3.4f\t %3.4f\t %3.4f\n', temp_matrix(n,:));
    end
    fclose(fid);
end






