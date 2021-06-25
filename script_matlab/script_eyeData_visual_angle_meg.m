

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


output_pool = cell(1,1);
counter = 1;
for ith = 1:length(final_matrix(:,1))
    cur_pool = final_matrix{ith, 6};
    target_pool = find(cellfun(@(x) ~isempty(x), cur_pool(:, 5)));
    for ith_target = 1:length(target_pool)
        output_pool{counter,1} = final_matrix{ith, 1};
        output_pool{counter,2} = final_matrix{ith, 2};
        ttarget_pool = cur_pool{target_pool(ith_target), 5};
        ititarget_pool = cur_pool{target_pool(ith_target), 4};
        output_pool{counter,3} = ttarget_pool;
        output_pool{counter,4} = cur_pool{target_pool(ith_target), 9};
        
        t_index_fix = find(cellfun(@(x) strcmp(x, 'ENDFIX'), ttarget_pool(:,3)));
        t_index_sca = find(cellfun(@(x) strcmp(x, 'ENDSACC'), ttarget_pool(:,3)));
        iti_index_fix = find(cellfun(@(x) strcmp(x, 'ENDFIX'), ititarget_pool(:,3)));
        iti_index_sca = find(cellfun(@(x) strcmp(x, 'ENDSACC'), ititarget_pool(:,3)));
        output_pool{counter, 5} = ititarget_pool(iti_index_fix, [6,7]);
        output_pool{counter, 6} = ititarget_pool(iti_index_sca, [6,7]);
        output_pool{counter, 7} = ttarget_pool(t_index_fix, [6,7]);
        output_pool{counter, 8} = ttarget_pool(t_index_sca, [6,7]);
        
        counter = counter + 1;
    end
end



%% separate b l r c
output_table2 = cell(1, 3); % left right back
ego_pool = {'left', 'right', 'back'};
counter = 1;
for ith_sub = 4:13
    for ith_sess = 1:4
        for ith_ego = 1:3
            temp_matrix = output_pool(find(cellfun(@(x, y, z) x==ith_sub & y==ith_sess & strcmp(z, ego_pool{ith_ego}), output_pool(:,1), output_pool(:,2), output_pool(:,4))), :);
            output_table2{counter, 1} = ith_sub;
            output_table2{counter, 2} = ith_sess;
            output_table2{counter, 3} = ego_pool{ith_ego};
            len_rows = length(temp_matrix(:,1));
            iti_fixation_coor = [];
            iti_saccade_axis = [];
            t_fixation_coor = [];
            t_saccade_axis = [];
            for ith = 1:len_rows
                if ~isempty(temp_matrix{ith, 5})
                    iti_fixation_coor = [iti_fixation_coor; temp_matrix{ith, 5}];
                end
                if ~isempty(temp_matrix{ith, 6})
                    iti_saccade_axis = [iti_saccade_axis; temp_matrix{ith, 6}];
                end
                if ~isempty(temp_matrix{ith, 7})
                    t_fixation_coor = [t_fixation_coor; temp_matrix{ith, 7}];
                end
                if ~isempty(temp_matrix{ith, 8})
                    t_saccade_axis = [t_saccade_axis; temp_matrix{ith, 8}];
                end
            end
            output_table2{counter, 4} = mean(cell2mat(iti_fixation_coor), 1);
            output_table2{counter, 5} = mean(cell2mat(iti_saccade_axis), 1);
            output_table2{counter, 6} = mean(cell2mat(t_fixation_coor), 1);
            output_table2{counter, 7} = mean(cell2mat(t_saccade_axis),1);
            counter = counter + 1;
        end
    end
end
output_table2(find(cellfun(@(x) isempty(x), output_table2(:,6))),:)=[];


output_table3 = cell(1, 3); % left right back
ego_pool = {'left', 'right', 'back'};
counter = 1;
for ith_sub = 4:13
    for ith_ego = 1:3
        sub_ego_pool = output_table2(find(cellfun(@(x,y) x==ith_sub & strcmp(y, ego_pool{ith_ego}), output_table2(:, 1), output_table2(:, 3))), :);
        
        len_row = length(sub_ego_pool(:,1));
        temp_fixaton = nan(len_row,2);
        temp_saccade = nan(len_row,2);
        for ith_temp = 1:len_row
            fix_x = sub_ego_pool{ith_temp,6}(1) -1024/2;
            fix_y = sub_ego_pool{ith_temp,6}(2) -768/2;
            sac_x = sub_ego_pool{ith_temp,6}(1) -1024/2;
            sac_y = sub_ego_pool{ith_temp,6}(2) -768/2;
            bl_fix_x = sub_ego_pool{ith_temp,4}(1) -1024/2;
            bl_fix_y = sub_ego_pool{ith_temp,4}(2) -768/2;
            bl_sac_x = sub_ego_pool{ith_temp,5}(1) -1024/2;
            bl_sac_y = sub_ego_pool{ith_temp,5}(2) -768/2;
            temp_fixaton(ith_temp, :) = [fix_x-bl_fix_x, fix_y-bl_fix_y];
            temp_saccade(ith_temp, :) = [sac_x-bl_sac_x, sac_y-bl_sac_y];
        end
        mean_fix = mean(temp_fixaton,1);
        mean_sac = mean(temp_saccade,1);
        output_table3{counter, 1} = ith_sub;
        output_table3{counter, 2} = ith_ego;
        output_table3{counter, 3} = mean_fix(1);
        output_table3{counter, 4} = mean_fix(2);
        output_table3{counter, 5} = mean_sac(1);
        output_table3{counter, 6} = mean_sac(2);
        counter = counter + 1;
    end
end

output_table3 = cell2mat(output_table3);

output_table3(:,[3,5]) = output_table3(:,[3,5]) * (419/1024);
output_table3(:,[4,6]) = output_table3(:,[4,6]) * (315/768);
output_table3(:,3:6) = radtodeg(atan(output_table3(:,3:6)/751));

% 1, 2, 3 left right back

fid = fopen('/Users/bo/Desktop/meg_fix_sac_visual_angle_back_1s_with_correction.txt','w');
for n = 1:length(output_table3(:,1))
    fprintf(fid,'%1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\n', output_table3(n,:));
end
fclose(fid);







