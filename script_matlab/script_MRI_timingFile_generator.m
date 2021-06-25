



%% full trials
clear;
mypath = '/Users/bo/Documents/data_yuji_lab/data_fmri/';

for sub = 8:26
    mypath_rawData_record = [mypath,'/NAYA', num2str(sub),  '/sub_', num2str(sub), '_formal_rawdata_t.txt'];
    mypath_time_record =   [mypath,'/NAYA', num2str(sub),  '/sub_', num2str(sub), '_formal_Time_record_t.txt'];
    timeFile = fopen(mypath_time_record,'rt');
    table_time = textscan(timeFile,'%s %s %s %s %s %s %s %s %s %s %s %s %s'); %13
    fclose(timeFile);
    rawFile = fopen(mypath_rawData_record,'rt');
    table_rawdata = textscan(rawFile,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s'); % 18
    fclose(rawFile);

    % combine map & direction 11 12 13 14 21 22 23 24....
    durwalking=8;
    for sss = 1:4
        % trial
        session_index_trialtype = find(cellfun(@(x,y) strcmp(x,  num2str(sss)) & ~strcmp( y, 'NA'), table_time{1,13}, table_time{1,12}));
        timepoint_cell = table_time{1,5}(session_index_trialtype);
        timepoint_mat = cellfun(@(x) str2double(x), timepoint_cell);
        the_onset =  timepoint_mat.*0.001 - 6;
        % trialtype
        temp_cell = table_rawdata(:,[2,3]);
        cell_map = temp_cell{1};
        cell_direction = temp_cell{2};
        temp_mat = strcat(cell2mat(cell_map), cell2mat(cell_direction));
        temp_mat = temp_mat(session_index_trialtype,:);
        temp_mat = str2num(temp_mat);
        combined_mat = [the_onset temp_mat];
        session_ith = unique(combined_mat(:,2));
        for trialtype = 1:12
            trialtype_ith = find(combined_mat(:,2) == session_ith(trialtype));
            onset_list = combined_mat(trialtype_ith,1);
            temp_txt = nan(length(onset_list),3);
            temp_txt(:,1) = onset_list;
            temp_txt(:,2) = durwalking;
            temp_txt(:,3) = 1;
            fid = fopen([mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/map_direction', num2str(trialtype), '_', num2str(durwalking), 's.txt'],'w');
            for print_it = 1:length(temp_txt(:,1))
                fprintf(fid,'%f\t  %d\t  %d\n',temp_txt(print_it,:));
            end
            fclose(fid);
        end
    end
    

    % 36 trialtype walking 
    key_walking='walking';
    dur_walking=8;
    for sss = 1:4
        %trialtype
        session_index_trialtype = find(cellfun(@(x,y) strcmp( x,  num2str(sss)) & ~strcmp( y, 'NA'), table_time{1,13}, table_time{1,12}));
        timepoint_cell = table_time{1,5}(session_index_trialtype);
        timepoint_mat = cellfun(@(x) str2double(x), timepoint_cell);
        the_onset =  timepoint_mat.*0.001 - 6;
        for trialtype = 1:36
            temp_txt = nan(1,3);
            temp_txt(1,1) = the_onset(trialtype);
            temp_txt(1,2) = dur_walking;
            temp_txt(1,3) = 1;
            fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/', key_walking, num2str(dur_walking),'s_',num2str(trialtype) ,'.txt'],'w');
            fprintf(fid,'%f\t  %d\t  %d\n',temp_txt);
            fclose(fid);
        end
    end
    
    
    % 36 trialtype walking 4s 1
    key_walking='walking1_';
    dur_walking=4;
    for sss = 1:4
        %trialtype
        session_index_trialtype = find(cellfun(@(x,y) strcmp( x,  num2str(sss)) & ~strcmp( y, 'NA'), table_time{1,13}, table_time{1,12}));
        timepoint_cell = table_time{1,5}(session_index_trialtype);
        timepoint_mat = cellfun(@(x) str2double(x), timepoint_cell);
        the_onset =  timepoint_mat.*0.001 - 6;
        for trialtype = 1:36
            temp_txt = nan(1,3);
            temp_txt(1,1) = the_onset(trialtype);
            temp_txt(1,2) = dur_walking;
            temp_txt(1,3) = 1;
            fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/', key_walking, num2str(dur_walking),'s_',num2str(trialtype) ,'.txt'],'w');
            fprintf(fid,'%f\t  %d\t  %d\n',temp_txt);
            fclose(fid);
        end
    end
    
    % 36 trialtype walking 4s 2
    key_walking='walking2_';
    dur_walking=4;
    for sss = 1:4
        %trialtype
        session_index_trialtype = find(cellfun(@(x,y) strcmp( x,  num2str(sss)) & ~strcmp( y, 'NA'), table_time{1,13}, table_time{1,12}));
        timepoint_cell = table_time{1,5}(session_index_trialtype);
        timepoint_mat = cellfun(@(x) str2double(x), timepoint_cell);
        the_onset =  timepoint_mat.*0.001 - 6 + 4;
        for trialtype = 1:36
            temp_txt = nan(1,3);
            temp_txt(1,1) = the_onset(trialtype);
            temp_txt(1,2) = dur_walking;
            temp_txt(1,3) = 1;
            fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/', key_walking, num2str(dur_walking),'s_',num2str(trialtype) ,'.txt'],'w');
            fprintf(fid,'%f\t  %d\t  %d\n',temp_txt);
            fclose(fid);
        end
    end
    
    
    % 36 trialtype facing 
    key_facing='facing';
    dur_facing=4;
    for sss = 1:4
        %trialtype
        session_index_trialtype = find(cellfun(@(x,y) strcmp( x,  num2str(sss)) & ~strcmp( y, 'NA'), table_time{1,13}, table_time{1,12}));
        timepoint_cell = table_time{1,7}(session_index_trialtype);
        timepoint_mat = cellfun(@(x) str2double(x), timepoint_cell);
        the_onset =  timepoint_mat.*0.001 - 6;
        for trialtype = 1:36
            temp_txt = nan(1,3);
            temp_txt(1,1) = the_onset(trialtype);
            temp_txt(1,2) = dur_facing;
            temp_txt(1,3) = 1;
            fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/', key_facing, num2str(dur_facing),'s_',num2str(trialtype) ,'.txt'],'w');
            fprintf(fid,'%f\t  %d\t  %d\n',temp_txt);
            fclose(fid);
        end
    end
    
    % 36 trialtype facing 2s
    key_facing='facing';
    dur_facing=2;
    for sss = 1:4
        %trialtype
        session_index_trialtype = find(cellfun(@(x,y) strcmp( x,  num2str(sss)) & ~strcmp( y, 'NA'), table_time{1,13}, table_time{1,12}));
        timepoint_cell = table_time{1,7}(session_index_trialtype);
        timepoint_mat = cellfun(@(x) str2double(x), timepoint_cell);
        the_onset =  timepoint_mat.*0.001 - 6;
        for trialtype = 1:36
            temp_txt = nan(1,3);
            temp_txt(1,1) = the_onset(trialtype);
            temp_txt(1,2) = dur_facing;
            temp_txt(1,3) = 1;
            fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/', key_facing, num2str(dur_facing),'s_',num2str(trialtype) ,'.txt'],'w');
            fprintf(fid,'%f\t  %d\t  %d\n',temp_txt);
            fclose(fid);
        end
    end
    
    % 36 trialtype facing noise
    key_facing='facing_noise';
    dur_facing=2;
    for sss = 1:4
        %trialtype
        session_index_trialtype = find(cellfun(@(x,y) strcmp( x,  num2str(sss)) & ~strcmp( y, 'NA'), table_time{1,13}, table_time{1,12}));
        timepoint_cell = table_time{1,8}(session_index_trialtype);
        timepoint_mat = cellfun(@(x) str2double(x), timepoint_cell);
        the_onset =  timepoint_mat.*0.001 - 6;
        for trialtype = 1:36
            temp_txt = nan(1,3);
            temp_txt(1,1) = the_onset(trialtype);
            temp_txt(1,2) = dur_facing;
            temp_txt(1,3) = 1;
            fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/', key_facing, num2str(dur_facing),'s_',num2str(trialtype) ,'.txt'],'w');
            fprintf(fid,'%f\t  %d\t  %d\n',temp_txt);
            fclose(fid);
        end
    end
    
    % 36 trialtype facing 8s
    key_facing='facing';
    dur_facing=8;
    for sss = 1:4
        %trialtype
        session_index_trialtype = find(cellfun(@(x,y) strcmp( x,  num2str(sss)) & ~strcmp( y, 'NA'), table_time{1,13}, table_time{1,12}));
        timepoint_cell = table_time{1,7}(session_index_trialtype);
        timepoint_mat = cellfun(@(x) str2double(x), timepoint_cell);
        the_onset =  timepoint_mat.*0.001 - 6;
        for trialtype = 1:36
            temp_txt = nan(1,3);
            temp_txt(1,1) = the_onset(trialtype);
            temp_txt(1,2) = dur_facing;
            temp_txt(1,3) = 1;
            fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/', key_facing, num2str(dur_facing),'s_',num2str(trialtype) ,'.txt'],'w');
            fprintf(fid,'%f\t  %d\t  %d\n',temp_txt);
            fclose(fid);
        end
    end
    
    % 36 trialtype targeting 
    key_targeting='targeting';
    dur_targeting=4;
    for sss = 1:4
        %trialtype
        session_index_trialtype = find(cellfun(@(x,y) strcmp( x,  num2str(sss)) & ~strcmp( y, 'NA'), table_time{1,13}, table_time{1,12}));
        timepoint_cell = table_time{1,9}(session_index_trialtype);
        timepoint_mat = cellfun(@(x) str2double(x), timepoint_cell);
        the_onset =  timepoint_mat.*0.001 - 6;
        for trialtype = 1:36
            temp_txt = nan(1,3);
            temp_txt(1,1) = the_onset(trialtype);
            temp_txt(1,2) = dur_targeting;
            temp_txt(1,3) = 1;
            fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/', key_targeting, num2str(dur_targeting), 's_',num2str(trialtype) ,'.txt'],'w');
            fprintf(fid,'%f\t  %d\t  %d\n',temp_txt);
            fclose(fid);
        end
    end
    
    % 36 trialtype targeting 2s
    key_targeting='targeting';
    dur_targeting=2;
    for sss = 1:4
        %trialtype
        session_index_trialtype = find(cellfun(@(x,y) strcmp( x,  num2str(sss)) & ~strcmp( y, 'NA'), table_time{1,13}, table_time{1,12}));
        timepoint_cell = table_time{1,9}(session_index_trialtype);
        timepoint_mat = cellfun(@(x) str2double(x), timepoint_cell);
        the_onset =  timepoint_mat.*0.001 - 6;
        for trialtype = 1:36
            temp_txt = nan(1,3);
            temp_txt(1,1) = the_onset(trialtype);
            temp_txt(1,2) = dur_targeting;
            temp_txt(1,3) = 1;
            fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/', key_targeting, num2str(dur_targeting), 's_',num2str(trialtype) ,'.txt'],'w');
            fprintf(fid,'%f\t  %d\t  %d\n',temp_txt);
            fclose(fid);
        end
    end
    
    % 36 trialtype targeting
    key_targeting='targeting_noise';
    dur_targeting=2;
    for sss = 1:4
        %trialtype
        session_index_trialtype = find(cellfun(@(x,y) strcmp( x,  num2str(sss)) & ~strcmp( y, 'NA'), table_time{1,13}, table_time{1,12}));
        timepoint_cell = table_time{1,10}(session_index_trialtype);
        timepoint_mat = cellfun(@(x) str2double(x), timepoint_cell);
        the_onset =  timepoint_mat.*0.001 - 6;
        for trialtype = 1:36
            temp_txt = nan(1,3);
            temp_txt(1,1) = the_onset(trialtype);
            temp_txt(1,2) = dur_targeting;
            temp_txt(1,3) = 1;
            fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/', key_targeting, num2str(dur_targeting), 's_',num2str(trialtype) ,'.txt'],'w');
            fprintf(fid,'%f\t  %d\t  %d\n',temp_txt);
            fclose(fid);
        end
    end
    
    %% md trials
    temp_pool = cellfun(@(x) str2double(x), table_time{1,12});
    hm_ith = find(isnan(temp_pool));
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13}(hm_ith));
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 2}(hm_ith)) * 0.001 - 6 ;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 7}(hm_ith)) * 0.001 -6;
    temp_table = [hm_session hm_onset hm_end];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2) ;
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/md_trial.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    

    
    
    
    %map 
    hm_ith = find(isnan(cellfun(@(x) str2double(x), table_time{1,12})));
    map_col =  cellfun(@(x) str2double(x), table_rawdata{1, 2});
    map_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 5}) * 0.001 - 6;
    the_end = cellfun(@(x) str2double(x), table_time{1, 7}) * 0.001 - 6 ;
    hd_end = cellfun(@(x) str2double(x), table_time{1, 6}) * 0.001 - 6 ;
    hd_col = ones(length(map_col),1);
    hd_col(hm_ith) = 0;
    temp_table = [map_session the_onset the_end map_col hd_col];
    temp_table(find(temp_table(:,5)==0),3)=hd_end(find(temp_table(:,5)==0));
    map_pool = unique(temp_table(:,4));
    %% map1
    temp_table_map1 = temp_table(find(temp_table(:,4)==map_pool(1)),:);
    for sss = 1:4
        sess_ith = find(temp_table_map1(:,1)==sss);
        pre_onset =  temp_table_map1(sess_ith,2) ;
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_map1(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/map1.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% map2
    temp_table_map2 = temp_table(find(temp_table(:,4)==map_pool(2)),:);
    for sss = 1:4
        sess_ith = find(temp_table_map2(:,1)==sss);
        pre_onset =  temp_table_map2(sess_ith,2) ;
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_map2(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/map2.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% map3
    temp_table_map3 = temp_table(find(temp_table(:,4)==map_pool(3)),:);
    for sss = 1:4
        sess_ith = find(temp_table_map3(:,1)==sss);
        pre_onset =  temp_table_map3(sess_ith,2) ;
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_map3(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/map3.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    % direction
    hm_ith = find(isnan(cellfun(@(x) str2double(x), table_time{1,12})));
    direction_ith =  cellfun(@(x) str2double(x), table_rawdata{1, 3});
    direction_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 5}) * 0.001 - 6;
    the_end = cellfun(@(x) str2double(x), table_time{1, 7}) * 0.001 - 6 ;
    hd_end = cellfun(@(x) str2double(x), table_time{1, 6}) * 0.001 - 6 ;
    hd_col = ones(length(direction_ith),1);
    hd_col(hm_ith) = 0;
    temp_table = [direction_session the_onset the_end direction_ith hd_col];
    temp_table(find(temp_table(:,5)==0),3)=hd_end(find(temp_table(:,5)==0));
    %% direction1
    temp_table_direction1 = temp_table(find(temp_table(:,4)==1),:);
    for sss = 1:4
        sess_ith = find(temp_table_direction1(:,1)==sss);
        pre_onset =  temp_table_direction1(sess_ith,2) ;
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_direction1(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/direction1.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% direction2
    temp_table_direction2 = temp_table(find(temp_table(:,4)==2),:);
    for sss = 1:4
        sess_ith = find(temp_table_direction2(:,1)==sss);
        pre_onset =  temp_table_direction2(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_direction2(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/direction2.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% direction3
    temp_table_direction3 = temp_table(find(temp_table(:,4)==3),:);
    for sss = 1:4
        sess_ith = find(temp_table_direction3(:,1)==sss);
        pre_onset =  temp_table_direction3(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_direction3(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/direction3.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% direction4
    temp_table_direction4 = temp_table(find(temp_table(:,4)==4),:);
    for sss = 1:4
        sess_ith = find(temp_table_direction4(:,1)==sss);
        pre_onset =  temp_table_direction4(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_direction4(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/direction4.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    facingidle_ith =  cellfun(@(x) str2double(x), table_rawdata{1, 8});
    facingidle_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 7}) * 0.001 - 6;
    the_end = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001 - 6 ;
    temp_table = [facingidle_session the_onset the_end facingidle_ith];
    %% facingidle1
    temp_table_facingidle1 = temp_table(find(temp_table(:,4)==1),:);
    for sss = 1:4
        sess_ith = find(temp_table_facingidle1(:,1)==sss);
        pre_onset =  temp_table_facingidle1(sess_ith,2) ;
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_facingidle1(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facingidle1.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% facingidle2
    temp_table_facingidle2 = temp_table(find(temp_table(:,4)==2),:);
    for sss = 1:4
        sess_ith = find(temp_table_facingidle2(:,1)==sss);
        pre_onset =  temp_table_facingidle2(sess_ith,2) ;
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_facingidle2(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facingidle2.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% facingidle3
    temp_table_facingidle3 = temp_table(find(temp_table(:,4)==3),:);
    for sss = 1:4
        sess_ith = find(temp_table_facingidle3(:,1)==sss);
        pre_onset =  temp_table_facingidle3(sess_ith,2) ;
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_facingidle3(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facingidle3.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    facingidle_ith =  cellfun(@(x) str2double(x), table_rawdata{1, 8});
    facingidle_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 7}) * 0.001 - 6;
    the_end = cellfun(@(x) str2double(x), table_time{1, 8}) * 0.001 - 6 ;
    temp_table = [facingidle_session the_onset the_end facingidle_ith];
    %% facingidle1 2s
    temp_table_facingidle1 = temp_table(find(temp_table(:,4)==1),:);
    for sss = 1:4
        sess_ith = find(temp_table_facingidle1(:,1)==sss);
        pre_onset =  temp_table_facingidle1(sess_ith,2) ;
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_facingidle1(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facingidle1_2s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% facingidle2 2s
    temp_table_facingidle2 = temp_table(find(temp_table(:,4)==2),:);
    for sss = 1:4
        sess_ith = find(temp_table_facingidle2(:,1)==sss);
        pre_onset =  temp_table_facingidle2(sess_ith,2) ;
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_facingidle2(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facingidle2_2s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% facingidle3 2s
    temp_table_facingidle3 = temp_table(find(temp_table(:,4)==3),:);
    for sss = 1:4
        sess_ith = find(temp_table_facingidle3(:,1)==sss);
        pre_onset =  temp_table_facingidle3(sess_ith,2) ;
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_facingidle3(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facingidle3_2s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    

    %% target 4s
    target_ith =  cellfun(@(x) str2double(x), table_rawdata{1, 11});
    target_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001-6 ;
    the_end = cellfun(@(x) str2double(x), table_time{1, 11}) * 0.001-6 ;
    temp_table = [target_session the_onset the_end target_ith];
    %% target1
    temp_table_target1 = temp_table(find(temp_table(:,4)==1),:);
    for sss = 1:4
        sess_ith = find(temp_table_target1(:,1)==sss);
        pre_onset =  temp_table_target1(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_target1(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/target1.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% target2
    temp_table_target2 = temp_table(find(temp_table(:,4)==2),:);
    for sss = 1:4
        sess_ith = find(temp_table_target2(:,1)==sss);
        pre_onset =  temp_table_target2(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_target2(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/target2.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% target3
    temp_table_target3 = temp_table(find(temp_table(:,4)==3),:);
    for sss = 1:4
        sess_ith = find(temp_table_target3(:,1)==sss);
        pre_onset =  temp_table_target3(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_target3(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/target3.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    
    
    target_ith =  cellfun(@(x) str2double(x), table_rawdata{1, 11});
    target_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001-6 ;
    the_end = cellfun(@(x) str2double(x), table_time{1, 10}) * 0.001-6 ;
    temp_table = [target_session the_onset the_end target_ith];
    %% target1 2s
    temp_table_target1 = temp_table(find(temp_table(:,4)==1),:);
    for sss = 1:4
        sess_ith = find(temp_table_target1(:,1)==sss);
        pre_onset =  temp_table_target1(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_target1(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/target1_2s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% target2 2s
    temp_table_target2 = temp_table(find(temp_table(:,4)==2),:);
    for sss = 1:4
        sess_ith = find(temp_table_target2(:,1)==sss);
        pre_onset =  temp_table_target2(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_target2(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/target2_2s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% target3 2s
    temp_table_target3 = temp_table(find(temp_table(:,4)==3),:);
    for sss = 1:4
        sess_ith = find(temp_table_target3(:,1)==sss);
        pre_onset =  temp_table_target3(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_target3(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/target3_2s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    % cue
    hm_ith = find(isnan(cellfun(@(x) str2double(x), table_time{1,12})));
    cue_type =  cellfun(@(x) str2double(x), table_rawdata{1, 12});
    cue_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 11}) * 0.001 - 6;
    the_end = cellfun(@(x) str2double(x), table_time{1, 12}) * 0.001 - 6 ;
    temp_table = [cue_session the_onset the_end cue_type];
    %% cue1  123
    temp_table_cue_123 = temp_table(find(temp_table(:,4)==123),:);
    for sss = 1:4
        sess_ith = find(temp_table_cue_123(:,1)==sss);
        pre_onset =  temp_table_cue_123(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_cue_123(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        if(isempty(error_txt))
            'found empty'
        end
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/cue_123.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% cue2  132
    temp_table_cue_132 = temp_table(find(temp_table(:,4)==132),:);
    for sss = 1:4
        sess_ith = find(temp_table_cue_132(:,1)==sss);
        pre_onset =  temp_table_cue_132(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_cue_132(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        if(isempty(error_txt))
            'found empty'
        end
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/cue_132.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% cue3  213
    temp_table_cue_213 = temp_table(find(temp_table(:,4)==213),:);
    for sss = 1:4
        sess_ith = find(temp_table_cue_213(:,1)==sss);
        pre_onset =  temp_table_cue_213(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_cue_213(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        if(isempty(error_txt))
            'found empty'
        end
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/cue_213.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% cue4  231
    temp_table_cue_231 = temp_table(find(temp_table(:,4)==231),:);
    for sss = 1:4
        sess_ith = find(temp_table_cue_231(:,1)==sss);
        pre_onset =  temp_table_cue_231(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_cue_231(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        if(isempty(error_txt))
            'found empty'
        end
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/cue_231.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% cue5  312
    temp_table_cue_312 = temp_table(find(temp_table(:,4)==312),:);
    for sss = 1:4
        sess_ith = find(temp_table_cue_312(:,1)==sss);
        pre_onset =  temp_table_cue_312(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_cue_312(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        if(isempty(error_txt))
            'found empty'
        end
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/cue_312.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% cue6  321
    temp_table_cue_321 = temp_table(find(temp_table(:,4)==321),:);
    for sss = 1:4
        sess_ith = find(temp_table_cue_321(:,1)==sss);
        pre_onset =  temp_table_cue_321(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table_cue_321(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        if(isempty(error_txt))
            'found empty'
        end
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/cue_321.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    %% response stick
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 12}) * 0.001 - 6;
    md_onset = cellfun(@(x) str2double(x), table_time{1, 7}) * 0.001 - 6;
    temp_table = [hm_session hm_onset md_onset];
    md_ith = find(isnan(temp_table(:,2)));
    temp_table(md_ith,2)= temp_table(md_ith,3);
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        error_txt = [pre_onset zeros(length(pre_onset(:,1)),1) ones(length(pre_onset(:,1)),1)];
        fid = fopen([mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/response_stick.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    %% response period
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 11}) * 0.001-6 ;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 12}) * 0.001-6 ;
    temp_table = [hm_session hm_onset hm_end];
    temp_table(find(isnan(temp_table(:,3))),:)=[];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/response_period.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    %% tnoise+response
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 10}) * 0.001 - 6;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 12}) * 0.001 - 6;
    temp_table = [hm_session hm_onset hm_end    hm_end- hm_onset];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        error_txt(find(isnan(error_txt(:,1))),:)=[];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/tnoise_response.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end    
    
    % egocentric 4s
    egocentric_list =  table_rawdata{1, 10};
    egocentric_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001 - 6;
    the_end = cellfun(@(x) str2double(x), table_time{1, 11}) * 0.001 - 6 ;
    temp_table = cell(length(egocentric_list), 4);
    for n = 1:length(egocentric_list)
        temp_table{n,1} = egocentric_session(n);
        temp_table{n,2} = egocentric_list{n};
        temp_table{n,3} = the_onset(n);
        temp_table{n,4} = the_end(n);
    end
    %% ego1 left
    temp_table_ego1 = temp_table(cellfun(@(x) strcmp(x, 'left'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_ego1(:,1)));
        pre_onset =  cell2mat(temp_table_ego1(sess_ith,3));
        pre_end = cell2mat(temp_table_ego1(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/ego_left4s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% ego2 right
    temp_table_ego2 = temp_table(cellfun(@(x) strcmp(x, 'right'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_ego2(:,1)));
        pre_onset =  cell2mat(temp_table_ego2(sess_ith,3));
        pre_end = cell2mat(temp_table_ego2(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/ego_right4s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% ego3 back
    temp_table_ego3 = temp_table(cellfun(@(x) strcmp(x, 'back'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_ego3(:,1)));
        pre_onset =  cell2mat(temp_table_ego3(sess_ith,3));
        pre_end = cell2mat(temp_table_ego3(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/ego_back4s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    % egocentric 2s
    egocentric_list =  table_rawdata{1, 10};
    egocentric_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001 - 6;
    the_end = cellfun(@(x) str2double(x), table_time{1, 10}) * 0.001 - 6 ;
    temp_table = cell(length(egocentric_list), 4);
    for n = 1:length(egocentric_list)
        temp_table{n,1} = egocentric_session(n);
        temp_table{n,2} = egocentric_list{n};
        temp_table{n,3} = the_onset(n);
        temp_table{n,4} = the_end(n);
    end
    %% ego1 left
    temp_table_ego1 = temp_table(cellfun(@(x) strcmp(x, 'left'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_ego1(:,1)));
        pre_onset =  cell2mat(temp_table_ego1(sess_ith,3));
        pre_end = cell2mat(temp_table_ego1(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/ego_left2s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% ego2 right
    temp_table_ego2 = temp_table(cellfun(@(x) strcmp(x, 'right'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_ego2(:,1)));
        pre_onset =  cell2mat(temp_table_ego2(sess_ith,3));
        pre_end = cell2mat(temp_table_ego2(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/ego_right2s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% ego3 back
    temp_table_ego3 = temp_table(cellfun(@(x) strcmp(x, 'back'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_ego3(:,1)));
        pre_onset =  cell2mat(temp_table_ego3(sess_ith,3));
        pre_end = cell2mat(temp_table_ego3(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/ego_back2s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    % egocentric 2s noise period
    egocentric_list =  table_rawdata{1, 10};
    egocentric_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 10}) * 0.001 - 6;
    the_end = cellfun(@(x) str2double(x), table_time{1, 11}) * 0.001 - 6 ;
    temp_table = cell(length(egocentric_list), 4);
    for n = 1:length(egocentric_list)
        temp_table{n,1} = egocentric_session(n);
        temp_table{n,2} = egocentric_list{n};
        temp_table{n,3} = the_onset(n);
        temp_table{n,4} = the_end(n);
    end
    %% ego1 left
    temp_table_ego1 = temp_table(cellfun(@(x) strcmp(x, 'left'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_ego1(:,1)));
        pre_onset =  cell2mat(temp_table_ego1(sess_ith,3));
        pre_end = cell2mat(temp_table_ego1(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/ego_left2s_tn.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% ego2 right
    temp_table_ego2 = temp_table(cellfun(@(x) strcmp(x, 'right'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_ego2(:,1)));
        pre_onset =  cell2mat(temp_table_ego2(sess_ith,3));
        pre_end = cell2mat(temp_table_ego2(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/ego_right2s_tn.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% ego3 back
    temp_table_ego3 = temp_table(cellfun(@(x) strcmp(x, 'back'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_ego3(:,1)));
        pre_onset =  cell2mat(temp_table_ego3(sess_ith,3));
        pre_end = cell2mat(temp_table_ego3(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/ego_back2s_tn.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    
    
    % button 
    button_list =  table_rawdata{1, 17};
    button_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 11}) * 0.001 - 6;
    the_end = cellfun(@(x) str2double(x), table_time{1, 12}) * 0.001 - 6 ;
    temp_table = cell(length(button_list), 4);
    for n = 1:length(button_list)
        temp_table{n,1} = button_session(n);
        temp_table{n,2} = button_list{n};
        temp_table{n,3} = the_onset(n);
        temp_table{n,4} = the_end(n);
    end
    %% button 1
    temp_table_button1 = temp_table(cellfun(@(x) strcmp(x, '1'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_button1(:,1)));
        pre_onset =  cell2mat(temp_table_button1(sess_ith,3));
        pre_end = cell2mat(temp_table_button1(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/button1.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% button 2
    temp_table_button2 = temp_table(cellfun(@(x) strcmp(x, '2'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_button2(:,1)));
        pre_onset =  cell2mat(temp_table_button2(sess_ith,3));
        pre_end = cell2mat(temp_table_button2(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/button2.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% button 3
    temp_table_hand3 = temp_table(cellfun(@(x) strcmp(x, '3'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_hand3(:,1)));
        pre_onset =  cell2mat(temp_table_hand3(sess_ith,3));
        pre_end = cell2mat(temp_table_hand3(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/button3.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    % angle 
    angle_list =  table_rawdata{1, 15};
    angle_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 7}) * 0.001 - 6;
    the_end = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001 - 6 ;
    temp_table = cell(length(angle_list), 4);
    for n = 1:length(angle_list)
        temp_table{n,1} = angle_session(n);
        temp_table{n,2} = angle_list{n};
        temp_table{n,3} = the_onset(n);
        temp_table{n,4} = the_end(n);
    end
    %% angle1 left_45
    temp_table_angle1 = temp_table(cellfun(@(x) strcmp(x, 'left_45°'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_angle1(:,1)));
        pre_onset =  cell2mat(temp_table_angle1(sess_ith,3));
        pre_end = cell2mat(temp_table_angle1(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/angle1.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% angle2 left_135
    temp_table_angle2 = temp_table(cellfun(@(x) strcmp(x, 'left_135°'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_angle2(:,1)));
        pre_onset =  cell2mat(temp_table_angle2(sess_ith,3));
        pre_end = cell2mat(temp_table_angle2(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/angle2.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% angle3 right_45
    temp_table_angle3 = temp_table(cellfun(@(x) strcmp(x, 'right_45°'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_angle3(:,1)));
        pre_onset =  cell2mat(temp_table_angle3(sess_ith,3));
        pre_end = cell2mat(temp_table_angle3(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/angle3.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% angle4 right_135
    temp_table_angle4 = temp_table(cellfun(@(x) strcmp(x, 'right_135°'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_angle4(:,1)));
        pre_onset =  cell2mat(temp_table_angle4(sess_ith,3));
        pre_end = cell2mat(temp_table_angle4(sess_ith,4));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/angle4.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end    
    
    
    % facing period left / right     big / small angle
    angle_list =  table_rawdata{1, 15};
    angle_list_direction = cell(length(angle_list), 1);
    angle_list_angle = cell(length(angle_list), 1);
    for n = 1:length(angle_list)
        temp_split = split(angle_list{n,1}, '_');
        if strcmp(temp_split, 'NA')
            angle_list_direction{n,1} = NaN;
            angle_list_angle{n,1} = NaN;
        else
            angle_list_direction{n,1} = temp_split{1};
            angle_list_angle{n,1} = temp_split{2};
        end
    end
    angle_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    the_onset = cellfun(@(x) str2double(x), table_time{1, 7}) * 0.001 - 6;
    the_end = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001 - 6 ;
    temp_table = cell(length(angle_list), 4);
    for n = 1:length(angle_list)
        temp_table{n,1} = angle_session(n);
        temp_table{n,2} = angle_list_direction{n};
        temp_table{n,3} = angle_list_angle{n};
        temp_table{n,4} = the_onset(n);
        temp_table{n,5} = the_end(n);
    end
    %% f left
    temp_table_left = temp_table(cellfun(@(x) strcmp(x, 'left'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_left(:,1)));
        pre_onset =  cell2mat(temp_table_left(sess_ith,4));
        pre_end = cell2mat(temp_table_left(sess_ith,5));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facing_left.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% f right
    temp_table_right = temp_table(cellfun(@(x) strcmp(x, 'right'), temp_table(:,2)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_right(:,1)));
        pre_onset =  cell2mat(temp_table_right(sess_ith,4));
        pre_end = cell2mat(temp_table_right(sess_ith,5));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facing_right.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% 45 degree
    temp_table_angle45 = temp_table(cellfun(@(x) strcmp(x, '45°'), temp_table(:,3)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_angle45(:,1)));
        pre_onset =  cell2mat(temp_table_angle45(sess_ith,4));
        pre_end = cell2mat(temp_table_angle45(sess_ith,5));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facing_angle45.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    %% 135 degree
    temp_table_angle135 = temp_table(cellfun(@(x) strcmp(x, '135°'), temp_table(:,3)), :);
    for sss = 1:4
        sess_ith = find(cellfun(@(x) x==sss, temp_table_angle135(:,1)));
        pre_onset =  cell2mat(temp_table_angle135(sess_ith,4));
        pre_end = cell2mat(temp_table_angle135(sess_ith,5));
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facing_angle135.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    
    %% blank+baseline period
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 2}) * 0.001 - 6;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 5}) * 0.001 - 6;
    temp_table = [hm_session hm_onset hm_end    hm_end- hm_onset];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/blank+baseline_period.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end    
    

    
    %% baseline period
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 4}) * 0.001 - 6;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 5}) * 0.001 - 6;
    temp_table = [hm_session hm_onset hm_end    hm_end- hm_onset];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/baseline_period.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    %% walking period
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 5}) * 0.001 - 6;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 7}) * 0.001 - 6;
    temp_table = [hm_session hm_onset hm_end    hm_end- hm_onset];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/walking_period.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    %% facing period
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 7}) * 0.001 - 6;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001 - 6;
    temp_table = [hm_session hm_onset hm_end hm_end- hm_onset];
    md_ith = find(isnan(temp_table(:,3)));
    temp_table(md_ith,:) = [];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facing_period.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    %% facing period
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 5}) * 0.001 - 6;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001 - 6;
    temp_table = [hm_session hm_onset hm_end hm_end- hm_onset];
    md_ith = find(isnan(temp_table(:,3)));
    temp_table(md_ith,:) = [];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facing_period.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    %% walking + facing period 12s
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 5}) * 0.001 - 6;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001 - 6;
    temp_table = [hm_session hm_onset hm_end hm_end- hm_onset];
    md_ith = find(isnan(temp_table(:,3)));
    temp_table(md_ith,:) = [];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/walking_facing_12s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    %% facing period 8s
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 7}) * 0.001 - 6;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 11}) * 0.001 - 6;
    temp_table = [hm_session hm_onset hm_end hm_end- hm_onset];
    md_ith = find(isnan(temp_table(:,3)));
    temp_table(md_ith,:) = [];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([  mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/facing_period8s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    
    %% target period
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001 - 6;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 11}) * 0.001 - 6;
    temp_table = [hm_session hm_onset hm_end hm_end- hm_onset];
    md_ith = find(isnan(temp_table(:,3)));
    temp_table(md_ith,:) = [];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/target_period.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end
    
    %% target period 2s
    hm_session =  cellfun(@(x) str2double(x), table_time{1, 13});
    hm_onset = cellfun(@(x) str2double(x), table_time{1, 9}) * 0.001 - 6;
    hm_end = cellfun(@(x) str2double(x), table_time{1, 11}) * 0.001 - 6;
    temp_table = [hm_session hm_onset hm_end hm_end- hm_onset];
    md_ith = find(isnan(temp_table(:,3)));
    temp_table(md_ith,:) = [];
    for sss = 1:4
        sess_ith = find(temp_table(:,1)==sss);
        pre_onset =  temp_table(sess_ith,2);
        pre_onset(find(pre_onset<0)) = 0;
        pre_end = temp_table(sess_ith,3);
        dura = pre_end - pre_onset - 2;
        error_txt = [pre_onset dura ones(length(dura),1)];
        fid = fopen([mypath,  'NAYA',num2str(sub), '/fsl/s', num2str(sss), '/target_period2s.txt'],'w');
        for n = 1:length(error_txt(:,1))
            fprintf(fid,'%f\t  %f\t  %d\n',error_txt(n,:));
        end
        fclose(fid);
    end

    
end



