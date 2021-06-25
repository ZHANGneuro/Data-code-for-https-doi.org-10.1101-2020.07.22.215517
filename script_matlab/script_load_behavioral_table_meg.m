function [table_rawdata, table_time] = extract_behavioral_table_meg (global_path, sub)

path_rawData_record = [global_path,'/NAYA', num2str(sub),  '/sub_', num2str(sub), '_formal_rawdata_t.txt'];
path_time_record =   [global_path,'/NAYA', num2str(sub),  '/sub_', num2str(sub), '_formal_Time_record_t.txt'];


timeFile = fopen(path_time_record,'rt');
table_time = textscan(timeFile,'%s %s %s %s %s %s %s %s %s %s %s %s %s');
fclose(timeFile);
rawFile = fopen(path_rawData_record,'rt');
table_rawdata = textscan(rawFile,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s');
fclose(rawFile);




