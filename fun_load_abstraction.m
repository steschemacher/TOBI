function [abstr] = fun_load_abstraction(abstr_name,path_output,wellen_start_line,wellen_end_line)

% function to load abstractions from WaSiM

abstraction_data = importdata(strcat(path_output,'\',abstr_name));

abstr = struct('ShQ',[],'basin_rule',[],'S0',[],'abstraction',abstraction_data.data(wellen_start_line:wellen_end_line,5));
