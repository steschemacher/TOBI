function [basin] = fun_load_basin(SQ_name,data_name,path,path_output)

% function to load basins from WaSiM

speicher = importdata(strcat(path_output,'\',data_name));
SQ = importdata(strcat(path,'\',SQ_name));
ShQ(:,1) = SQ(:,1);
ShQ(:,3) = SQ(:,2);

basin_rule = speicher.data(:,:);
S0 = speicher.data(1,7);
basin = struct('ShQ',ShQ,'basin_rule',basin_rule,'S0',S0,'abstraction',[]);