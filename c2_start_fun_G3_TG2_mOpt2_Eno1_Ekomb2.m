Gebiet = 3;
TG = 2;
mOpt = 2;
E_no = 1;
E_komb = 2;
E = [1;2;3;4;5;6;7];
path = 'PFAD';
path1 = 'PFAD CODE';
path2 = 'PFAD RESULTS';
addpath(path);
addpath(path1);
addpath(path2);
cd(path2);
diary('G3_TG2_mOpt2_Eno1_Ekomb2.txt');
disp('diary started.')
try 
   disp('Gebiet 3, TG 2 gestartet:')
   c2_fun_start_runs(Gebiet,TG,mOpt,E_no,E_komb,E,path2)
catch ME 
   cd(strcat(path2,'\Sim_overview')) 
   while exist('open_status_142.txt','file')==2 
       pause(5*rand(1)) 
   end 
   save('open_status_142.txt') 
   load('status_142.mat'); 
       status_142(2) = -1; 
       status_142(4:5) = now; 
   save('status_142.mat','status_142'); 
   delete('open_status_142.txt') 
   msg = getReport(ME); 
   disp(msg); 
   

 disp('Execution not possible: G3_TG2_mOpt2_Eno1_Ekomb2') 
   diary off 
   pause(5) 
   quit
end