%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Berechnung von Retention und Abflussrouting im Flussnetz       %%%
%%%                          << TOBas >>                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input:
%   - Information about the grid
%       - grid-name
%       - distance between gages [m]
%   - Tanalys Output directory (with all gages)
%       - PUR (gages)
%   - Tanalys Output directory (with one gage at the basin outlet)
%       - DGM (digital elevation model)
%       - DEP (river depth)
%       - EZG (catchment + subcatchments)
%       - FZS (flow times) + only with the 1st gage!!!!
%       - FLK (flow directions)
%       - gridasci.exe
%   - WaSiM Output directory (with all gages)
%       - qgko*.stat
%       - qges*.stat
%       - gwin*.stat
% required functions:
%   - Flusspunkt.m (class with propertis and functions of the river flowpoints)
%   - Pegel.m (class with properties and functions of the gages)
%   - fun_routing.m (recursive function for discharge routing)
%   - fun_retention.m (function for basin retention including outlet optimization)
%   - fun_load_basin.m (function to import basin information from WaSiM)
%   - fun_load_abstraction.m (function to import abstractions from WaSiM)
%   - fun_include_basins.m (function import existing basins in the routing structure)
%   - fun_include_abstraction.m (function to include abstractions in the routing st.)
%   - fun_S_gesamt.m (function to sum up the retention volume in a basin)
% output: efficacy of the basin at every gage


function f2__wirksamkeitsanalyse(Gebiet,TG,...
    E_no_opt,Beckenkombination_speichern,...
    mOpt,opt_parameter,...
    path_allg,path_becken,...
    method,beobachtungspunkte)



%% Daten einlesen

% Pfaddefinitionen
path_data = strcat(path_allg,'\Inputdaten\',sprintf('G%d_TG%d_input',Gebiet,TG));
if mOpt==1
    E_opt_komb = E_no_opt(1,3);
    E_opt = E_no_opt(1,2);
    path_results = strcat(path_allg,'\Results\',sprintf('G%d_TG%d_mOpt%d_komb%d',Gebiet,TG,mOpt,abs(E_opt_komb)));
else
    path_results = strcat(path_allg,'\Results\',sprintf('G%d_TG%d_mOpt%d',Gebiet,TG,mOpt));
end
cd(path_allg)
if exist(path_results,'dir')~=7
    mkdir(path_results)
end


cd(path_data)
% räumliche Daten
load('river.mat')
% Zeitreihen
load('ereignis_real.mat')    % ereignisse, Q_mean, HHQ
% Routing-Struktur
load('routing_struct.mat')

% Reale Beckenstandorte
cd(path_becken)
if mOpt==1
    if E_no_opt(1,3)==-1
        load(sprintf('G%d_beckenpunkte_select1.mat',Gebiet))
        beckenpunkte_select2 = beckenpunkte_select1;
    else
        load(sprintf('G%d_beckenpunkte_select2.mat',Gebiet)) 
    end
elseif mOpt==10
    load(sprintf('G%d_beckenpunkte_select1.mat',Gebiet))
    beckenpunkte_select2 = beckenpunkte_select1;
    if E_no_opt(1,3)==-3
        load(sprintf('G%d_TG%d_beckenpunkte_LfU_fiktiv.mat',Gebiet,TG))
        beckenpunkte_select2 = beckenpunkte_LfU(opt_parameter);
        if real_fiktiv==2
            cd(path_results)
            path_results = strcat(path_results,'\fiktiv');
            if exist(path_results,'dir')~=7
            mkdir(path_results);
            end
        end
    end
%     if E_no_opt(1,3)==4        
        path_results = strcat(path_results,'\ganglinien_superposition_v2');         
        mkdir(path_results)
        pfad_original_data = strcat(path_allg,'\Results\',sprintf('G%d_TG%d_mOpt2',Gebiet,TG));
%     end
else
    load(sprintf('G%d_beckenpunkte_select2.mat',Gebiet))
end


%% Definition of Gage Properties

[pegel,max_gage] = f2_1__define_gage_properties(grids,river,A_ezg,river_points,Q_mean,HHQ_Ereignisse,routing);

clear river river_points HHQ_Ereignisse
clear routing grids


%% Selection of basin combinations
cd(path_results)
if mOpt==1
    if E_no_opt(1,2)==1
        if exist(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG),'file')==2
            if exist(sprintf('Beckenkombination%d_Ekomb%d.mat',E_opt,abs(E_opt_komb)),'file')==2
                load(sprintf('Beckenkombination%d_Ekomb%d.mat',E_opt,abs(E_opt_komb)));
            else
                load(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG))
            end
        else
            [Beckenkombination,pegel_runde,pegel_becken] = f2_2__select_basins(pegel,Gebiet,TG,mOpt,beckenpunkte_select2,opt_parameter);
            save(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG),'Beckenkombination','pegel_runde','pegel_becken')
        end
    else
        while exist(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG),'file')~=2
            disp('Waiting for basin combination of event 1 to read.')
            pause(60)
        end
        if exist(sprintf('Beckenkombination%d_Ekomb%d.mat',E_opt,abs(E_opt_komb)),'file')==2
            load(sprintf('Beckenkombination%d_Ekomb%d.mat',E_opt,abs(E_opt_komb)));
        else
            load(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG))
        end
    end
elseif mOpt==2
    if exist(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG),'file')==2
        load(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG))
    else
        [Beckenkombination,pegel_runde,pegel_becken] = f2_2__select_basins(pegel,Gebiet,TG,mOpt,beckenpunkte_select2,opt_parameter);
        save(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG),'Beckenkombination','pegel_runde','pegel_becken')
    end
elseif mOpt~=13
    
    cd(pfad_original_data)
    files = dir(strcat(sprintf('Eno%d_Eopt%d_kombE%d',E_no_opt(1,:)),'*','_opt_beckenkombination_',sprintf('%4.2f',opt_parameter(1)),'*','.mat'));
    load(files.name,'beckenkombination_opt');
    cd(path_results)
     [Beckenkombination,pegel_runde,pegel_becken] = f2_2__select_basins(pegel,Gebiet,TG,mOpt,beckenpunkte_select2,beckenkombination_opt);
        save(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG),'Beckenkombination','pegel_runde','pegel_becken')
else
     cd(path_results)
     [Beckenkombination,pegel_runde,pegel_becken] = f2_2__select_basins(pegel,Gebiet,TG,mOpt,beckenpunkte_select2,[]);
        save(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG),'Beckenkombination','pegel_runde','pegel_becken')

end

% Anzahl der Simulationsruns (Einzelspeicherdateien)
if mOpt==1
    anzahl_sim = size(E_no_opt,1)*length(Beckenkombination);
elseif mOpt==2
    anzahl_sim = size(E_no_opt,1);
end
n_anzahl_sim = 0;
%% Berechnung der Wirksamkeiten

start_global = now;
% E_index = 0;
for E_index = 1:size(E_no_opt,1)
    
    if mOpt==2
        n_anzahl_sim = n_anzahl_sim+1;
    end
    
%     E_index = E_index+1;
    clearvars -except E_index start_global Beckenkombination pegel_runde pegel_becken pegel max_gage path_data path_results Gebiet TG ...
        Beckenkombination_speichern mOpt opt_parameter path_allg path_becken method beobachtungspunkte E_no_opt beckenpunkte_select2 n_anzahl_sim anzahl_sim ...
        pfad_original_data

    E_no = E_no_opt(E_index,1);
    E_opt = E_no_opt(E_index,2);
    E_opt_komb = E_no_opt(E_index,3);
    
    % bestehende Beckenkombinationen und Drosselöffnungen einladen
    if mOpt~=2

        % Drosselweite berechnen (0) oder einlesen (1)
        if E_no==E_opt
            feste_Drosselweite=0;
        else
            feste_Drosselweite=1;
        end
    end
    
    % Ganglinien + Richards-Datei einlesen
    cd(path_data)
    load('ereignis_real.mat')
    load('richards.mat')
    Ereignis_neu(E_no) = Ereignis(E_no);
    Ereignis = Ereignis_neu;
%     Ereignis = rmfield(Ereignis,{'qgko_spec','qges_date','qges_spec','gwin_date','gwin_spec'});
    richards_neu(E_no,:) = richards(E_no,:);
    richards = richards_neu;
    clear richards_neu Ereignis_neu
    clear qgko_date
    qgko_date = Ereignis(E_no).qgko_date;
    
    
    
    if mOpt==2
        cd(path_results)
        cd(path_results)
        datei = dir(sprintf('Eno%d_Eopt%d_kombE%d_*noBasins_opt_beckenkombination_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,opt_parameter));
        if isempty(datei)==0        
%         if exist(sprintf('Eno%d_Eopt%d_kombE%d_%02dnoBasins_opt_beckenkombination_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,size(Beckenkombination,1),opt_parameter))
            disp(strcat(sprintf('Eno%d_Eopt%d_kombE%d_%02dnoBasins_opt_beckenkombination_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,size(Beckenkombination,1),opt_parameter),' already exists.'))
                        
%             clear Beckenkombination
%             datei = dir(sprintf('Eno%d_Eopt%d_kombE%d_*noBasins_opt_beckenkombination_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,opt_parameter));
%             load(datei(end).name)
%             Beckenkombination(1).becken(:,1) = beckenkombination_opt;
%             Beckenkombination(1).becken_real(:,1) = beckenkombination_opt;
                        
%             cd(strcat(path_allg,'\Results\Sim_overview'))
%             
%             while exist('open_characteristics.txt','file')==2
%                 pause(5*rand(1))
%             end
%             save('open_characteristics.txt','Gebiet','-ascii')
%             load('characteristics.mat')
%             delete('open_characteristics.txt')
%             zeile = characteristics(ismember(characteristics(:,2:6),[Gebiet,TG,mOpt,1,E_opt_komb],'rows'),1);
%             
%             while exist(sprintf('open_status_%d.txt',zeile),'file')==2
%                 pause(5*rand(1))
%             end
%             save(sprintf('open_status_%d.txt',zeile),'Gebiet','-ascii')
%             load(sprintf('status_%d.mat',zeile))
%             status = eval(sprintf('status_%d', zeile));
%             status(2) = 1;    % status: läuft oder nicht
%             status(3) = Beckenkombination_speichern(2)/Beckenkombination_speichern(1)*100+n_anzahl_sim/(anzahl_sim*Beckenkombination_speichern(1))*100; % Fotschritt 8%)
%             status(5) = now;    % Zeit
%             tmp = ['status_', num2str(zeile),'= status'];
%             eval(tmp)
%             clear tmp
%             save(sprintf('status_%d.mat',zeile),sprintf('status_%d',zeile))
%             delete(sprintf('open_status_%d.txt',zeile))
%             
            continue
        else
            disp(strcat('Optimization of simulation ',sprintf('Eno%d_Eopt%d_kombE%d_parameter_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f',E_no,E_opt,E_opt_komb,opt_parameter),' starts now.'))
        end
        if exist(sprintf('Eoptkomb%d_beckenkombination_drosselweite_alleopt.mat',E_opt_komb))
            load(sprintf('Eoptkomb%d_beckenkombination_drosselweite_alleopt.mat',E_opt_komb))
            if E_no==E_opt && E_opt==E_opt_komb
                n_opt_komb = length(Beckenkombination_final)+1;
                
                cd(path_results)
                datei = dir(sprintf('Eno%d_Eopt%d_kombE%d_*noBasins_opt_beckenkombination_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,opt_parameter));
                if isempty(datei)==0
%                 if exist(sprintf('Eno%d_Eopt%d_kombE%d_%02dnoBasins_opt_beckenkombination_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,size(Beckenkombination,1),opt_parameter))
                    clear Beckenkombination
                    datei = dir(sprintf('Eno%d_Eopt%d_kombE%d_*noBasins_opt_beckenkombination_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,opt_parameter));
                    load(datei(end).name)                    
                    Beckenkombination(1).becken(:,1) = beckenkombination_opt;
                    Beckenkombination(1).becken_real(:,1) = beckenkombination_opt;
                else 
                    clear Beckenkombination
                    if exist(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG),'file')==2
                        load(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG))
                    end
                end
            else
                opt_parameter_tmp = vertcat(Beckenkombination_final.opt_parameter);
                n_opt_komb = find(ismember(opt_parameter_tmp,opt_parameter,'rows')==1);% findrows length(Beckenkombination_final);
                clear Beckenkombination
                
                Beckenkombination(1).becken(:,1) = Beckenkombination_final(n_opt_komb).becken_real;
                Beckenkombination(1).becken_real(:,1) = Beckenkombination_final(n_opt_komb).becken_real;
                if E_no~=E_opt
                    Beckenkombination(1).drosselweite(:,1,E_opt) = Beckenkombination_final(n_opt_komb).drosselweite(:,E_opt);
                end
                pegel_runde{1} = TG;
            end
        else
            n_opt_komb = 1;
        end
        if E_no==E_opt
            feste_Drosselweite=0;
        else
            feste_Drosselweite=1;
        end  
    elseif mOpt==10
        cd(pfad_original_data)
        load(sprintf('Eoptkomb%d_beckenkombination_drosselweite_alleopt.mat',E_opt_komb))
        if E_no==E_opt && E_opt==E_opt_komb
            n_opt_komb = length(Beckenkombination_final)+1;
            
            cd(pfad_original_data)
            datei = dir(sprintf('Eno%d_Eopt%d_kombE%d_*noBasins_opt_beckenkombination_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,opt_parameter));
            if isempty(datei)==0
                %                 if exist(sprintf('Eno%d_Eopt%d_kombE%d_%02dnoBasins_opt_beckenkombination_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,size(Beckenkombination,1),opt_parameter))
                clear Beckenkombination
                datei = dir(sprintf('Eno%d_Eopt%d_kombE%d_*noBasins_opt_beckenkombination_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,opt_parameter));
                load(datei(end).name)
                if isempty(beckenkombination_opt)==1
                    files = dir(sprintf('Eno%d_Eopt%d_kombE%d_*noBasins_opt_lokal_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,opt_parameter));
                    load(files.name,'lokal_opt');
                    beckenkombination_opt = lokal_opt(:,2)';                    
                end
                Beckenkombination(1).becken(:,1) = beckenkombination_opt;
                Beckenkombination(1).becken_real(:,1) = beckenkombination_opt;
            else
                clear Beckenkombination
                if exist(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG),'file')==2
                    load(sprintf('BeckenkombinationBasis_G%d_TG%d.mat',Gebiet,TG))
                end
            end
        else
            opt_parameter_tmp = vertcat(Beckenkombination_final.opt_parameter);
            n_opt_komb = find(ismember(opt_parameter_tmp,opt_parameter,'rows')==1);% findrows length(Beckenkombination_final);
            clear Beckenkombination
            
            Beckenkombination(1).becken(:,1) = Beckenkombination_final(n_opt_komb).becken_real;
            Beckenkombination(1).becken_real(:,1) = Beckenkombination_final(n_opt_komb).becken_real;
            if E_no~=E_opt
                Beckenkombination(1).drosselweite(:,1,E_opt) = Beckenkombination_final(n_opt_komb).drosselweite(:,E_opt);
            end
            pegel_runde{1} = TG;
        end
        cd(path_results)
    end
    if size(Beckenkombination)>1
        keyboard
    end
    
    for i_no_of_basins = 1:size(Beckenkombination,1)
        if mOpt==1
        n_anzahl_sim = n_anzahl_sim+1;
        end
%         no_of_basins = beckenanzahlen(i_no_of_basins);
        
        clear Wirksamkeit_Beckenkombination info_sim info_beckenkombination becken_in_TG becken_pos_real
        clear lokal regional wirksamkeitsverlauf daten_ges zielfunktion optimierungsfunktion
        zeit_start = now;
        
        % Überprüfung ob die Simulationen bereits durchgeführt wurden
        if mOpt==1
        cd(path_results)
%         fprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_zielfunktion.mat',E_no,E_opt,abs(E_opt_komb),i_no_of_basins)
        if exist(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_zielfunktion.mat',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'file')==2
            disp(strcat(sprintf('Eno%d_Eopt%d_Ekomb%d__%02dnoBasins_zielfunktion.mat',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),' already exists.'))
            
            cd(strcat(path_allg,'\Results\Sim_overview'))
            
            while exist('open_characteristics.txt','file')==2
                pause(5*rand(1))
            end
            save('open_characteristics.txt','Gebiet','-ascii')
            load('characteristics.mat')
            delete('open_characteristics.txt')
            zeile = characteristics(ismember(characteristics(:,2:6),[Gebiet,TG,mOpt,E_no,E_opt_komb],'rows'),1);
            
            while exist(sprintf('open_status_%d.txt',zeile),'file')==2
                pause(5*rand(1))
            end
            save(sprintf('open_status_%d.txt',zeile),'Gebiet','-ascii')
            load(sprintf('status_%d.mat',zeile))
            status = eval(sprintf('status_%d', zeile));
            status(2) = 1;    % status: läuft oder nicht
            status(3) = n_anzahl_sim/anzahl_sim*100; % Fotschritt 8%)
            status(5) = now;    % Zeit
            tmp = ['status_', num2str(zeile),'= status'];
            eval(tmp)
            clear tmp
            save(sprintf('status_%d.mat',zeile),sprintf('status_%d',zeile))
            delete(sprintf('open_status_%d.txt',zeile))
            
%             cd(strcat(path_allg,'\Results'))
%             
%             count = 0; err_count = 0;
%             while count == err_count
%                 try     load('simulation_overview.mat','characteristics')
%                 catch;   err_count = err_count + 1;     end
%                 count = count + 1; pause(3)
%             end
%             
%             zeile = characteristics(ismember(characteristics(:,2:6),[Gebiet,TG,mOpt,E_no,E_opt_komb],'rows'),1);
%             count = 0; err_count = 0;
%             while count == err_count
%                 try     load('simulation_overview.mat',sprintf('status_%d',zeile))
%                 catch;   err_count = err_count + 1;     end
%                 count = count + 1; pause(3)
%             end
%             
% %             load('simulation_overview.mat','characteristics')
% %             zeile = find(ismember(characteristics(:,2:6),[Gebiet,TG,mOpt,E_no,E_opt_komb],'rows')==1);
% %             load('simulation_overview.mat',sprintf('status_%d',zeile))
%             status = eval(sprintf('status_%d', zeile));
%             status(2) = 1;    % status: läuft oder nicht
%             status(3) = n_anzahl_sim/anzahl_sim*100; % Fotschritt 8%)
%             status(5) = now;    % Zeit
%             tmp = ['status_', num2str(zeile),'= status'];
%             eval(tmp)
%             clear tmp
% %             assignin ('base',['status_' num2str(zeile)], status);
%             save('simulation_overview','-append',sprintf('status_%d',zeile))
            
%             m = matfile('simulation_overview.mat','Writable',true);
%             zeile = find(m.ismember(simulation_overview(:,1:5),[Gebiet,TG,mOpt,E_opt,E_opt_komb],'rows')==1);
%             m.simulation_overview(zeile,7) = n_anzahl_sim/anzahl_sim*100;
%             m.simulation_overview(zeile,8) = 1;
%             m.simulation_overview(zeile,10) = now;
%             save('simulation_overview.mat','simulation_overview')
            continue
        end
        fprintf('Eno%d_Eopt%d_Ekomb%d_%d: ',E_no,E_opt,abs(E_opt_komb),i_no_of_basins)
        end

        
        % Pegelnummern aus vorab-Simulation ersetzen
        if mOpt ~=2
            Beckenkombination(i_no_of_basins,1).becken_real = Beckenkombination(i_no_of_basins,1).becken;
        elseif mOpt == 2
            if E_no==E_opt && E_opt==E_opt_komb
                for runde = 1:size(Beckenkombination(i_no_of_basins,1).becken,2)
                    becken_pos_real = Beckenkombination(i_no_of_basins,1).becken(:,runde);
                    becken_pos_real = becken_pos_real(isnan(becken_pos_real)==0);
                    if i_no_of_basins>1 && isempty(becken_pos_real)==0
                        for pegel_beckenkomb_vorab = [becken_pos_real(find(becken_pos_real<0))]
                            becken_pos_real = setdiff(becken_pos_real,pegel_beckenkomb_vorab);
                            becken_pos_real = vertcat(becken_pos_real,[beckenkomb_vorab_opt{(-1)*pegel_beckenkomb_vorab}]');
                        end
                    end
                    Beckenkombination(i_no_of_basins,1).becken_real(:,runde) = 0;
                    Beckenkombination(i_no_of_basins,1).becken_real(1:length(becken_pos_real),runde) = becken_pos_real;
                end
                Beckenkombination(i_no_of_basins,1).becken_real = unique(Beckenkombination(i_no_of_basins,1).becken_real','rows')';
            end
        end
        max_becken_pos = size(Beckenkombination(i_no_of_basins,1).becken_real,1);
        
        for runde = 1:size(Beckenkombination(i_no_of_basins,1).becken_real,2)          
            
            %% Preparation of gages and basin positions
            
            %%% Delete values of previous run and add new discharges
            [pegel] = f2_3__clean_Pegel(pegel,E_no,Ereignis,richards,max_gage);
            
            %%% Beckenummern und Pegelnummern der Kombination bestimmen
            becken_pos_real = Beckenkombination(i_no_of_basins,1).becken_real(:,runde);
            becken_pos_real = becken_pos_real(becken_pos_real~=0);
            
            becken_pos = nan(size(becken_pos_real,1),2);
            if mOpt~=13
                for i = 1:size(becken_pos_real,1)
                    zeile = find([beckenpunkte_select2.no_becken] == becken_pos_real(i));
                    becken_pos(i,1) = beckenpunkte_select2(zeile).no_pegel;  % Pegelnummer
                    becken_pos(i,2) = becken_pos_real(i);   % Beckennummer
                end
            else
                becken_pos(:,1) = becken_pos_real;
                becken_pos(:,2) = becken_pos_real;
            end
            
            if isempty(becken_pos)==0
                
                %%% Beckeneigenschaften festlegen (Beckendefinition + Drosselweite)
                [pegel,startpegel] = f2_4__set_basin_properties(pegel,beckenpunkte_select2,becken_pos,Beckenkombination(i_no_of_basins),E_no,E_opt,feste_Drosselweite,runde,existing_structures,Gebiet,mOpt);
                
                
                %% Discharge Routing
                %%% basin retention, routing and wave superposition
                f2_5__retention_routing(pegel,startpegel,method);
                
                %%% save the diameter of the outflow
                if mOpt~=13 && feste_Drosselweite ==0
                    for j = 1:size(becken_pos,1)
                        no_pegel = becken_pos(j,1);
                        zeile_becken = find([Beckenkombination(i_no_of_basins).becken_real(:,runde)]==becken_pos(j,2));
                        if zeile_becken>size(becken_pos,1)
                            disp('Fehler bei der Becken-Drosselweite-Zuordnung.')
                        end
                        spalte = find([pegel(no_pegel).beckenparameter.beckennummer]==becken_pos(j,2));
                        Beckenkombination(i_no_of_basins).drosselweite(zeile_becken,runde,E_opt) = pegel(no_pegel).beckenparameter(spalte).A_L;
                    end
                end
                
            else
                for p = 1:length(pegel)
                    pegel(p).inflow{1}(:,2) = pegel(p).inflow_wasim(:,3);
                    pegel(p).outflow{1}(:,2) = pegel(p).inflow_wasim(:,3);
                end
            end
            
            
            %% Calculate efficiencies of basin combination
            [info_sim(runde,:),info_beckenkombination(runde,:,:),zielfunktion(runde,:,:),daten_ges(runde,:,:),lokal(runde,:,:),regional(runde,:,:)] = ...
                f2_6__results_efficiency(qgko_date,becken_pos,max_becken_pos,beobachtungspunkte,beckenpunkte_select2,pegel,mOpt);
            
            fprintf('%d,',runde)
            
            if mOpt==10
            elseif mOpt==13
                keyboard
            elseif mOpt==11 || mOpt==12
                [pegel] = f2_7__validierung_WaSiM(Gebiet,TG,E_no,becken_pos,beckenpunkte_select2,pegel,runde,i_no_of_basins,path_allg);
                keyboard
            end
            
        end
        
        %% Zielfunktion
        if mOpt==2
            
            for p = pegel_runde{i_no_of_basins}
                f_sa = opt_parameter(1);
                f_fuellung = opt_parameter(2);
                f_anzahl = opt_parameter(3);
                f_flaeche = opt_parameter(4);
                f_dammvol = opt_parameter(5);
                
                max_SA = max(squeeze(zielfunktion(:,p,1)));
                min_SA = min(squeeze(zielfunktion(:,p,1)));
                if max_SA~=min_SA
                    optimierungsfunktion(:,1) = f_sa*((squeeze(zielfunktion(:,p,4)))+...
                        (1-(squeeze(zielfunktion(:,p,1))-min_SA)/(max_SA-min_SA)))...
                        + f_fuellung * (1 - squeeze(zielfunktion(:,p,17)))...
                        + f_anzahl * squeeze(zielfunktion(:,p,18))...
                        + f_flaeche * squeeze(zielfunktion(:,p,19))...
                        + f_dammvol * squeeze(zielfunktion(:,p,20));
                else
                    optimierungsfunktion(:,1) = f_sa*((squeeze(zielfunktion(:,p,4))))...
                        + f_fuellung * (1 - squeeze(zielfunktion(:,p,17)))...
                        + f_anzahl * squeeze(zielfunktion(:,p,18))...
                        + f_flaeche * squeeze(zielfunktion(:,p,19))...
                        + f_dammvol * squeeze(zielfunktion(:,p,20));
                end
                
%                 max_fuellung = max(squeeze(zielfunktion(:,p,15)));
%                 min_fuellung = min(squeeze(zielfunktion(:,p,15)));
%                 optimierungsfunktion(:,1) = f_sa*(max_SA - (squeeze(zielfunktion(:,p,1))-min_SA))/(max_SA-min_SA)...
%                     + f_fuellung * (max_fuellung - (squeeze(zielfunktion(:,p,15))-min_fuellung))/(max_fuellung-min_fuellung)...
%                     + f_anzahl * squeeze(zielfunktion(:,p,16))...
%                     + f_flaeche * squeeze(zielfunktion(:,p,17))...
%                     + f_dammvol * squeeze(zielfunktion(:,p,18));
%                 optimierungsfunktion(:,1) = f_sa*(max_SA-(squeeze(zielfunktion(:,p,1)-min_SA)/(max_SA-min_SA)))...
%                     + f_fuellung * (1-squeeze(zielfunktion(:,p,15)))...
%                     + f_anzahl * squeeze(zielfunktion(:,p,16))...
%                     + f_flaeche * squeeze(zielfunktion(:,p,17))...
%                     + f_dammvol * squeeze(zielfunktion(:,p,18));
                
                komb = find(optimierungsfunktion==min(min(optimierungsfunktion)));
                
                best_becken_kombination = pegel_becken{1,p}(ismember(pegel_becken{1,p},Beckenkombination(i_no_of_basins).becken_real(:,komb)));
                
                beckenkomb_vorab_opt{p} = best_becken_kombination;
                if i_no_of_basins==size(Beckenkombination,1)
                    becken_kombination_final = Beckenkombination(i_no_of_basins).becken_real(:,komb);
                    beckenkombination_opt_all{p} = becken_kombination_final(becken_kombination_final~=0)';
                    drosselweite_opt_all{p} = Beckenkombination(i_no_of_basins).drosselweite(becken_kombination_final~=0,komb,E_opt)';
                end
            end
            if i_no_of_basins==size(Beckenkombination,1)
                beckenkombination_opt = beckenkombination_opt_all{TG};
                drosselweite_opt = drosselweite_opt_all{TG};
                if isempty(beckenkombination_opt)==0
                    optimierungsfunktion_opt = [optimierungsfunktion(komb),f_sa,f_fuellung,f_anzahl,f_flaeche,f_dammvol,...
                        1-squeeze(zielfunktion(komb,TG,1))/max(squeeze(zielfunktion(:,TG,1))),1-squeeze(zielfunktion(komb,TG,15)),...
                        squeeze(zielfunktion(komb,TG,16)),squeeze(zielfunktion(komb,TG,17)),squeeze(zielfunktion(komb,TG,18)),...
                        E_no,E_opt,E_opt_komb];
                else
                    optimierungsfunktion_opt = [NaN,f_sa,f_fuellung,f_anzahl,f_flaeche,f_dammvol,...
                        NaN,NaN,NaN,NaN,NaN,...
                        E_no,E_opt,E_opt_komb];
                end
                if isempty(komb)==1 
                    komb=1;
                end
            zielfunktion_opt = squeeze(zielfunktion(komb,:,:));
            lokal_opt = squeeze(lokal(komb,:,:));
            regional_opt = squeeze(regional(komb,:,:));
            daten_ges_opt = squeeze(daten_ges(komb,:,:));
            
            
            if isempty(beckenkombination_opt)==1 && length(Beckenkombination)==1
                if E_no==E_opt && E_no==E_opt_komb
                    Beckenkombination_final(n_opt_komb).becken_real = Beckenkombination(1).becken_real;
                end
                Beckenkombination_final(n_opt_komb).drosselweite(:,E_opt) = squeeze(Beckenkombination(1).drosselweite(:,runde,E_opt));
            else
                if E_no==E_opt && E_no==E_opt_komb
                    Beckenkombination_final(n_opt_komb).becken_real = beckenkombination_opt';
                end
                Beckenkombination_final(n_opt_komb).drosselweite(:,E_opt) = drosselweite_opt;
            end
            if E_no==E_opt && E_no==E_opt_komb
                Beckenkombination_final(n_opt_komb).opt_parameter = opt_parameter;
            end
            end
        end
        
        %% Save Results
        
        cd(path_results)
        if mOpt==1
            % Beckenkombination: basin numbers and throttle diameter
            if E_no == E_opt && E_no==E_opt_komb
                save(sprintf('Beckenkombination%d_Ekomb%d.mat',E_opt,abs(E_opt_komb)),'Beckenkombination')
            end
            
            % analysis files for basin combination efficiency
            save(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_info_sim',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'info_sim','-v7.3')
            save(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_info_beckenkombination',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'info_beckenkombination','-v7.3')
            save(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_zielfunktion',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'zielfunktion','-v7.3')
            save(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_lokal',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'lokal','-v7.3')
            save(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_regional',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'regional','-v7.3')
            save(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_daten_ges',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'daten_ges','-v7.3')
            
            cd(strcat(path_allg,'\Results\Sim_overview'))
            
            while exist('open_characteristics.txt','file')==2
                pause(5*rand(1))
            end
            save('open_characteristics.txt','Gebiet','-ascii')
            load('characteristics.mat')
            delete('open_characteristics.txt')
            zeile = characteristics(ismember(characteristics(:,2:6),[Gebiet,TG,mOpt,E_no,E_opt_komb],'rows'),1);
            
            while exist(sprintf('open_status_%d.txt',zeile),'file')==2
                pause(5*rand(1))
            end
            save(sprintf('open_status_%d.txt',zeile),'Gebiet','-ascii')
            load(sprintf('status_%d.mat',zeile))
            status = eval(sprintf('status_%d', zeile));
            status(2) = 1;    % status: läuft oder nicht
            status(3) = n_anzahl_sim/anzahl_sim*100; % Fotschritt 8%)
            status(5) = now;    % Zeit
            tmp = ['status_', num2str(zeile),'= status'];
            eval(tmp)
            clear tmp
            save(sprintf('status_%d.mat',zeile),sprintf('status_%d',zeile))
            delete(sprintf('open_status_%d.txt',zeile))
            
            %         cd(strcat(path_allg,'\Results'))
            %
            %         count = 0; err_count = 0;
            %         while count == err_count
            %             try     load('simulation_overview.mat','characteristics')
            %             catch;   err_count = err_count + 1;     end
            %             count = count + 1; pause(3)
            %         end
            %
            %         zeile = characteristics(ismember(characteristics(:,2:6),[Gebiet,TG,mOpt,E_no,E_opt_komb],'rows'),1);
            %         count = 0; err_count = 0;
            %         while count == err_count
            %             try     load('simulation_overview.mat',sprintf('status_%d',zeile))
            %             catch;   err_count = err_count + 1;     end
            %             count = count + 1; pause(3)
            %         end
            %
            % %         load('simulation_overview.mat','characteristics')
            % %         zeile = find(ismember(characteristics(:,2:6),[Gebiet,TG,mOpt,E_no,E_opt_komb],'rows')==1);
            % %         load('simulation_overview.mat',sprintf('status_%d',zeile))
            %         status = eval(sprintf('status_%d', zeile));
            %         status(2) = 1;    % status: läuft oder nicht
            %         status(3) = n_anzahl_sim/anzahl_sim*100; % Fotschritt 8%)
            %         status(5) = now;    % Zeit
            %         tmp = ['status_', num2str(zeile),'= status'];
            %         eval(tmp)
            %         clear tmp
            % %         assignin ('base',['status_' num2str(zeile)], status);
            %         save('simulation_overview','-append',sprintf('status_%d',zeile))
            
            %         m = matfile('simulation_overview.mat','Writable',true);
            %         zeile = find(ismember(m.simulation_overview(:,1:5),[Gebiet,TG,mOpt,E_opt,E_opt_komb],'rows')==1);
            %         m.simulation_overview(zeile,7) = n_anzahl_sim/anzahl_sim*100;
            %         m.simulation_overview(zeile,8) = 1;
            %         m.simulation_overview(zeile,10) = now;
            %         save('simulation_overview.mat','simulation_overview')
            
        elseif mOpt==2 && i_no_of_basins==size(Beckenkombination,1)
            if E_no == E_opt
                save(sprintf('Beckenkombination%d_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_opt,opt_parameter),'Beckenkombination')
            end
            save(sprintf('Eno%d_Eopt%d_kombE%d_%02dnoBasins_opt_beckenkombination_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,i_no_of_basins,opt_parameter),'beckenkombination_opt','optimierungsfunktion_opt','zielfunktion_opt','-v7.3')
            save(sprintf('Eno%d_Eopt%d_kombE%d_%02dnoBasins_opt_lokal_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,i_no_of_basins,opt_parameter),'lokal_opt','-v7.3')
            save(sprintf('Eno%d_Eopt%d_kombE%d_%02dnoBasins_opt_regional_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,i_no_of_basins,opt_parameter),'regional_opt','-v7.3')
            save(sprintf('Eno%d_Eopt%d_kombE%d_%02dnoBasins_opt_daten_ges_%4.2f_%4.2f_%4.2f_%4.2f_%4.2f.mat',E_no,E_opt,E_opt_komb,i_no_of_basins,opt_parameter),'daten_ges_opt','-v7.3')
            
            save(sprintf('Eoptkomb%d_beckenkombination_drosselweite_alleopt.mat',E_opt_komb),'Beckenkombination_final','-v7.3')
            
            disp('Data saved.')
            
        elseif mOpt==10
            % Beckenkombination: basin numbers and throttle diameter
            if E_no == E_opt && E_no==E_opt_komb
                save(sprintf('Beckenkombination%d_Ekomb%d.mat',E_opt,abs(E_opt_komb)),'Beckenkombination')
            end
            
            veraenderte_pegel = [];
            pegel_becken = becken_pos(:,1);
            for i = 1:length(pegel_becken)
            veraenderte_pegel = [veraenderte_pegel,pegel_becken(i),pegel(pegel_becken(i)).neighbors.alle_nachfolger];
            end
            veraenderte_pegel = unique(veraenderte_pegel);
            
            for p_no = 1:length(veraenderte_pegel)
                ganglinien_qgko(:,p_no) = pegel(veraenderte_pegel(p_no)).inflow_wasim(:,2);
                ganglinien_inflow(:,p_no) = pegel(veraenderte_pegel(p_no)).inflow{end}(:,2);
                ganglinien_outflow(:,p_no) = pegel(veraenderte_pegel(p_no)).outflow{end}(:,2);
            end
            
            for i_b_no = 1:size(info_beckenkombination,2)
                b_no = info_beckenkombination(1,i_b_no,1);
                ganglinien_storage(:,i_b_no) = pegel(b_no).becken.S(:,2);                
            end
            
            % analysis files for basin combination efficiency
            namen = dir(sprintf('Eno%d_Eopt%d_method%d*',E_no,E_opt,method));
            save(sprintf('Eno%d_Eopt%d_method%d_mOpt2_%4.2f_v2.mat',E_no,E_opt,method,opt_parameter(1)),...
                'info_sim','info_beckenkombination',...
                'zielfunktion','lokal','regional','daten_ges','opt_parameter',...
                'veraenderte_pegel','ganglinien_qgko','ganglinien_inflow',...
                'ganglinien_outflow','ganglinien_storage','-v7.3')
%             save(sprintf('Eno%d_Eopt%d_method%d_%d',E_no,E_opt,method,length(namen)+1),'info_sim','info_beckenkombination',...
%                 'zielfunktion','lokal','regional','daten_ges','opt_parameter',...
%                 'veraenderte_pegel','ganglinien_qgko','ganglinien_inflow','ganglinien_outflow','-v7.3')
%             save(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_info_beckenkombination',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'info_beckenkombination','-v7.3')
%             save(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_zielfunktion',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'zielfunktion','-v7.3')
%             save(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_lokal',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'lokal','-v7.3')
%             save(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_regional',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'regional','-v7.3')
%             save(sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_daten_ges',E_no,E_opt,abs(E_opt_komb),i_no_of_basins),'daten_ges','-v7.3')
            
        end
        % Disp der Rechenzeit
        zeit_ende = now;
        dauer_gesamt = (zeit_ende-zeit_start)*24*3600;
        fprintf('wirksamkeit_Eno%d_Eopt%d_%dnoBasins, Gesamtdauer: %8.2f s\n',E_no,E_opt,i_no_of_basins,dauer_gesamt)
        
    end
end

ende_global = now;
dauer_gesamt = (ende_global-start_global)*24*3600;
fprintf('Gesamtdauer: %8.2f s\n',dauer_gesamt)