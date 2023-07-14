%%% Import and save catchment data %%%

for Gebiet = [3]
    if Gebiet==1
        TGs = [1,3,9];
    elseif Gebiet==2
        TGs = [1,2,6];
    elseif Gebiet==3
        TGs = 2;
    elseif Gebiet==4
        TGs = 1;
    end
    for TG = TGs
        clearvars -except Gebiet TGs TG
        
        read_grids = 0;
        read_ereignisse = 0;    % Ereignisse, Q_mean, HHQ
        read_ereignisse_real = 1;
        read_richards = 0;
        read_routing = 0;
        
        %% Pfaddefinitionen
        path_APs = 'PFAD TOBI';
        
        pfad_all = 'PFAD ERGEBNIS' ;
        path_allg = strcat(path_APs,'\Inputdaten');
        
        pfad_kalib = 'PFAD KALIBRIERUNG' ;
        
        path_data = strcat(path_allg,'\',sprintf('G%d_TG%d_input',Gebiet,TG));
        mkdir(path_data);
        grid_length = 100;
        
        if Gebiet == 1
            A_ges = 1167.9;   % Einzugsgebietsgröße (immer Gesamtgebiet)
            grid_name = 'main100';
            
            pfad = strcat(pfad_all,'Main');
            %%% Abhängig vom Teilgebiet
            if TG == 1
                E_name = {'Ereignis_HQ_adv_G1_210.0_12','Ereignis_HQ_adv_G1_320.0_16',...
                    'Ereignis_HQ_adv_G1_460.0_17','Ereignis_HQ_adv_G1_581.0_20','Ereignis_HQ_adv_G1_732.0_21',...
                    'Ereignis_HQ_konv_G1_210.0_29','Ereignis_HQ_konv_G1_320.0_22'};
            elseif TG == 3
                E_name = {'Ereignis_HQ_adv_G3_85.0_15','Ereignis_HQ_adv_G3_130.0_15',...
                    'Ereignis_HQ_adv_G3_180.0_17','Ereignis_HQ_adv_G3_228.0_17','Ereignis_HQ_adv_G3_289.0_18',...
                    'Ereignis_HQ_konv_G3_85.0_4','Ereignis_HQ_konv_G3_130.0_6'};
            elseif TG == 9
                E_name = {'Ereignis_HQ_adv_G9_70.0_14','Ereignis_HQ_adv_G9_90.0_15',...
                    'Ereignis_HQ_adv_G9_110.0_16','Ereignis_HQ_adv_G9_134.0_7','Ereignis_HQ_adv_G9_164.0_5',...
                    'Ereignis_HQ_konv_G9_70.0_8','Ereignis_HQ_konv_G9_90.0_10'};
            else
                error('Falsche Kombination aus Gebiet und Teilgebiet gewählt.')
            end
            
            path_tanalys_allgages = strcat(path_APs,'PFAD TANALYS' );          % tanalys output for the simulation (with all additional gages)
            path_tanalys_gage1 = strcat(path_APs,'PFAD TANALYS 1 PEGEL');  % tanalys output for analysis (only with the main gage)
            
        elseif Gebiet == 2
        end
        
        
        %% Grid-Analyse (Flussnetz)
        
        if read_grids == 1
            
            %%% loading input grids (output ascii-files from tanalys)
            
            % copy grids to working directory
            cd(path_tanalys_gage1)
            copyfile(sprintf('%s.dgm',grid_name),path_data);
            copyfile(sprintf('%s.dep',grid_name),path_data);
            copyfile(sprintf('%s.ezg',grid_name),path_data);
            copyfile(sprintf('%s.fzs',grid_name),path_data);
            copyfile(sprintf('%s.flk',grid_name),path_data);
            copyfile('gridasci.exe',path_data);
            cd(path_tanalys_allgages)
            copyfile(sprintf('%s.pur',grid_name),path_data);
            
            % conversion of grids into ascii-files
            cd(path_data)
            [~,~] = system(sprintf('gridasci.exe %s.dgm dgm.txt',grid_name));
            [~,~] = system(sprintf('gridasci.exe %s.dep dep.txt',grid_name));
            [~,~] = system(sprintf('gridasci.exe %s.ezg ezg.txt',grid_name));
            [~,~] = system(sprintf('gridasci.exe %s.fzs fzs.txt',grid_name));
            [~,~] = system(sprintf('gridasci.exe %s.flk flk.txt',grid_name));
            [~,~] = system(sprintf('gridasci.exe %s.pur pur.txt',grid_name));
            
            % loading input grids (output ascii-files from tanalys)
            cd(path_data)
            DGM = importdata('dgm.txt','	',6);
            str_temp = char(DGM.textdata(1,1));
            ncols = str2num(str_temp(6:length(str_temp)));
            clear str_temp
            str_temp = char(DGM.textdata(2,1));
            nrows = str2num(str_temp(6:length(str_temp)));
            clear str_temp
            DGM = DGM.data;
            DEP = importdata('dep.txt','	',6);
            DEP = DEP.data;
            EZG = importdata('ezg.txt','	',6);
            EZG = EZG.data;
            FZS = importdata('fzs.txt','	',6);
            FZS = FZS.data;
            FLK = importdata('flk.txt','	',6);
            FLK = FLK.data;
            PUR = importdata('pur.txt','	',6);
            PUR = PUR.data;
            
            
            %%% analysis of the river net
            RIV = DEP;
            RIV(EZG<0) = -1;        % delete rivers outside the catchment
            RIV(RIV==-9999) = 0;    % set nodata-values to 0
            RIV(RIV>0) = 1;         % set all rivers to 1
            DEP(EZG<0) = 0;
            DEP(DEP==-9999) = 0;
            PUR(PUR==-9999) = 0;    % set nodata-values to 0
            
            % struct with all necessary grids for characterizing a flow point
            grids = struct('DGM',DGM,'RIV',RIV,'EZG',EZG,'FZS',FZS,'FLK',FLK,'PUR',PUR);
            
            %%% pool flow points in river sections
            X_CORD = repmat(1:ncols,nrows,1);       % matric with x-coordinates
            Y_CORD = repmat((1:nrows)',1,ncols);    % matrix with y-coordinates
            % table with river points and coordinates
            river_points(:,2:4) = [reshape(X_CORD,[],1) reshape(Y_CORD,[],1) reshape(RIV,[],1)];
            river_points = river_points(river_points(:,4)==1,1:3);
            river_points(:,1) = 1:length(river_points);
            % define flow points and their upstream and downstream neighbors
            for i = 1:length(river_points)
                flusspunkte(i) = Flusspunkt(i,grids,river_points);
                flusspunkte(i).fun_nachfolger
            end
            for i = 1:length(flusspunkte)
                flusspunkte.fun_vorgaenger(i);
            end
            
            flussnetz = river_points;   % vector with all river points
            i = 0;
            % loop to define river sections
            while isempty(flussnetz)==0 % while not all flow points are aggregated to a river section
                i = i+1;
                point = flussnetz(1,1); % selection of any point out of flow points which haven't been aggregated yet
                number = i;             % number of the river section
                
                gages = [];
                % define flow point with [number x-coordinate y-coordinate elevation flow_time]
                flowpoints(1,:) = [point flusspunkte(point).coordinates(1,1) flusspunkte(point).coordinates(1,2) flusspunkte(point).elevation flusspunkte(point).fliesszeit];
                % check whether there is a gage on the flow point
                pur = PUR(flusspunkte(point).coordinates(1,2),flusspunkte(point).coordinates(1,1));
                if pur>0 % if there is a gage -> add the point to the vector gages
                    gages(end+1,:) = [pur flusspunkte(point).coordinates(1,1) flusspunkte(point).coordinates(1,2) flusspunkte(point).elevation flusspunkte(point).fliesszeit];
                end
                
                % Vorgaenger
                point = flusspunkte(flussnetz(1,1)).vorgaenger; % find upstream neighbor
                while length(point)==1 % while there is one (no additional inflow, not uppermost point of the stream)
                    % define flow point and check whether there is a gage on the flow point
                    flowpoints(end+1,:) = [point flusspunkte(point).coordinates(1,1) flusspunkte(point).coordinates(1,2) flusspunkte(point).elevation flusspunkte(point).fliesszeit];
                    pur = PUR(flusspunkte(point).coordinates(1,2),flusspunkte(point).coordinates(1,1));
                    if pur>0
                        gages(end+1,:) = [pur flusspunkte(point).coordinates(1,1) flusspunkte(point).coordinates(1,2) flusspunkte(point).elevation flusspunkte(point).fliesszeit];
                    end
                    point = flusspunkte(point).vorgaenger; % find new upstream neighbor
                end
                
                % Nachfolger
                point = flusspunkte(flussnetz(1,1)).nachfolger; % find downstream neighbor
                if isempty(point)==0    % if a downstream neighbor exists (point is not the catchment outlet)
                    while length(flusspunkte(point).vorgaenger)==1 % while the upstream neighbor has only one neighbor (last point of a river section is the one after an inflow)
                        % define flow point and check whether there is a gage on the flow point
                        flowpoints(end+1,:) = [point flusspunkte(point).coordinates(1,1) flusspunkte(point).coordinates(1,2) flusspunkte(point).elevation flusspunkte(point).fliesszeit];
                        pur = PUR(flusspunkte(point).coordinates(1,2),flusspunkte(point).coordinates(1,1));
                        if pur>0
                            gages(end+1,:) = [pur flusspunkte(point).coordinates(1,1) flusspunkte(point).coordinates(1,2) flusspunkte(point).elevation flusspunkte(point).fliesszeit];
                        end
                        point = flusspunkte(point).nachfolger; % find new downstream neighbor
                        if isempty(point)==1
                            break
                        end
                    end
                end
                
                % if river section has no gages
                if isempty(gages)
                    gages = NaN;
                end
                
                flowpoints = sortrows(flowpoints,5); % sort flow points after flow time
                flussnetz = setdiff(flussnetz(:,1),flowpoints(:,1));    % delete points from vector flowpoints in vector flussnetz (new vector flussnetz)
                
                % adapt table "river_points"
                river_points(ismember(river_points(:,1),flowpoints(:,1))==1,4) = i; % col 4: number of the river section
                if isnan(gages)==0  % if gages exist in the river section i
                    [Lia,Locb] = ismember(river_points(:,2:3),gages(:,2:3),'rows');
                    river_points(Lia==1,5) = gages(Locb(Lia~=0),1);                     % col 5: gage number (if gage exists)
                    river_points(river_points(:,5)==0,5)=NaN;                           % col 5: "NaN" (if no gage exists)
                end
                
                % define struct "river" for the river section i
                river(i) = struct('number',number,'flowpoints',flowpoints,'gages',gages);
                clear flowpoints gage
            end
            
            % grids + Gewässerstruktur speichern
            cd(path_data)
            save('river.mat','river','grids','river_points')
            fprintf('river.mat in %s gespeichert.\n',path_data)
        end
        
        
        
        %% Load Discharges from WaSiM Output
        
        if read_ereignisse == 1
            
            for E_no = 1:length(E_name)
                
                path_wasim = strcat(pfad,'\',E_name{E_no});
                path_wasim_output = strcat(path_wasim,'\Output');
                
                clear qgko qgko_date qgko_abs qgko_spec
                clear qges qges_date qges_abs qges_spec
                clear gwin gwin_date gwin_abs gwin_spec
                clear qinf qinf_date qinf_abs qinf_spec
                clear datum A_ezg_mat A_teilezg_mat
                
                cd(path_wasim_output)
                %%% qgko (routed discharge of the subcatchment)
                qgko = importdata(sprintf('qgko%s.stat',grid_name));
                
                % convert date (string -> vector)
                datum(:,1) = str2double(qgko.textdata(4:size(qgko.textdata,1)-89,1));
                datum(:,2) = str2double(qgko.textdata(4:size(qgko.textdata,1)-89,2));
                datum(:,3) = str2double(qgko.textdata(4:size(qgko.textdata,1)-89,3));
                datum(:,4) = str2double(qgko.textdata(4:size(qgko.textdata,1)-89,4));
                % convert date (vector -> numeric value)
                qgko_date(:,1) = datenum(datum(:,1),datum(:,2),datum(:,3),datum(:,4),0,0);
                
                % areas of all subcatchments
                A_ezg(1,:) = qgko.data(1,:).*A_ges;
                
                % routed discharge
                qgko_spec(:,:) = qgko.data(2:size(qgko.data,1)-88,:);   % specific discharge [l/s/km^2]
                A_ezg_mat(:,:) = repmat(A_ezg(1,:),size(qgko_spec,1),1);
                qgko_abs = qgko_spec.*A_ezg_mat/3.6;          % absolute discharge [m^3/s]
                
                %%% qges (generated runoff of the subcatchment area)
                qges = importdata(sprintf('qges%s.stat',grid_name));
                
                % convert date (string -> vector)
                datum(:,1) = str2double(qges.textdata(4:size(qges.textdata,1),1));
                datum(:,2) = str2double(qges.textdata(4:size(qges.textdata,1),2));
                datum(:,3) = str2double(qges.textdata(4:size(qges.textdata,1),3));
                datum(:,4) = str2double(qges.textdata(4:size(qges.textdata,1),4));
                % convert date (vector -> numeric value)
                qges_date(:,1) = datenum(datum(:,1),datum(:,2),datum(:,3),datum(:,4),0,0);
                
                % areas of all subcatchments (not including the areas of upsteam subcatchments)
                A_teilezg(1,:) = qges.data(1,1:end-1).*A_ges;
                
                
                % routed discharge
                qges_spec(:,:) = qges.data(2:size(qges.data,1),1:end-1); % specific discharge [l/s/km^2]
                A_teilezg_mat(:,:) = repmat(A_teilezg(1,:),size(qges_spec,1),1);
                qges_abs = qges_spec.*A_teilezg_mat/3.6;       % absolute discharge [m^3/s]
                
                %%% qges (generated runoff of the subcatchment area)
                gwin = importdata(sprintf('gwin%s.stat',grid_name));
                
                % convert date (string -> vector)
                datum(:,1) = str2double(gwin.textdata(4:size(gwin.textdata,1),1));
                datum(:,2) = str2double(gwin.textdata(4:size(gwin.textdata,1),2));
                datum(:,3) = str2double(gwin.textdata(4:size(gwin.textdata,1),3));
                datum(:,4) = str2double(gwin.textdata(4:size(gwin.textdata,1),4));
                % convert date (vector -> numeric value)
                gwin_date(:,1) = datenum(datum(:,1),datum(:,2),datum(:,3),datum(:,4),0,0);
                
                % areas of all subcatchments (not including the areas of upsteam subcatchments)
                A_teilezg(1,:) = gwin.data(1,1:end-1).*A_ges;
                
                
                % routed discharge
                gwin_spec(:,:) = gwin.data(2:size(gwin.data,1),1:end-1); % specific discharge [l/s/km^2]
                A_teilezg_mat(:,:) = repmat(A_teilezg(1,:),size(gwin_spec,1),1);
                gwin_abs = gwin_spec.*A_teilezg_mat/3.6;       % absolute discharge [m^3/s]
                
                
                %%% Speichern der Daten in Struct
                Ereignis(E_no).qgko_date = qgko_date;
                Ereignis(E_no).qgko_abs = qgko_abs;
                Ereignis(E_no).qgko_spec = qgko_spec;
                
                Ereignis(E_no).qges_date = qges_date;
                Ereignis(E_no).qges_abs = qges_abs;
                Ereignis(E_no).qges_spec = qges_spec;
                
                Ereignis(E_no).gwin_date = gwin_date;
                Ereignis(E_no).gwin_abs = gwin_abs;
                Ereignis(E_no).gwin_spec = gwin_spec;
                
                % Bestehende Becken berücksichtigen
                if Gebiet==2
                    path_wasim_SQ = 'PFAD BESTEHENDE BECKEN/SEEN';
                    existing_structures(E_no,9)=fun_load_basin('sp1_teg.txt','sp1_teg.ddp,',path_wasim_SQ,path_wasim_output);
                    existing_structures(E_no,10)=fun_load_basin('sp2_sch.txt','sp2_sch.ddp,',path_wasim_SQ,path_wasim_output);
                    existing_structures(E_no,11)=fun_load_basin('sp3_hrb.txt','sp3_hrb.ddp,',path_wasim_SQ,path_wasim_output);
                    existing_structures(E_no,12)=fun_load_basin('sp4_hrb.txt','sp4_hrb.ddp,',path_wasim_SQ,path_wasim_output);
                else
                    existing_structures = [];
                end
                
            end
            
            %%% MQ einlesen
            pfad_stat_values = strcat(pfad_all,'\analysis');
            cd(pfad_stat_values)
            
            load(sprintf('G%d_TG%d_stat_values.mat',Gebiet,TG));
            Q_mean = stat_values_aneich(1,2:end);
            
            %%% HHQ der Ereignisse berechnen
            for i = 1:7
                HHQ_Ereignisse_tmp(i,:) = max(Ereignis(i).qgko_abs,[],1);
            end
            HHQ_Ereignisse = max(HHQ_Ereignisse_tmp,[],1);
            
            
            % Ereignisse speichern
            cd(path_data)
            save('ereignis.mat','A_teilezg','A_ezg','Ereignis','existing_structures','Q_mean','HHQ_Ereignisse')
            fprintf('ereignis.mat in %s gespeichert.\n',path_data)

        end
        
        %% reale Ereignisse aus der Kalibrierungszeitreihe
        
        if read_ereignisse_real == 1
            
            cd(pfad_kalib)
            load(sprintf('G%d_output_TOBI.mat',Gebiet),'qgko_TG_m3s','qdir_ZG_mm','qifl_ZG_mm','qbas_ZG_mm','date')
            
            cd(path_data)
            load('ereignis.mat','A_ezg','A_teilezg')
            
            cd('PFAD KALIBRIERUNG')
                
            ereignis_zeitraeume
            
            A_teilezg_mat(:,:) = repmat(A_teilezg(1,:),size(qdir_ZG_mm,1),1);
            qgko_abs = qgko_TG_m3s;
            qdir_abs = qdir_ZG_mm.*A_teilezg_mat/3.6;       % absolute discharge [m^3/s]
            qifl_abs = qifl_ZG_mm.*A_teilezg_mat/3.6;       % absolute discharge [m^3/s]
            qbas_abs = qbas_ZG_mm.*A_teilezg_mat/3.6;       % absolute discharge [m^3/s]
                
            for E_no = 1:length(ereignisse_plot)
                
                zeile1 = find(date == datenum(mat_beginn(ereignisse_plot(E_no),:)));
                zeile2 = find(date == datenum(mat_ende(ereignisse_plot(E_no),:)));
                zeilen = zeile1:zeile2;
                
                
                %%% Speichern der Daten in Struct
                Ereignis(E_no).qgko_date = date(zeilen,1);
                Ereignis(E_no).qgko_abs = qgko_abs(zeilen,:);
                
                Ereignis(E_no).qdir_date = date(zeilen,1);
                Ereignis(E_no).qdir_abs = qdir_abs(zeilen,:);
                
                Ereignis(E_no).qifl_date = date(zeilen,1);
                Ereignis(E_no).qifl_abs = qifl_abs(zeilen,:);
                
                Ereignis(E_no).qbas_date = date(zeilen,1);
                Ereignis(E_no).qbas_abs = qbas_abs(zeilen,:);
                
                Ereignis(E_no).gwin_date = date(zeilen,1);
                Ereignis(E_no).gwin_abs = zeros(size(Ereignis(E_no).qbas_abs));
                
                                
                % Bestehende Becken berücksichtigen
                if Gebiet==2
                    path_wasim_SQ = 'PFAD BESTEHENDE BECKEN/SEEN';
                    existing_structures(E_no,9)=fun_load_basin('sp1_teg.txt','sp1_teg.ddp,',path_wasim_SQ,path_wasim_output);
                    existing_structures(E_no,10)=fun_load_basin('sp2_sch.txt','sp2_sch.ddp,',path_wasim_SQ,path_wasim_output);
                    existing_structures(E_no,11)=fun_load_basin('sp3_hrb.txt','sp3_hrb.ddp,',path_wasim_SQ,path_wasim_output);
                    existing_structures(E_no,12)=fun_load_basin('sp4_hrb.txt','sp4_hrb.ddp,',path_wasim_SQ,path_wasim_output);
                else
                    existing_structures = [];
                end
                
            end
            
            %%% MQ einlesen
            pfad_stat_values = strcat(pfad_all,'\analysis');
            cd(pfad_stat_values)
            
            load(sprintf('G%d_TG%d_stat_values.mat',Gebiet,TG));
            Q_mean = stat_values_aneich(1,2:end);
            
            %%% HHQ der Ereignisse berechnen
            for i = 1:length(ereignisse_plot)
                HHQ_Ereignisse_tmp(i,:) = max(Ereignis(i).qgko_abs,[],1);
            end
            HHQ_Ereignisse = max(HHQ_Ereignisse_tmp,[],1);
            
            
            % Ereignisse speichern
            cd(path_data)
            save('ereignis_real.mat','A_teilezg','A_ezg','Ereignis','existing_structures','Q_mean','HHQ_Ereignisse')
            fprintf('ereignis_real.mat in %s gespeichert.\n',path_data)

        end
        
        
        %% read Richards-Datei
        if read_richards == 1
            
            for E_no = 1:length(E_name)
                
                path_wasim = strcat(pfad,'\',E_name{E_no});
                path_wasim_regress = strcat(path_wasim,'\Regress');
                
                % Datei einlesen
                fid=fopen(strcat(path_wasim_regress,'\','storage_richards.ftz_0'));
                txt =textscan(fid,'%s','delimiter','\n');
                rows = txt{1,1};
                fclose(fid);
                
                % Zeilen
                zeilen_begin_fluss = find(contains(rows,'(tracer/flow index and tributary subbasin)'))+1;
                zeilen_ende_fluss = find(contains(rows,'qvzuq[subchannel 1]'))-2;
                zeilen_begin_qvzuq = find(contains(rows,'qvzuq[subchannel 1]'))-1;
                zeilen_ende_qvzuq = find(contains(rows,'flow_mm[1]'))-1;
                
                % fluss-array speichern
                for i = 1:length(zeilen_begin_fluss)
                    leerzeichen = find(rows{zeilen_begin_fluss(i)-1}==' ');
                    nummer = str2num(rows{zeilen_begin_fluss(i)-1}(leerzeichen(1)+1:leerzeichen(2)-1));
                    richards(E_no,nummer).no = nummer;
                    
                    % Fluss
                    fluss_cell = rows(zeilen_begin_fluss(i):zeilen_ende_fluss(i));
                    data_tmp = sprintf('%s', fluss_cell{:});
                    leerzeichen = find(data_tmp==' ');
                    data_tmp2(1) = str2num(data_tmp(1:leerzeichen(1)-1));
                    for j = 1:length(leerzeichen)-1
                        data_tmp2(j+1) = str2num(data_tmp(leerzeichen(j)+1:leerzeichen(j+1)-1));
                    end
                    richards(E_no,nummer).fluss = data_tmp2';
                    clear data_tmp2 data_tmp leerzeichen
                    
                    qvzuq_cell = rows(zeilen_begin_qvzuq(i):zeilen_ende_qvzuq(i));
                    for j = 1:2:length(qvzuq_cell)
                        leerzeichen = find(qvzuq_cell{j}==' ');
                        
                        data_tmp(1) = str2num(qvzuq_cell{j}(1:leerzeichen(1)-1));
                        for k = 1:length(leerzeichen)-5
                            data_tmp(k+1) = str2num(qvzuq_cell{j}(leerzeichen(k)+1:leerzeichen(k+1)-1));
                        end
                        richards(E_no,nummer).Qaus((j+1)/2,:) = data_tmp;
                        
                    end
                    
                end
            end
            
            % Routing-Parameter speichern
            cd(path_data)
            save('richards.mat','richards')
            fprintf('richards.mat in %s gespeichert.\n',path_data)
        end
        
        
        %% read Routing-Struktur
        
        if read_routing == 1
            
            %%% Controlfile einlesen und Routing-Struktur selektieren
            pfad_control = strcat(pfad,'\',E_name{1},'\Control');
            
            fid=fopen(strcat(pfad_control,'\controlfile.txt'));
            txt =textscan(fid,'%s','delimiter','\n');
            rows = txt{1,1};
            fclose(fid);
            
            start_routing = find(contains(rows,'33                 #  timeoffset')==1)+1;
            ende_routing = find(contains(rows,'#   _                 _                     ___               _ _')==1)-3;
            routing_struktur = rows(start_routing:ende_routing);
            
            %%% Controlfile analysieren
            for i = 1:length(routing_struktur)
                if contains(routing_struktur{i},' OL ')==1
                    text = routing_struktur{i};
                    index = regexp(text,'OL')+3;
                    index2 = regexp(text(index:end),'(')-3;
                    nummer = str2double(text(index:index+index2));
                    
                    i1 = regexp(text,'kh')+3;
                    i2 = i1 + regexp(text(i1:end),',')-1;
                    i2 = i2(1);
                    routing(nummer).kh = str2double(text(i1:i2));
                    
                    i1 = regexp(text,'kv')+3;
                    i2 = i1 + regexp(text(i1:end),',')-1;
                    i2 = i2(1);
                    routing(nummer).kv = str2double(text(i1:i2));
                    
                    i1 = regexp(text,'Bh')+3;
                    i2 = i1 + regexp(text(i1:end),',')-1;
                    i2 = i2(1);
                    routing(nummer).Bh = str2double(text(i1:i2));
                    
                    i1 = regexp(text,'Bv')+3;
                    i2 = i1 + regexp(text(i1:end),',')-1;
                    i2 = i2(1);
                    routing(nummer).Bv = str2double(text(i1:i2));
                    
                    i1 = regexp(text,'Th')+3;
                    i2 = i1 + regexp(text(i1:end),',')-1;
                    i2 = i2(1);
                    routing(nummer).Th = str2double(text(i1:i2));
                    
                    i1 = regexp(text,'Mh')+3;
                    i2 = i1 + regexp(text(i1:end),',')-1;
                    i2 = i2(1);
                    routing(nummer).Mh = str2double(text(i1:i2));
                    
                    i1 = regexp(text,'Mv')+3;
                    i2 = i1 + regexp(text(i1:end),',')-1;
                    i2 = i2(1);
                    routing(nummer).Mv = str2double(text(i1:i2));
                    
                    i1 = regexp(text,'I')+3;
                    i2 = i1 + regexp(text(i1:end),',')-1;
                    i2 = i2(1);
                    routing(nummer).I = str2double(text(i1:i2));
                    
                    i1 = regexp(text,'L=')+2;
                    i2 = i1 + regexp(text(i1:end),',')-1;
                    i2 = i2(1);
                    routing(nummer).L = str2double(text(i1:i2));
                    
                    i1 = regexp(text,'AE')+3;
                    i2 = i1 + regexp(text(i1:end),')')-2;
                    i2 = i2(1);
                    routing(nummer).AE = str2double(text(i1:i2));
                end
            end
            
            for nummer = 2:length(routing)
            %%% Routing-Parameter berechnen
            dt = 1; % time step in h
            
            % Wertebereich der wq-Tabelle
            routing(nummer).qmax = 2000;
            routing(nummer).qmin = 0.1;
            intzahl = 90;
            
            dq = log(routing(nummer).qmax/routing(nummer).qmin)/intzahl;
            lnq = log(routing(nummer).qmax)+dq;
            
            q_in_tmp = lnq-(intzahl+1)*dq:dq:lnq-dq;
            routing(nummer).q_in = exp(q_in_tmp);
            routing(nummer).q_mm = .0036*routing(nummer).q_in;
            Q_in = routing(nummer).q_in*routing(nummer).AE/1000;
            
            
            for i = length(Q_in):-1:1
                T_v(i) = 0;
                v_v(i) = 0;
                Q_v(i) = 0;
                R_v(i) = 0; % Vorlang anfangs nicht überflutet
                
                Q_h(i) = Q_in(i);
                
                % Ermitteln von v, T, R für Hauptbett
                if i==length(Q_in)
                    v_h(i) = 1;
                else
                    v_h(i) = v_h(i+1);
                end
                
                [v_h(i),T_h(i),R_h(i)] = fun_IterateManning(Q_h(i),routing(nummer).Bh,v_h(i),...
                    routing(nummer).Mh,routing(nummer).I,50,.005);
                
                
                if T_h(i)>routing(nummer).Th
                    v_v(i) = v_h(i);
                    k=0;
                    
                    
                    Th_alt(i) = T_h(i);
                    Qh_alt(i) = Q_h(i);
                    
                    T_h(i) = (T_h(i) + T_v(i) + routing(nummer).Th)/2;
                    T_v(i) = T_h(i) - routing(nummer).Th;
                    
                    if Th_alt(i) > T_h(i) % Vorland bekommt mehr Wasser
                        Qdiffh = (Th_alt(i)-T_h(i))*routing(nummer).Bh*v_h(i);
                        Qdiffv = (Th_alt(i)-T_h(i))*routing(nummer).Bv*v_v(i);
                        
                        if Qdiffh<Qdiffv; Q_h(i)=Q_h(i)-Qdiffh; else; Q_h(i)=Q_h(i)-Qdiffv; end
                        Q_v(i) = Q_v(i)+(Qh_alt(i)-Q_h(i));
                    else % Hauptbett bekommt mehr Wasser
                        Qdiffh = (T_h(i)-Th_alt(i))*routing(nummer).Bh*v_h(i);
                        Qdiffv = (T_h(i)-Th_alt(i))*routing(nummer).Bv*v_v(i);
                        
                        if Qdiffh<Qdiffv; Q_h(i)=Q_h(i)+Qdiffh; else; Q_h(i)=Q_h(i)+Qdiffv; end
                        Q_v(i) = Q_v(i)-(Q_h(i)-Qh_alt(i));
                    end
                    
                    [v_h(i),T_h(i),R_h(i)] = fun_IterateManning(Q_h(i),routing(nummer).Bh,v_h(i),...
                        routing(nummer).Mh,routing(nummer).I,50,.005);
                    [v_v(i),T_v(i),R_v(i)] = fun_IterateManning(Q_v(i),routing(nummer).Bv,v_v(i),...
                        routing(nummer).Mv,routing(nummer).I,50,.005);
                    
                    while abs(T_h(i)-T_v(i)-routing(nummer).Th)>0.01 && k<49
                        
                        Th_alt(i) = T_h(i);
                        Qh_alt(i) = Q_h(i);
                        
                        T_h(i) = (T_h(i) + T_v(i) + routing(nummer).Th)/2;
                        T_v(i) = T_h(i) - routing(nummer).Th;
                        
                        if Th_alt(i) > T_h(i) % Vorland bekommt mehr Wasser
                            Qdiffh = (Th_alt(i)-T_h(i))*routing(nummer).Bh*v_h(i);
                            Qdiffv = (Th_alt(i)-T_h(i))*routing(nummer).Bv*v_v(i);
                            
                            if Qdiffh<Qdiffv; Q_h(i)=Q_h(i)-Qdiffh; else; Q_h(i)=Q_h(i)-Qdiffv; end
                            Q_v(i) = Q_v(i)+(Qh_alt(i)-Q_h(i));
                        else % Hauptbett bekommt mehr Wasser
                            Qdiffh = (T_h(i)-Th_alt(i))*routing(nummer).Bh*v_h(i);
                            Qdiffv = (T_h(i)-Th_alt(i))*routing(nummer).Bv*v_v(i);
                            
                            if Qdiffh<Qdiffv; Q_h(i)=Q_h(i)+Qdiffh; else; Q_h(i)=Q_h(i)+Qdiffv; end
                            Q_v(i) = Q_v(i)-(Q_h(i)-Qh_alt(i));
                        end
                        
                        [v_h(i),T_h(i),R_h(i)] = fun_IterateManning(Q_h(i),routing(nummer).Bh,v_h(i),...
                            routing(nummer).Mh,routing(nummer).I,50,.005);
                        [v_v(i),T_v(i),R_v(i)] = fun_IterateManning(Q_v(i),routing(nummer).Bv,v_v(i),...
                            routing(nummer).Mv,routing(nummer).I,50,.005);
                        
                        k = k+1;
                    end
                    
                end
                
                v_l(i) = v_v(i)*Q_v(i)/(Q_v(i)+Q_h(i)) + v_h(i)*Q_h(i)/(Q_v(i)+Q_h(i));
                t_Trans(i) = routing(nummer).L/v_l(i)/3600;      % in [h]
                
            end
            
            routing(nummer).vl = real(v_l);
            routing(nummer).vv = real(v_v);
            routing(nummer).vh = real(v_h);
            routing(nummer).dh = real(T_h);
            routing(nummer).dv = real(T_v);
            routing(nummer).t_Trans = real(t_Trans);
            routing(nummer).Q_h = real(Q_h);
            routing(nummer).Q_v = real(Q_v);
            routing(nummer).Q_in = real(Q_in);
            routing(nummer).QvzuQ = real(Q_v./(Q_v+Q_h));
            routing(nummer).n = ceil(routing(nummer).t_Trans(1)-routing(nummer).t_Trans(end));
            routing(nummer).delta = floor(routing(nummer).t_Trans./routing(nummer).n*60);
            end
            
            %%% Routing-Struct speichern
            cd(path_data)
            save('routing_struct.mat','routing')
            fprintf('routing_struct.mat in %s gespeichert.\n',path_data)
            
        end
        
    end
end