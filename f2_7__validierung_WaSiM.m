function [pegel] = f2_7__validierung_WaSiM(Gebiet,TG,E_no,becken_pos,beckenpunkte_select2,pegel,runde,i_no_of_basins,path_allg)

pfad_ereignisse = 'PFAD SIMULATIONSERGEBNISSE';

if Gebiet == 1
    A_ges = 1167.9;
    if TG ==1
        Ereignis_name = {'Ereignis_HQ_adv_G1_82.0_8','Ereignis_HQ_adv_G1_115.0_15',...
            'Ereignis_HQ_adv_G1_170.0_17','Ereignis_HQ_adv_G1_216.0_19','Ereignis_HQ_adv_G1_274.0_15',...
            'Ereignis_HQ_konv_G1_82.0_12','Ereignis_HQ_konv_G1_115.0_9'};      
    end
    path_wasim_basis = strcat(pfad_ereignisse,'Main');
    grid_name = 'main100';
elseif Gebiet == 2
    A_ges = 392.17; 
    if TG ==1
        Ereignis_name = {'Ereignis_HQ_adv_G1_82.0_8','Ereignis_HQ_adv_G1_115.0_15',...
            'Ereignis_HQ_adv_G1_170.0_17','Ereignis_HQ_adv_G1_216.0_19','Ereignis_HQ_adv_G1_274.0_15',...
            'Ereignis_HQ_konv_G1_82.0_12','Ereignis_HQ_konv_G1_115.0_9'};      
    elseif TG==2
        Ereignis_name = {'Ereignis_HQ_adv_G2_58.0_15','Ereignis_HQ_adv_G2_84.0_14',...
            'Ereignis_HQ_adv_G2_120.0_19','Ereignis_HQ_adv_G2_150.0_14','Ereignis_HQ_adv_G2_187.0_13',...
            'Ereignis_HQ_konv_G2_58.0_9','Ereignis_HQ_konv_G2_84.0_11'};      
    elseif TG==6
        Ereignis_name = {'Ereignis_HQ_adv_G6_30.0_8','Ereignis_HQ_adv_G6_42.0_16',...
            'Ereignis_HQ_adv_G6_55.0_17','Ereignis_HQ_adv_G6_67.0_13','Ereignis_HQ_adv_G6_83.0_14',...
            'Ereignis_HQ_konv_G6_30.0_9','Ereignis_HQ_konv_G6_42.0_8'};      
    end
    path_wasim_basis = strcat(pfad_ereignisse,'Mangfall');
    grid_name = 'm100';
elseif Gebiet == 3
    A_ges = 104;
    Ereignis_name = {'Ereignis_HQ_adv_G2_17.5_19','Ereignis_HQ_adv_G2_27.0_14',...
                'Ereignis_HQ_adv_G2_38.0_16','Ereignis_HQ_adv_G2_48.0_19','Ereignis_HQ_adv_G2_60.6_15',...
                'Ereignis_HQ_konv_G2_17.5_18','Ereignis_HQ_konv_G2_27.0_12'};
    path_wasim_basis = strcat(pfad_ereignisse,'Glonn');
    grid_name = 'glonn100';
elseif Gebiet == 4
    A_ges = 91.1;
    Ereignis_name = {'Ereignis_HQ_adv_G1_18.4_15','Ereignis_HQ_adv_G1_27.0_14',...
        'Ereignis_HQ_adv_G1_38.3_16','Ereignis_HQ_adv_G1_49.4_18','Ereignis_HQ_adv_G1_63.2_19',...
        'Ereignis_HQ_konv_G1_18.4_4','Ereignis_HQ_konv_G1_27.0_3'};
    path_wasim_basis = strcat(pfad_ereignisse,'Otterbach');
    grid_name = 'otter100';
end

path_wasim = strcat(path_wasim_basis,'\',Ereignis_name{E_no});
cd(strcat(path_allg,'\WaSiM'))
ordner = dir('Valid*');
pfad_valid = strcat(path_allg,'\WaSiM\',sprintf('Valid%d_%d_%d',length(ordner)+1,E_no,i_no_of_basins));
pfad_valid_control = strcat(pfad_valid,'\Control');
mkdir(pfad_valid_control);
if Gebiet == 1
    copyfile(strcat(path_wasim,'\Control\controlfile.txt'),...
    strcat(pfad_valid_control,'\controlfile.txt'));
else
    copyfile(strcat(path_wasim,'\Control\controlfile.txt'),...
    strcat(pfad_valid_control,'\controlfile.txt'));
end
copyfile(strcat(path_wasim,'\Control\vcomp120.dll'),...
    strcat(pfad_valid_control,'\vcomp120.dll'));
copyfile(strcat(path_wasim,'\Control\wasimvzo64.exe'),...
    strcat(pfad_valid_control,'\wasimvzo64.exe'));
copyfile(strcat(path_wasim,'\Regress\storage_richards.ftz_0'),...
    strcat(pfad_valid_control,'\storage_richards.ftz_0'));

pfad_conv = 'PFAD WASIM FILES';
addpath(pfad_conv)
path_work = strcat(path_allg,'\grids\',sprintf('G%d_TG%d',Gebiet,TG));
mkdir(path_work)

% Output-Ordner erstellen
pfad_valid_output = strcat(pfad_valid,'\Output');
mkdir(pfad_valid_output);

%% Becken-Grids erstellen
% Basisgrid einlesen (über dgm + ezg)
DGM = fun_GRIDimport(pfad_conv,strcat(path_wasim,'\Input'),path_work,grid_name,'dgm');
EZG = fun_GRIDimport(pfad_conv,strcat(path_wasim,'\Input'),path_work,grid_name,'ezg');
if Gebiet == 2
    LAK = fun_GRIDimport(pfad_conv,strcat(path_wasim,'\Input'),path_work,grid_name,'lak');
    MAXPOND = fun_GRIDimport(pfad_conv,strcat(path_wasim,'\Input'),path_work,grid_name,'maxpond');
else
    % Grids definieren
    LAK = EZG;  % LAK (Nummern der Becken auf Koordinaten, wenn im Gebiet --> ansonsten auf Pegelpunkt)
    LAK.data(isnan(LAK.data)==0) = NaN; % alles was kein Pond ist = NaN
    
    MAXPOND = EZG;  % MAXPOND (wie LAK, max Tiefe immer 20 m)
    MAXPOND.data(LAK.data~=0) = 0;  % alles was kein Pond ist = 0
    
    POND = EZG;
    POND.data(LAK.data~=0) = 0;  % alles was kein Pond ist = 0
end

if Gebiet == 2
    becken_exist = 4;
else
    becken_exist = 0;
end

for i = 1:length(becken_pos(:,1))
    no_becken = becken_pos(i,2);
    no_pegel = becken_pos(i,1);
    
    B_col = round((beckenpunkte_select2([beckenpunkte_select2.no_becken]==no_becken).coords(1) - DGM.xll)/100);
    B_row = round((beckenpunkte_select2([beckenpunkte_select2.no_becken]==no_becken).coords(2) - DGM.yll)/100);
    
    if isnan(LAK.data(B_row,B_col))==1
        if EZG.data(B_row,B_col) == no_pegel
            LAK.data(B_row,B_col) = i+becken_exist;
            MAXPOND.data(B_row,B_col) = 20;
        else
            if isnan(LAK.data(pegel(no_pegel).gage.coordinates(2),pegel(no_pegel).gage.coordinates(1)))==1 && ...
                    EZG.data(pegel(no_pegel).gage.coordinates(2),pegel(no_pegel).gage.coordinates(1))==no_pegel
                LAK.data(pegel(no_pegel).gage.coordinates(2),pegel(no_pegel).gage.coordinates(1)) = i+becken_exist;
                MAXPOND.data(pegel(no_pegel).gage.coordinates(2),pegel(no_pegel).gage.coordinates(1)) = 20;
            else
                [X,Y] = find(EZG.data(pegel(no_pegel).gage.coordinates(2)-1:pegel(no_pegel).gage.coordinates(2)+1,...
                pegel(no_pegel).gage.coordinates(1)-1:pegel(no_pegel).gage.coordinates(1)+1)==no_pegel & ...
                isnan(LAK.data(pegel(no_pegel).gage.coordinates(2)-1:pegel(no_pegel).gage.coordinates(2)+1,...
                pegel(no_pegel).gage.coordinates(1)-1:pegel(no_pegel).gage.coordinates(1)+1))==1,1);
                LAK.data(pegel(no_pegel).gage.coordinates(2)-2+X,pegel(no_pegel).gage.coordinates(1)-2+Y) = i+becken_exist;
                MAXPOND.data(pegel(no_pegel).gage.coordinates(2)-2+X,pegel(no_pegel).gage.coordinates(1)-2+Y) = 20;
            end
        end
    else
        if isnan(LAK.data(pegel(no_pegel).gage.coordinates(2),pegel(no_pegel).gage.coordinates(1)))==1 && ...
                EZG.data(pegel(no_pegel).gage.coordinates(2),pegel(no_pegel).gage.coordinates(1))==no_pegel
            LAK.data(pegel(no_pegel).gage.coordinates(2),pegel(no_pegel).gage.coordinates(1)) = i+becken_exist;
            MAXPOND.data(pegel(no_pegel).gage.coordinates(2),pegel(no_pegel).gage.coordinates(1)) = 20;
        else
            [X,Y] = find(EZG.data(pegel(no_pegel).gage.coordinates(2)-1:pegel(no_pegel).gage.coordinates(2)+1,...
                pegel(no_pegel).gage.coordinates(1)-1:pegel(no_pegel).gage.coordinates(1)+1)==no_pegel & ...
                isnan(LAK.data(pegel(no_pegel).gage.coordinates(2)-1:pegel(no_pegel).gage.coordinates(2)+1,...
                pegel(no_pegel).gage.coordinates(1)-1:pegel(no_pegel).gage.coordinates(1)+1))==1,1);
            LAK.data(pegel(no_pegel).gage.coordinates(2)-2+X,pegel(no_pegel).gage.coordinates(1)-2+Y) = i+becken_exist;
            MAXPOND.data(pegel(no_pegel).gage.coordinates(2)-2+X,pegel(no_pegel).gage.coordinates(1)-2+Y) = 20;
        end
        
    end
    
end

% Grids exportieren (in Control-Ordner)
fun_GRIDexport(pfad_conv,path_work,LAK,grid_name,'lak')
fun_GRIDexport(pfad_conv,path_work,MAXPOND,grid_name,'maxpond')
if Gebiet ~=2
    fun_GRIDexport(pfad_conv,path_work,POND,strcat('pond',grid_name),'grd')
end

copyfile(strcat(path_work,'\',grid_name,'.lak'),...
    strcat(pfad_valid_control,'\',grid_name,'.lak'));
copyfile(strcat(path_work,'\',grid_name,'.maxpond'),...
    strcat(pfad_valid_control,'\',grid_name,'.maxpond'));

	if Gebiet ~= 2
if exist(strcat(path_wasim,'\Regress\','pond',grid_name,'.grd'),'file')==0
    copyfile(strcat(path_work,'\','pond',grid_name,'.grd'),...
    strcat(path_wasim,'\Regress','\pond',grid_name,'.grd'));
end
end

%% Controlfile anpassen
fid=fopen(strcat(pfad_valid_control,'\controlfile.txt'));
txt =textscan(fid,'%s','delimiter','\n');
rows = txt{1,1};
fclose(fid);

rows_neu = rows;

% Input- und Output-Ordner anpassen
rows_neu{71} = sprintf('$set $mainpath	     = %s/',path_wasim);
rows_neu{80} = sprintf('$set $DefaultOutputDirectory = %s/',pfad_valid_output);

zeile_richards = find(contains(rows_neu,'$outpath//storage_richards.ftz    # if readgrids = 1')==1);
rows_neu{zeile_richards} = strcat(pfad_valid_control,'/storage_richards.ftz');

% Zeilen zum Einlesen der Grids anpassen (Gesamtanzahl anpassen!)
if Gebiet ==2
    zeile_standardGrids1 = find(contains(rows_neu,'$inpath//$lake_grid')==1);
    zeile_standardGrids2 = find(contains(rows_neu,'$inpath//$grid//.maxpond')==1,1);
    % rows_neu{zeile_standardGrids1+1} = '18';
    rows_neu{zeile_standardGrids1} = strcat(pfad_valid_control,'\','$grid//.lak 			lake_codes 				  0    # grid with a unique code for each lake');
    rows_neu{zeile_standardGrids2} = strcat(pfad_valid_control,'\','$grid//.maxpond 			max_ponding_storage        0    # grid with height of small dams around the fields for water ponding (in m). 0 if no ponding occurs');
else
    zeile_standardGrids1 = find(contains(rows_neu,'[standard_grids]')==1);
    zeile_standardGrids2 = find(contains(rows_neu,'#$inpath//')==1,2);
    rows_neu{zeile_standardGrids1+1} = '18';
    rows_neu{zeile_standardGrids2(1)} = strcat(pfad_valid_control,'\','$grid//.lak 			lake_codes 				  0    # grid with a unique code for each lake');
    rows_neu{zeile_standardGrids2(2)} = strcat(pfad_valid_control,'\','$grid//.maxpond 			max_ponding_storage        0    # grid with height of small dams around the fields for water ponding (in m). 0 if no ponding occurs');


% Zeilen bei Lake-Modell anpassen (ganzes Lake-Modell einfügen, weil
% auskommentiert)
zeile_lakemodel = find(contains(rows_neu,'[lake_model]')==1);
rows_neu{zeile_lakemodel+1} = '0';
rows_neu{zeile_lakemodel+2} = '1                            # method for recalculating DHM,';
rows_neu{zeile_lakemodel+3} = '# 1 = do not change the DHM, it refects already the ground surface of the lakes, ;';
rows_neu{zeile_lakemodel+4} = '# 2 = use mean_pond_grid to calculate dhm corrections';
rows_neu{zeile_lakemodel+5} = '# max_pond_grid will be used for mapping the cells pond content to a lake during model runs - so the lake level may well rise above the normal surface';
rows_neu{zeile_lakemodel+6} = '0.1  # Albedo_OpenWater (will be used only, when the pond is filled with water when calculating potential evaporation -> otherwise, the normal landuse for this cell is referenced for this parameter)';
rows_neu{zeile_lakemodel+7} = '#				10   # rsc for water (usage as above)';
rows_neu{zeile_lakemodel+8} = '0.4  # z0 for water (usage as above) ';
rows_neu{zeile_lakemodel+9} = '#				10.0 # LAI_OpenWater (usage as above)';
rows_neu{zeile_lakemodel+10} = '#				1.0  # VCF_OpenWater (usage as above)';
rows_neu{zeile_lakemodel+11} = '0 #$readgrids                   # readgrid code 0 do not read, 1 = read grids --> ';
rows_neu{zeile_lakemodel+12} = '# if 0, the initial valte for the POND-grid as Volume of Lakes and Reservoirs is set by V0 from the routing description, ';
rows_neu{zeile_lakemodel+13} = '# if readgrids=1, no initialization in done (POND-Grid is read in) but the Vakt-Value is set by the various grids';
end
% ponds in unsatzon-model aktivieren
zeile_unsatzon = find(contains(rows_neu,'[unsatzon_model]')==1);
% rows_neu{zeile_unsatzon(end)+5} = '1  # controlling surface storage in ponds:       0 = no ponds,       1 = using ponds for surface storage (pond depth as standard grid needed -> height of dams oround fields)';

for i = 1:length(becken_pos(:,1))
pegel_geo_temp(i,1) = beckenpunkte_select2([beckenpunkte_select2.no_becken]==becken_pos(i,2)).elev_dam;
pegel_geo_temp(i,2) = i+becken_exist;
pegel_geo_temp(i,3) = becken_pos(i,2);
pegel_geo_temp(i,4) = becken_pos(i,1);
end

pegel_geo_temp = sortrows(pegel_geo_temp,1);
beckennummer = pegel_geo_temp(:,3)';
pegelnummer = pegel_geo_temp(:,4)';
becken_kurznummer = pegel_geo_temp(:,2)';

% Becken in Routing-Struktur einfügen
for i = 1:length(pegelnummer)
    pegel_no = pegelnummer(i);
    becken_no = beckennummer(i);
    zeile_becken = find(contains(rows_neu,sprintf('TG %d ',pegel_no))==1);
    if isempty(zeile_becken)==0
        rows_neu(zeile_becken+3:length(rows_neu)+1) = rows_neu(zeile_becken+2:length(rows_neu));
        rows_neu{zeile_becken+2} = sprintf('and SP %d (file = sp%d_%d_%d_hrb.ddp , V0=0 , C0 =  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0, dTmin=60)',...
            becken_kurznummer(i),becken_kurznummer(i),becken_no,pegel_no);
    else
        error('Not possible')
        zeile_becken = find(contains(rows_neu,sprintf('OL %d ',pegel_no))==1);
        rows_neu{zeile_becken+2:length(rows_neu)+2} = rows_neu{zeile_becken:length(rows_neu)};   
        rows_neu{zeile_becken} = sprintf('TG %d (AE=%9.7f, AErel=1.0)',pegel_no,pegel(1,pegel_no).Aezg(pegel_no));
        rows_neu{zeile_becken+1} = sprintf('from SP %d (file = sp%d_%d_hrb.ddp , V0=0 , C0 =  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0, dTmin=60)',...
            becken_kurznummer(i),becken_kurznummer(i),becken_no,becken_no,pegel_no);
            
    end  
end

for i = 1:length(pegelnummer)
    pegel_no = pegelnummer(i);
    becken_no = beckennummer(i);
    % Becken-V-Q-Beziehung einfügen
    
    for k = 1:length(pegel(1,pegel_no).beckenparameter)
        zeile_becken = find(strcmp(rows_neu,'[multilayer_landuse] ')==1);
        zeile_becken = zeile_becken-22;
        rows_neu(zeile_becken+length(pegel(1,pegel_no).beckenparameter(k).ShQ)+3:length(rows_neu)+length(pegel(1,pegel_no).beckenparameter(k).ShQ)+3) ...
            = rows_neu(zeile_becken:length(rows_neu));
        rows_neu{zeile_becken+length(pegel(1,pegel_no).beckenparameter(k).ShQ)+2} = '';
        rows_neu{zeile_becken-1} = sprintf('[abstraction_rule_reservoir_%d]',becken_kurznummer(i));
        rows_neu{zeile_becken} = sprintf('%d',length(pegel(1,pegel_no).beckenparameter(k).ShQ));
        for j = 1:length(pegel(1,pegel_no).beckenparameter(k).ShQ)
            rows_neu{zeile_becken+j} = sprintf('%10.4f %10.4f',abs(pegel(1,pegel_no).beckenparameter(k).ShQ(j,[1,3])));
        end
        rows_neu{zeile_becken+j+1} = '';
    end
end

% for i = 1:length(Wirksamkeit_Beckenkombination(runde).pegelnummer)
%     pegel_no = Wirksamkeit_Beckenkombination(runde).pegelnummer(i);
%     becken_no = Wirksamkeit_Beckenkombination(runde).beckennummer(i);
%     % Becken-V-Q-Beziehung einfügen
%     zeile_becken = find(strcmp(rows_neu,'[multilayer_landuse] ')==1);
%     zeile_becken = zeile_becken-22;
%     diff_tmp = 0;
%     for k = 1:length(pegel(1,pegel_no).beckenparameter)
%         rows_neu(zeile_becken+length(pegel(1,pegel_no).beckenparameter(k).ShQ)-diff_tmp+3:length(rows_neu)+length(pegel(1,pegel_no).beckenparameter(k).ShQ)-diff_tmp+3) ...
%             = rows_neu(zeile_becken:length(rows_neu));
%         rows_neu{zeile_becken+length(pegel(1,pegel_no).beckenparameter(k).ShQ)-diff_tmp+2} = '';
%         rows_neu{zeile_becken-1} = sprintf('[abstraction_rule_reservoir_%d]',becken_no);
%         rows_neu{zeile_becken} = sprintf('%d',length(pegel(1,pegel_no).beckenparameter(k).ShQ)-diff_tmp);
%         for i = 1:length(pegel(1,pegel_no).beckenparameter(k).ShQ)-diff_tmp
% %             if i==1
% %                 rows_neu{zeile_becken+i} = sprintf('%10.4f %10.4f',0.001,pegel(1,pegel_no).beckenparameter(k).ShQ(i+diff_tmp+1,[3]));
% %             elseif i==2
% %                 rows_neu{zeile_becken+i} = sprintf('%10.4f %10.4f',1.000,pegel(1,pegel_no).beckenparameter(k).ShQ(i+diff_tmp,[3]));
% %             else
%                 rows_neu{zeile_becken+i-1} = sprintf('%10.4f %10.4f',pegel(1,pegel_no).beckenparameter(k).ShQ(i+diff_tmp,[1,3]));
% %             end
%         end
%         rows_neu{zeile_becken+i+1} = '';
%     end
% end

fid=fopen(strcat(pfad_valid_control,'\','controlfile_neu.txt'),'wt');
fprintf(fid,'%s\n',rows_neu{:});
fclose(fid);


%% Richards-Datei anpassen
clear rows

fid=fopen(strcat(pfad_valid_control,'\storage_richards.ftz_0'));
txt =textscan(fid,'%s','delimiter','\n');
rows = txt{1,1};
fclose(fid);

sum_tmp = ones(1,length(pegelnummer));

while sum(sum_tmp)>0
    k = find(sum_tmp>0,1);
    clear zeilen_becken
    zeilen_becken = find(contains(rows,sprintf('%d (routing subbasin with tributaries ',pegelnummer(k)))==1);
    if length(zeilen_becken)>1
        for i=1:length(zeilen_becken)
            zeile_str = char(rows{zeilen_becken(i)});
            if strcmp(zeile_str(1:length(num2str(pegelnummer(k)))),string(pegelnummer(k)))==1
                zeilen_becken_neu = zeilen_becken(i);
                break
            end
        end
    else 
        zeilen_becken_neu = zeilen_becken;
    end
    zeile_reservoir = find(contains(rows(zeilen_becken_neu:end),'0 (number of reservoirs)')==1,1)+zeilen_becken_neu-1;
    loc_becken = find([pegelnummer]==pegelnummer(k));
    anzahl_becken = length(loc_becken);
    rows_temp = rows;
    rows{zeile_reservoir} = sprintf('%d (number of reservoirs)',anzahl_becken);
    for n = 1:anzahl_becken
    rows{zeile_reservoir+n} = sprintf('0.0  (reservoir %d)',beckennummer(loc_becken(n)));
    end
    sum_tmp(loc_becken) = 0;

    rows(zeile_reservoir+n+1:length(rows_temp)+n) = rows_temp(zeile_reservoir+1:length(rows_temp));
    
end

fid=fopen(strcat(pfad_valid_control,'\storage_richards.ftz_0'),'wt');
fprintf(fid,'%s\n',rows{:});
fclose(fid);

%% Modell starten
cd(pfad_valid_control)
tic
[A,B] = system('wasimvzo64.exe controlfile_neu.txt','-echo');
toc

%% Ergebnisse importieren

    cd(pfad_valid_output)
    %%% qgko (routed discharge of the subcatchment)
    if exist(sprintf('qgko%s.stat',grid_name),'file')==2
        qgko = importdata(sprintf('qgko%s.stat',grid_name));
        qges = importdata(sprintf('qges%s.stat',grid_name));
    else
        Wirksamkeit_Beckenkombination.gof = [];
        return
    end

    clear datum qgko_date_becken
%         % convert date (string -> vector)
        datum(:,1) = str2double(qgko.textdata(4:size(qgko.textdata,1)-89,1));
        datum(:,2) = str2double(qgko.textdata(4:size(qgko.textdata,1)-89,2));
        datum(:,3) = str2double(qgko.textdata(4:size(qgko.textdata,1)-89,3));
        datum(:,4) = str2double(qgko.textdata(4:size(qgko.textdata,1)-89,4));
        % convert date (vector -> numeric value)
        qgko_date_becken(:,1) = datenum(datum(:,1),datum(:,2),datum(:,3),datum(:,4),0,0);

        % areas of all subcatchments
        A_ezg(1,:) = qgko.data(1,:).*A_ges;
        A_teilezg(1,:) = qges.data(1,:).*A_ges;

        % routed discharge
        qgko_spec(:,:) = qgko.data(2:size(qgko.data,1)-88,:);   % specific discharge [l/s/km^2]
        A_ezg_mat(:,:) = repmat(A_ezg(1,:),size(qgko_spec,1),1);
        qgko_abs_becken = qgko_spec.*A_ezg_mat/3.6;          % absolute discharge [m^3/s]

        qges_spec(:,:) = qges.data(2:size(qges.data,1),:);   % specific discharge [l/s/km^2]
        A_teilezg_mat(:,:) = repmat(A_teilezg(1,:),size(qges_spec,1),1);
        qges_abs_becken = qges_spec.*A_teilezg_mat/3.6;          % absolute discharge [m^3/s]
        
            if isempty(qgko_abs_becken)==1
                return
            end
            
        for p = 1:length(pegel)
            p;
             pegel(p).valid(:,1) = qgko_date_becken;
             pegel(p).valid(:,2) = qgko_abs_becken(:,p);
             pegel(p).valid(:,3) = qgko_abs_becken(:,p);
             pegel(p).valid(:,5) = qges_abs_becken(:,p);
        end
         
        % Beckendaten einlesen
        for i = 1:length(pegelnummer)
            pegel_no = pegelnummer(i);
            becken_no = beckennummer(i);
            hrb_name = dir(sprintf('sp*%d_%d_hrb.ddp',becken_no,pegel_no));
            speicher = importdata(strcat(pfad_valid_output,'\',hrb_name.name));
            pegel(pegel_no).valid(:,3) = speicher.data(:,10);
            pegel(pegel_no).valid(:,4) = speicher.data(:,7);
%             pegel(pegel_no).valid(:,6) = speicher.data(:,7);
        end
        if Gebiet==2
        
            speicher = importdata(strcat(pfad_valid_output,'\','sp1_teg.ddp,'));
            pegel(9).valid(:,3) = speicher.data(:,10);
            pegel(9).valid(:,4) = speicher.data(:,7);
            
            speicher = importdata(strcat(pfad_valid_output,'\','sp2_sch.ddp,'));
            pegel(10).valid(:,3) = speicher.data(:,10);
            pegel(10).valid(:,4) = speicher.data(:,7);
            
            speicher = importdata(strcat(pfad_valid_output,'\','sp3_hrb.ddp,'));
            pegel(11).valid(:,3) = speicher.data(:,10);
            pegel(11).valid(:,4) = speicher.data(:,7);
            
            speicher = importdata(strcat(pfad_valid_output,'\','sp4_hrb.ddp,'));
            pegel(12).valid(:,3) = speicher.data(:,10);
            pegel(12).valid(:,4) = speicher.data(:,7);

        end
    
   raus=1       
if raus==0

%% Berechnung der Wirksamkeitsverläufe
date_format = 'dd.mm.yyyy';

% regarded hydrograph time series
% wellen_start_line_becken = find(qgko_date_becken(:,1)==datenum(wellen_start_date,date_format));
% wellen_end_line_becken = find(qgko_date_becken(:,1)==datenum(wellen_end_date,date_format));
zeitreihe_becken = 1:size(qgko_date_becken,1);%wellen_start_line_becken:wellen_end_line_becken;

fun_PBIAS   = @(x,y) sum(x-y)*100/sum(x);
fun_NSE     = @(x,y) 1-sum((x-y).^2)/sum((x-mean(x)).^2);


for m = 1:length(Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf)
    startpegel_no = Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).nachfolger_beckenposition(1,1);
    for n = 1:length(Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).nachfolger_beckenposition(1,:))
    Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).wasim_abminderung(1,n) = ...
        (max(pegel(Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).nachfolger_beckenposition(1,n)).inflow_wasim(:,3))...
        -max(qgko_abs_becken(:,Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).nachfolger_beckenposition(1,n))))./...
        max(pegel(Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).nachfolger_beckenposition(1,n)).inflow_wasim(:,3));
    p = Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).nachfolger_beckenposition(1,n);
    gof_nse(p)= fun_NSE(qgko_abs_becken(zeitreihe_becken,p),pegel(p).outflow(:,2));
    gof_pbias(p)= fun_PBIAS(qgko_abs_becken(zeitreihe_becken,p),pegel(p).outflow(:,2));
    Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).gof_nse(1,n) = gof_nse(Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).nachfolger_beckenposition(1,n));
    Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).gof_pbias(1,n) = gof_pbias(Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).nachfolger_beckenposition(1,n));
    gof_dsa(p) =  Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).wasim_abminderung(1,n)- Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).nachfolger_abminderung(1,n);
    end
end

%%

for m = 1:length(Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf)
    Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).Dsa = abs(Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).nachfolger_abminderung-...
        Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).wasim_abminderung);
    Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).mean_gof_nse = mean(Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).gof_nse);
    Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).mean_gof_pbias = mean(Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).gof_pbias);
    Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).mean_Dsa = mean(Wirksamkeit_Beckenkombination(runde).wirksamkeitsverlauf(m).Dsa);
end

%% mittlere Gütekriterien über alle veränderten Beckenpunkte

p_m = findobj(pegel,'becken_yn',1,'-or','becken_vor_yn',1);

for i = 1:length(p_m)
    m = p_m(i).gage.no;
    sa_wasim_tmp(i) = (max(pegel(m).inflow_wasim(:,3))-max(qgko_abs_becken(:,m)))./max(pegel(m).inflow_wasim(:,3));
    sa_matlab_tmp(i) = (max(pegel(m).inflow_wasim(:,3))-max(pegel(m).outflow(:,2)))./max(pegel(m).inflow_wasim(:,3));
    dsa_tmp(i) = abs(sa_wasim_tmp(i)-sa_matlab_tmp(i));
    nse_tmp(i) = fun_NSE(qgko_abs_becken(:,m),pegel(m).outflow(:,2));
    pbias_tmp(i) = fun_PBIAS(qgko_abs_becken(:,m),pegel(m).outflow(:,2));   
end

Wirksamkeit_Beckenkombination(runde).gof.Dsa = (dsa_tmp);
Wirksamkeit_Beckenkombination(runde).gof.Dsa_mean = mean(dsa_tmp);
Wirksamkeit_Beckenkombination(runde).gof.Dsa_std = std(dsa_tmp);
Wirksamkeit_Beckenkombination(runde).gof.Dsa_var = var(dsa_tmp);
Wirksamkeit_Beckenkombination(runde).gof.Dsa_median = median(dsa_tmp);
Wirksamkeit_Beckenkombination(runde).gof.nse = mean(nse_tmp);
Wirksamkeit_Beckenkombination(runde).gof.pbias = mean(pbias_tmp);
end