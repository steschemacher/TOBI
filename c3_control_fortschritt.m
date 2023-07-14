%%% Datenkontrolle der Becken-Szenarien %%%
clearvars


path = 'PFAD';
path1 = strcat(path,'\TOBas');
path2 = strcat(path,'\Results');
path3 = strcat(path2,'\Sim_overview');
addpath(path);

name2 = {'daten_ges','info_beckenkombination','info_sim','lokal','regional','zielfunktion'};

E = [1:7]';

% Simulationsübersicht einlesen
cd(path3)
while exist('open_characteristics.txt','file')==2
    pause(5*rand(1))
end
save('open_characteristics.txt','E','-ascii')
load('characteristics.mat');
delete('open_characteristics.txt')

% Schleife über alle Simulationsblöcke
for sim_block = 1:size(characteristics,1)
    
    Gebiet = characteristics(sim_block,2);
    TG = characteristics(sim_block,3);
    mOpt = characteristics(sim_block,4);
    Eopt = characteristics(sim_block,5);
    Ekomb = abs(characteristics(sim_block,6));
    
    if mOpt == 1
        E_no_opt = [[Eopt,Eopt,Ekomb];...
            [setdiff(E,Eopt),repmat(Eopt,length(E)-1,1),repmat(Ekomb,length(E)-1,1)]];
    elseif mOpt == 2
        E_no_opt = [[Ekomb,Ekomb,Ekomb];...
            [setdiff(E,Ekomb),repmat(Ekomb,length(E)-1,2)];...
            [repmat(setdiff(E,Ekomb),1,2),repmat(Ekomb,length(E)-1,1)]];
    end
    
    
    if mOpt==1
        path4 = strcat(path2,'\',sprintf('G%d_TG%d_mOpt%d_komb%d',Gebiet,TG,mOpt,Ekomb));
    else
        path4 = strcat(path2,'\',sprintf('G%d_TG%d_mOpt%d',Gebiet,TG,mOpt));
    end
    cd(path4)
    
    % Berechnung der benötigten Rechenschritte
    if mOpt==1
        anzahl_bloecke = 10*size(E_no_opt,1);
    else
        anzahl_bloecke = 21*size(E_no_opt,1);
    end
    
    % Schleife über Eno, Eopt und Beckenanzahlen
    if mOpt==1
        schritte = 1:10;
    else
        i_all = [1:-.05:0];
        anteil1 = .15;
        anteil2 = (1-anteil1*2);
        schritte_val = fliplr(i_all.*anteil2);
        schritte = 1:length(schritte_val);
    end
    
    for n = 1:size(E_no_opt,1)
        E_no = E_no_opt(n,1);
        E_opt = E_no_opt(n,2);
        E_komb = E_no_opt(n,3);
        %     for E_no = E_NO
        %         for E_opt = EOPT
        for n_schritte = schritte
            
            % Kontrolle der Datenvollständigkeit (daten_ges, info_beckenkombination, info_sim, lokal, regional, zielfunktion)
            
            if mOpt==1
                name1 = sprintf('Eno%d_Eopt%d_komb%d_%02dnoBasins_',E_no,E_opt,E_komb,n_schritte);
                name3 = '';
            else
                name1 = sprintf('Eno%d_Eopt%d_kombE%d',E_no,E_opt,E_komb);
                name3 = sprintf('_%4.2f_',schritte_val(n_schritte));
            end
            for n_name2 = 1%:length(name2)
                if isempty(dir(strcat(name1,'*',name2{n_name2},name3,'*','.mat')))==0
                    daten_vollst{n_name2,sim_block}(n,n_schritte) = 1;
                else
                    daten_vollst{n_name2,sim_block}(n,n_schritte) = 0;
                end
            end
        end
%         end
    end
    
    % Berechnung des Rechenfortschritts
    rechenfortschritt(sim_block,1) = sum(sum(daten_vollst{1,sim_block}(:,:)))/anzahl_bloecke;
    fprintf('Simulationsblock %d: Gebiet %d, TG %d, mOpt %d, Eopt %d, Ekomb %d; Fortschritt: %4.2f \n',...
        sim_block,Gebiet,TG,mOpt,E_opt,Ekomb,rechenfortschritt(sim_block,1))
    
    
    
end