function [Beckenkombination,pegel_runde,pegel_becken] = f2_2__select_basins(pegel,Gebiet,TG,mOpt,beckenpunkte_select,opt_parameter)


if mOpt == 1
    
    becken_m = [];
    for i = 1:length(beckenpunkte_select)
        if any(TG==pegel(1, beckenpunkte_select(i).no_pegel).neighbors.alle_nachfolger)==1
            becken_m = [becken_m,beckenpunkte_select(i).no_becken];
        end
    end
    
    if length(becken_m)>10
        beckenanzahlen = ceil([length(becken_m)/10:floor(length(becken_m)/10):length(becken_m)]);
    else
        beckenanzahlen = 1:length(becken_m);
    end
    
    
    for i = 1:length(beckenanzahlen)
        no_of_basins = beckenanzahlen(i);
        
        clear wirksamkeit pegelstandort combinations Wirksamkeit_Beckenkombination
        
        possible_combinations = nchoosek(length(becken_m),no_of_basins)';
        if possible_combinations < opt_parameter
            % if no_of_basins==1
            combinations = nchoosek(becken_m,no_of_basins);
        else
            no_combinations = 1;
            while no_combinations<=opt_parameter
                combinations(no_combinations,:) = randsample(becken_m,no_of_basins);
                combinations = sort(combinations,2);
                combinations = unique(combinations,'rows');
                no_combinations = size(combinations,1)+1;
            end
        end
        Beckenkombination(i,1).becken = combinations';
    end
    
    pegel_runde = [];
    pegel_becken = [];
    
elseif mOpt==2
   %%
    clear Beckenkombination pegel_becken pegel_beckenanzahl pegel_delete pegel_delete_unterhalb pegel_no pegel_runde pegel_select pegel_vergessen pege_vergessen_beckenanzahl pegel_zusätzlich runde_einzelkombination runde_einzelkombination_combine
if Gebiet==1
    max_becken = 5;
else
    max_becken = 5;
end
min_becken = 2;
i = 1;

pegel_select = [TG,pegel(TG).neighbors.alle_vorgaenger];
pegel_ausschluss = setdiff(1:length(pegel),pegel_select);

pegel_becken{1,length(pegel)}=[];
pegel_beckenanzahl(1,length(pegel)) = 0;
pegel_no = 1:length(pegel);
for j = pegel_select
    pegel_becken{i,j} = [beckenpunkte_select(ismember([beckenpunkte_select.no_pegel],[j,[pegel(j).neighbors.alle_vorgaenger]])).no_becken];
    pegel_beckenanzahl(i,j) = length(pegel_becken{i,j});
%     pegel_no(1,j) = j;
end

% Pegelpunkte selektieren, die 3-5 Becken enthalten
pegel_runde{i} = pegel_no(1,setdiff(find(pegel_beckenanzahl(i,:)<=max_becken),find(pegel_beckenanzahl(i,:)<=min_becken)));

pegel_keinFortschritt = 1;  % dummy-Wert um die Schleife zu betreten
while isempty(pegel_keinFortschritt)==0
    % Pegelpunkte ausschließen, die in einem anderen Gebiet enthalten sind
    pegel_delete = [];
    for j = (pegel_runde{i})
        pegel_delete = [pegel_delete,[pegel(j).neighbors.alle_vorgaenger]];
    end
    pegel_delete = unique(pegel_delete);
    
    % Pegelpunkte und Becken für den jeweiligen Berechnungsschritt
    pegel_runde{i} = setdiff(pegel_runde{i},pegel_delete);
    
    for j = 1:length(pegel)
        pegel_becken{i+1,j} = pegel_becken{i,j};
        
        % Pegel ausschließen, die schon in vorheriger Runde ausgeschlossen
        % wurden, weil sie innerhalb des Gebiets eines anderen Pegels liegen
        if length(find(pegel_delete==j))==1
            pegel_becken{i+1,j} = [];
        end
        
        % Beckenkombinationen aus vorheriger Runde durch Pegelergebnis
        % ersetzen (*-1, um die Beckenkombination zuordnen zu können)
        if length(find(pegel_runde{i}==j))==1
            pegel_becken{i+1,j} = [];
        end
        
        pegel_beckenanzahl(i+1,j) = length(pegel_becken{i+1,j});
    end
    
    for j = 1:length(pegel)
        
        if length(find(pegel_runde{i}==j))==1
            % Beckenkombinationen aus Nachfolgerpegeln entfernen und
            % den Index des Pegelergebnisses ergänzen
            for k = setdiff([pegel(j).neighbors.alle_nachfolger],pegel_ausschluss)
                pegel_becken{i+1,k} = setdiff(pegel_becken{i,k},pegel_becken{i,j});
                pegel_becken{i+1,k} = [pegel_becken{i+1,k},j*(-1)];
                pegel_beckenanzahl(i+1,k) = length(pegel_becken{i+1,k});
            end
        end
        
    end
    
    pegel_keinFortschritt = find(pegel_beckenanzahl(i,:)-pegel_beckenanzahl(i+1,:)==0 & pegel_beckenanzahl(i,:)> max_becken);
    disp('Für Fortschritt kritische Pegelpunkte: ')
    disp(pegel_keinFortschritt)
    disp(pegel_beckenanzahl(i,[pegel_keinFortschritt]))

    if isempty(pegel_keinFortschritt)==0
        pegel_zusatz = [];
        for p = pegel_keinFortschritt
            if isempty([pegel(p).neighbors.alle_vorgaenger])==1
                pegel_zusatz = [pegel_zusatz,p];
                max_becken_neu = pegel_beckenanzahl(i,p);
                fprintf('Pegel %d mit %d Becken berücksichtigt.',p,max_becken_neu)
                pegel_runde{i} = pegel_no(1,setdiff(find(pegel_beckenanzahl(i,:)<=max_becken_neu),find(pegel_beckenanzahl(i,:)<=min_becken)));
            end
            for m = [pegel(p).neighbors.alle_vorgaenger]
                if pegel_beckenanzahl(i,m)==min_becken
                    pegel_zusatz = [pegel_zusatz,m];
                    fprintf('Pegel %d mit %d Becken berücksichtigt.',m,min_becken)
                end
            end
        end
        
        n=0;
        while isempty(pegel_zusatz)==1 && n<=2
            n=n+1;
            pegel_runde_neu = pegel_no(1,setdiff(find(pegel_beckenanzahl(i,:)<=max_becken+n),find(pegel_beckenanzahl(i,:)<=min_becken)));
            pegel_zusatz = setdiff(pegel_runde_neu,pegel_runde{i});
                fprintf('Mit %d Becken berücksichtigt: ',max_becken+n)
                disp(pegel_zusatz)
            pegel_runde{i} = pegel_runde_neu;
        end
        if isempty(pegel_zusatz)==1
            for p = pegel_keinFortschritt
                for m = [pegel(p).neighbors.alle_vorgaenger]
                    if pegel_beckenanzahl(i,m)==1
                        pegel_zusatz = [pegel_zusatz,m];
                        fprintf('Pegel %d mit %d Becken berücksichtigt.',m,1)
                    end
                end
            end
        end
        pegel_zusatz = unique(pegel_zusatz);
        if isempty(pegel_zusatz)==1
            disp('Problem1!!!')
        end
        
        pegel_runde{i} = [pegel_runde{i},pegel_zusatz];
        
        % Pegelpunkte ausschließen, die in einem anderen Gebiet enthalten sind
        pegel_delete = [];
        for j = (pegel_runde{i})
            pegel_delete = [pegel_delete,[pegel(j).neighbors.alle_vorgaenger]];
        end
        pegel_delete = unique(pegel_delete);
        
        % Pegelpunkte und Becken für den jeweiligen Berechnungsschritt
        pegel_runde{i} = setdiff(pegel_runde{i},pegel_delete);
    end
    
end

% clear becken_runde
% becken_runde(1,:) = pegel_beckenanzahl(i,[pegel_runde{i}]);
% k=0;
% for j = pegel_runde{i}
%     k=k+1;
%     pegel_nachfolger = [pegel(j).neighbors.alle_nachfolger];
%     becken_runde(2:length(pegel_nachfolger)+1,k) = pegel_beckenanzahl(i,pegel_nachfolger)';
% end
% 
% pegel_all = pegel_runde{i};
% % Pegel bestimmen, die unterhalb liegen
% pegel_delete_unterhalb = [];
% for j = [pegel_all]
%     pegel_delete_unterhalb = [pegel_delete_unterhalb,[pegel(j).neighbors.alle_nachfolger]];
% end
% pegel_delete_unterhalb = unique(pegel_delete_unterhalb);


%
while length(pegel_runde{i})>0
    %
    clear runde_einzelkombinationen runde_einzelkombinationen_combine kombination combinations
    % Bestimmung aller Kombinationen
    for j = 1:length(pegel_runde{i})
        j_val = pegel_runde{i}(j);
        for k = 0:length(pegel_becken{i,j_val})
            runde_einzelkombinationen{j,k+1} = nchoosek([pegel_becken{i,j_val}],k);
        end
    end
    
    % Kombination der (sich nicht überschneidenen) Kombinationen
    combinations = nan(size(runde_einzelkombinationen,2)*size(runde_einzelkombinationen,1),2^size(runde_einzelkombinationen,2));
    
    for j = 1:size(runde_einzelkombinationen,1)
        n_all = 0;
        for m = 1:size(runde_einzelkombinationen,2)
            if isempty(runde_einzelkombinationen{j,m})==1
%                 runde_einzelkombinationen{j,m} = NaN;
            end
            for n = 1:size(runde_einzelkombinationen{j,m},1)
                n_all = n_all+1;
                runde_einzelkombinationen_combine{j,n_all} = runde_einzelkombinationen{j,m}(n,:);
            end
        end
    end
    
    max_lines = 0;
    for n_all = 1:size(runde_einzelkombinationen_combine,2)
        kombination = [runde_einzelkombinationen_combine{:,n_all}];
%         if isempty(kombination)==0
        combinations(1:length(kombination),n_all) = kombination;
        max_lines = max(max_lines,size(combinations(isnan(combinations(:,n_all))==0,n_all),1));
%         end
    end
    
    Beckenkombination(i,1).becken = combinations(1:max_lines,:);
    %
    if isempty(find(pegel_runde{i}==TG))==0
        break
    end
    
    i = i+1;
    
    for j = 1:length(pegel)
        pegel_becken{i,j} = pegel_becken{i-1,j};
        
        % Pegel ausschließen, die schon in vorheriger Runde ausgeschlossen
        % wurden, weil sie innerhalb des Gebiets eines anderen Pegels liegen
        if length(find(pegel_delete==j))==1
            pegel_becken{i,j} = [];
        end
        
        % Beckenkombinationen aus vorheriger Runde durch Pegelergebnis
        % ersetzen (*-1, um die Beckenkombination zuordnen zu können)
        if length(find(pegel_runde{i-1}==j))==1
            pegel_becken{i,j} = [];
        end
        
        pegel_beckenanzahl(i,j) = length(pegel_becken{i,j});
    end
    
    for j = 1:length(pegel)
        
        if length(find(pegel_runde{i-1}==j))==1
            % Beckenkombinationen aus Nachfolgerpegeln entfernen und
            % den Index des Pegelergebnisses ergänzen
            for k = setdiff([pegel(j).neighbors.alle_nachfolger],pegel_ausschluss)
                pegel_becken{i,k} = setdiff(pegel_becken{i,k},pegel_becken{i-1,j});
                pegel_becken{i,k} = [pegel_becken{i,k},j*(-1)];
                pegel_beckenanzahl(i,k) = length(pegel_becken{i,k});
            end
        end
        
    end
    
    % Pegelpunkte selektieren, die 3-5 Becken enthalten
    pegel_runde{i} = pegel_no(1,setdiff(find(pegel_beckenanzahl(i,:)<=max_becken),find(pegel_beckenanzahl(i,:)<=min_becken)));
    max_exist_becken = max(pegel_beckenanzahl(i,:));
    if max_exist_becken<=min_becken
        pegel_runde{i} = pegel_no(1,find(pegel_beckenanzahl(i,:)==max_exist_becken));
        
    end
    
%     if length(pegel_runde{i})==1 && isempty([pegel(pegel_runde{i}).neighbors.alle_nachfolger])==1
%         
%     end


    
    pegel_keinFortschritt = 1;  % dummy-Wert um die Schleife zu betreten
    while isempty(pegel_keinFortschritt)==0
        % Pegelpunkte ausschließen, die in einem anderen Gebiet enthalten sind
        pegel_delete = [];
        for j = (pegel_runde{i})
            pegel_delete = [pegel_delete,[pegel(j).neighbors.alle_vorgaenger]];
        end
        pegel_delete = unique(pegel_delete);
        
        % Pegelpunkte und Becken für den jeweiligen Berechnungsschritt
        pegel_runde{i} = setdiff(pegel_runde{i},pegel_delete);
        
        for j = 1:length(pegel)
            pegel_becken{i+1,j} = pegel_becken{i,j};
            
            % Pegel ausschließen, die schon in vorheriger Runde ausgeschlossen
            % wurden, weil sie innerhalb des Gebiets eines anderen Pegels liegen
            if length(find(pegel_delete==j))==1
                pegel_becken{i+1,j} = [];
            end
            
            % Beckenkombinationen aus vorheriger Runde durch Pegelergebnis
            % ersetzen (*-1, um die Beckenkombination zuordnen zu können)
            if length(find(pegel_runde{i}==j))==1
                pegel_becken{i+1,j} = [];
            end
            
            pegel_beckenanzahl(i+1,j) = length(pegel_becken{i+1,j});
        end
            
        if isempty(find(pegel_runde{i}==TG))==0
            break
        end
        
        for j = 1:length(pegel)
            if length(find(pegel_runde{i}==j))==1
                % Beckenkombinationen aus Nachfolgerpegeln entfernen und
                % den Index des Pegelergebnisses ergänzen
                for k = setdiff([pegel(j).neighbors.alle_nachfolger],pegel_ausschluss)
                    pegel_becken{i+1,k} = setdiff(pegel_becken{i,k},pegel_becken{i,j});
                    pegel_becken{i+1,k} = [pegel_becken{i+1,k},j*(-1)];
                    pegel_beckenanzahl(i+1,k) = length(pegel_becken{i+1,k});
                end
            end
        end
        
        pegel_keinFortschritt = find(pegel_beckenanzahl(i,:)-pegel_beckenanzahl(i+1,:)==0 & pegel_beckenanzahl(i,:)> max_becken);
        disp('Für Fortschritt kritische Pegelpunkte: ')
        disp(pegel_keinFortschritt)
        disp(pegel_beckenanzahl(i,[pegel_keinFortschritt]))
        
        if isempty(pegel_keinFortschritt)==0
            pegel_zusatz = [];
            for p = pegel_keinFortschritt
                if isempty([pegel(p).neighbors.alle_vorgaenger])==1
                    pegel_zusatz = [pegel_zusatz,p];
                    max_becken_neu = pegel_beckenanzahl(i,p);
                    fprintf('Pegel %d mit %d Becken berücksichtigt.',p,max_becken_neu)
                    pegel_runde{i} = pegel_no(1,setdiff(find(pegel_beckenanzahl(i,:)<=max_becken_neu),find(pegel_beckenanzahl(i,:)<=min_becken)));
                end
                for m = [pegel(p).neighbors.alle_vorgaenger]
                    if pegel_beckenanzahl(i,m)==min_becken
                        pegel_zusatz = [pegel_zusatz,m];
                        fprintf('Pegel %d mit %d Becken berücksichtigt.',m,min_becken)
                    end
                end
            end
            
            
            n=0;
            while isempty(pegel_zusatz)==1 && n+max_becken<=7
                n=n+1;
                pegel_runde_neu = pegel_no(1,setdiff(find(pegel_beckenanzahl(i,:)<=max_becken+n),find(pegel_beckenanzahl(i,:)<=min_becken)));
                pegel_zusatz = setdiff(pegel_runde_neu,pegel_runde{i});
                fprintf('Mit %d Becken berücksichtigt: ',max_becken+n)
                disp(pegel_zusatz)
                pegel_runde{i} = pegel_runde_neu;

            end
            if isempty(pegel_zusatz)==1
                for p = pegel_keinFortschritt
                    for m = [pegel(p).neighbors.alle_vorgaenger]
                        if pegel_beckenanzahl(i,m)==1
                            pegel_zusatz = [pegel_zusatz,m];
                            fprintf('Pegel %d mit %d Becken berücksichtigt.',m,1)
                        end
                    end
                end
            end
            pegel_zusatz = unique(pegel_zusatz);
            if isempty(pegel_zusatz)==1
                fprintf('Problem2!!!, %d',pegel_keinFortschritt)
            end
            
            pegel_runde{i} = [pegel_runde{i},pegel_zusatz];
            
            % Pegelpunkte ausschließen, die in einem anderen Gebiet enthalten sind
            pegel_delete = [];
            for j = (pegel_runde{i})
                pegel_delete = [pegel_delete,[pegel(j).neighbors.alle_vorgaenger]];
            end
            pegel_delete = unique(pegel_delete);
            
            % Pegelpunkte und Becken für den jeweiligen Berechnungsschritt
            pegel_runde{i} = setdiff(pegel_runde{i},pegel_delete);
        end
    end
    

    
end

% prüfen, ob jedes Becken berücksichtigt wurde
for becken_no = [beckenpunkte_select(ismember([beckenpunkte_select.no_pegel],[pegel(TG).neighbors.alle_vorgaenger])).no_becken]
    for i = 1:length(Beckenkombination)
        anzahl(i) = length(find([Beckenkombination(i).becken]==becken_no));
    end
    if sum(anzahl)==0
        fprintf('Becken %d wird nicht berücksichtigt.\n',becken_no)
    end
end
beckenanzahlen = [];
    
elseif mOpt==10 || mOpt==11
    Beckenkombination(1,1).becken = opt_parameter';
    beckenanzahlen = length(opt_parameter);
    pegel_runde = [];
    pegel_becken = [];
elseif mOpt==12
    becken_m = [];
    pegel_m = [];
    for i = 1:length(beckenpunkte_select)
        if any(TG==pegel(beckenpunkte_select(i).no_pegel).neighbors.alle_nachfolger)==1
            if isempty([pegel(beckenpunkte_select(i).no_pegel).neighbors.vorgaenger])==0
                if isempty(find(pegel_m==[beckenpunkte_select(i).no_pegel]))==1
                    becken_m = [becken_m,beckenpunkte_select(i).no_becken];
                    pegel_m = [pegel_m,beckenpunkte_select(i).no_pegel];
                end
            end
        end
    end
    
    beckenanzahlen = opt_parameter;
    
    for i = 1:length(beckenanzahlen)
        no_of_basins = beckenanzahlen(i);
        
        clear wirksamkeit pegelstandort combinations Wirksamkeit_Beckenkombination
        
        possible_combinations = nchoosek(length(becken_m),no_of_basins)';
        if possible_combinations < opt_parameter
            % if no_of_basins==1
            combinations = nchoosek(becken_m,no_of_basins);
        else
            no_combinations = 1;
            while no_combinations<=opt_parameter
                combinations(no_combinations,:) = randsample(becken_m,no_of_basins);
                combinations = sort(combinations,2);
                combinations = unique(combinations,'rows');
                no_combinations = size(combinations,1)+1;
            end
        end
        Beckenkombination(i,1).becken = combinations';
    end
    
    pegel_runde = [];
    pegel_becken = [];
    
elseif mOpt == 13
    pegel_temp = findobj(pegel,'pegel_vor_yn',0); % basin at the gage but not before
    no_temp = [];
    for i = 1:length(pegel_temp)
        no_temp = [no_temp,pegel_temp(i).gage.no];
    end
    Beckenkombination(1).becken(:,1) = no_temp;
    
    beckenanzahlen = [];
    pegel_runde = [];
    pegel_becken = [];
    
end