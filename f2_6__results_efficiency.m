function [info_sim,info_beckenkombination,zielfunktion,daten_ges,lokal,regional] = f2_6__results_efficiency(qgko_date,becken_pos,max_becken_pos,beobachtungspunkte,beckenpunkte_select2,pegel,mOpt)

%%% Funktion zur Abspeicherung der Simulationsergebnisse
%
% info_sim:
%   1: Anfangszeitpunkt
%   2: Endzeitpunkt
%   3: Anzahl der Becken
%
% info_beckenkombination:
%   1: Pegelnummern
%   2: Beckennummern
%
% Zielfunktionen:
%   1: Scheitelabminderung regional
%   2: Scheiteländerung absolut, regional
%   3: Scheitelabfluss [m³/s]
%   4: Abflussanteil (Scheitelabfluss abgemindert/Scheitelabfluss original)
%   5: mittlere Beckenfüllung
%   6: Beckenanzahl
%   7: Überflutungsfläche
%   8: Flächenanteil
%   9: Dammvolumen
%   10: Damvolumen/EZG-Fläche
%   11: Beckenvolumen (potentiell HRB)
%   12: sV (HRB)
%   13: Beckenvolumen (genutzt ohne HWE)
%   14: sV (genutzt ohne HWE)
%   15: Beckenvolumen (genutzt mit HWE)
%   16: sV (genutzt mit HWE)
%   17: Beckenfüllung normiert (>1 --> 0)
%   18: Beckenzahl normiert (maximales n aller pot. Becken)
%   19: Überflutungsfläche normiert (maximales A aller pot. Becken)
%   20: Dammvolumen normiert (maximales V aller pot. Becken)
%   21: Beckenfüllung
%   22: Maximaler Abfluss (WaSiM)
%   23: zurückgehaltenes Volumen
%
% daten_ges:
%   1: Pegelnummer
%   2: TG-Größe
%   3: gesamtes gespeichertes Volumen MIT HWE
%   4: gesamtes gespeichertes Volumen OHNE HWE
%   5: potentielles gesamtes Volumen OHNE HWE
%   6: Scheitelabminderung
%   7: Anstiegsvolumen (Zufluss WaSiM)
%   8: Anstiegsvolumen (Zufluss MATLAB)
%   9: Anstiegsvolumen (Abfluss MATLAB)
%   10: Anstiegsdauer (Zufluss WaSiM)
%   11: Anstiegsdauer (Zufluss MATLAB)
%   12: Anstiegsdauer (Abfluss MATLAB)
%   13: Scheitel Zufluss WaSiM
%   14: Scheitel Zufluss MATLAB
%   15: Scheitel Abfluss MATLAB
%   16: Becken (ja/nein - 1/0)
%
% lokal:
% 	1: Pegelnummer
% 	2: Beckennummer
% 	3: xcoord
% 	4: ycoord
% 	5: Drosselöffnungsgröße
% 	6: maximaler Freispiegelabfluss
% 	7: maximaler Druckabfluss
% 	8: maximales Beckenvolumen (geplant)
% 	9: maximales Beckenvolumen (ist)
% 	10:	Hochwasserentlastung notwendig (1/0)
% 	11:	Becken gefüllt (1/0)
% 	12:	Füllung des Beckens (1: gefüllt)
% 	13:	Gebietsgröße
% 	14:	Beckenvolumen inkl. oberstrom (mit HWE)
% 	15:	Beckenvolumen inkl. oberstrom (ohne HWE)
% 	16:	pot. Beckenvolumen inkl. oberstrom
% 	17:	spec. Vol des Beckens inkl oberstrom (ohne HWE)
% 	18:	Becken oberstrom (1: ja, 0: nein)
% 	19:	Maximum (Zufluss WaSiM)
% 	20:	Maximum (Zufluss MATLAB)
% 	21:	Maximum (Abfluss MATLAB)
% 	22:	Abminderung am Becken
% 	23:	lokale Abminderung am Becken
% 	24:	Anteil der Becken-Scheitelreduktion am Gebietsauslassabfluss
% 	25:	Anstiegsvolumen (Zufluss WaSiM)
% 	26:	Anstiegsvolumen (Zufluss MATLAB)
% 	27:	Anstiegsvolumen (Abfluss MATLAB)
% 	28:	Spitze (QmaxF bis .9*Qmax)
% 	29:	Spitze (QmaxF bis .8*Qmax)
% 	30:	Spitze (QmaxF bis .7*Qmax)
% 	31:	Spitze (QmaxF bis .6*Qmax)
% 	32:	Spitze (QmaxF bis .5*Qmax)
% 	33:	Spitze (QmaxF bis .4*Qmax)
% 	34:	Spitze (QmaxF bis .3*Qmax)
% 	35:	Spitze (QmaxF bis .2*Qmax)
% 	36:	Spitze (QmaxF bis .1*Qmax)
% 	37:	Dauer (Anstieg Zufluss WaSiM)
% 	38:	Dauer (Anstieg Zufluss MATLAB)
% 	39:	Dauer (Anstieg Abfluss MATLAB)
% 	40:	Anstiegsvolumen (Zufluss WaSiM über Freispiegel)
% 	41:	Anstiegsvolumen (Zufluss MATLAB über Freispiegel)
% 	42:	Anstiegsvolumen (Abfluss MATLAB über Freispiegel)
% 	43:	Dauer (Anstieg Zufluss WaSiM über Freispiegel)
% 	44:	Dauer (Anstieg Zufluss MATLAB über Freispiegel)
% 	45:	Dauer (Anstieg Abfluss MATLAB über Freispiegel)
%   46: Einstaudauer
%   47: zurückgehaltenes Volumen
%
% regional:
%	1: Pegelnummer
%	2: Beckenanzahl im TG
%	3: Gebietsgröße
%	4: Beckenvolumen inkl. oberstrom (mit HWE)
%	5: Beckenvolumen inkl. oberstrom (ohne HWE)
%	6: pot. Beckenvolumen inkl. oberstrom
%	7: spez. Vol (mit HWE)
%	8: spez. Vol (ohne HWE)
%	9: spez. Vol (pot. Beckenvol)
%	10: Reihenanordnung


runde = 1;
if max_becken_pos==0
    max_becken_pos = 1;
end

%% Informationen zur Beckenkombination
info_sim(runde,1) = qgko_date(1,1);
info_sim(runde,2) = qgko_date(end,1);
info_sim(runde,3) = size(becken_pos,1);

info_beckenkombination(runde,max_becken_pos,1:2) = nan;
info_beckenkombination(runde,1:length(becken_pos(:,1)),1) = becken_pos(:,1);
info_beckenkombination(runde,1:length(becken_pos(:,2)),2) = becken_pos(:,2);

%% Zielfunktionen für Optimierung
bemessung = 3;

for n = 1:length(pegel)
    zielfunktion(runde,n,1) = (max(pegel(n).inflow_wasim(:,3))-max(pegel(n).outflow{end}(:,2)))/max(pegel(n).inflow_wasim(:,3));
    zielfunktion(runde,n,2) = max(pegel(n).inflow_wasim(:,3))-max(pegel(n).outflow{end}(:,2));
    zielfunktion(runde,n,3) = max(pegel(n).outflow{end}(:,2));
    zielfunktion(runde,n,4) = max(pegel(n).outflow{end}(:,2))/max(pegel(n).inflow_wasim(:,3));
    becken_oberhalb = becken_pos(ismember(becken_pos(:,1),[n,pegel(n).neighbors.alle_vorgaenger]),:);
    zielfunktion(runde,n,5) = mean([pegel(becken_oberhalb(:,1)).beckenfuellung]);
    zielfunktion(runde,n,6) = size(becken_oberhalb,1); % geändert am 3.5.2019 (voher wurde bei Kombinationen mit 1 Becken 2 Becken angenommen)
    if zielfunktion(runde,n,6)~=0
        zielfunktion(runde,n,7) = sum([beckenpunkte_select2(find(ismember([beckenpunkte_select2.no_becken],becken_oberhalb(:,2))==1)).flaeche_see_hdam]);
        zielfunktion(runde,n,8) = zielfunktion(runde,n,7)*1e-6/pegel(n).gage.area;
        zielfunktion(runde,n,9) = sum([beckenpunkte_select2(find(ismember([beckenpunkte_select2.no_becken],becken_oberhalb(:,2))==1)).volumen_dam_hdam]);
        zielfunktion(runde,n,10) = zielfunktion(runde,n,9)/pegel(n).gage.area;
        V_HRB_tmp = sum([beckenpunkte_select2(find(ismember([beckenpunkte_select2.no_becken],becken_oberhalb(:,2))==1)).V_HRB],2);
        if mOpt~=13
            zielfunktion(runde,n,11) = V_HRB_tmp(bemessung);
        else
            zielfunktion(runde,n,11) = 0;
        end
        zielfunktion(runde,n,12) = zielfunktion(runde,n,11)/pegel(n).gage.area*1e-3;
        zielfunktion(runde,n,13) = pegel(n).gage.S_gesamt(end,2);
        zielfunktion(runde,n,14) = pegel(n).gage.S_gesamt(end,2)/pegel(n).gage.area*1e-3;
        zielfunktion(runde,n,15) = pegel(n).gage.S_gesamt(end,1);
        zielfunktion(runde,n,16) = pegel(n).gage.S_gesamt(end,1)/pegel(n).gage.area*1e-3;
        zielfunktion(runde,n,21) = zielfunktion(runde,n,15)./zielfunktion(runde,n,11);
    else
        zielfunktion(runde,n,5) = 1;
        zielfunktion(runde,n,7) = 0;
        zielfunktion(runde,n,8) = 0;
        zielfunktion(runde,n,9) = 0;
        zielfunktion(runde,n,10) = 0;
        zielfunktion(runde,n,11) = 0;
        zielfunktion(runde,n,12) = 0;
        zielfunktion(runde,n,13) = 0;
        zielfunktion(runde,n,14) = 0;
        zielfunktion(runde,n,15) = 0;
        zielfunktion(runde,n,16) = 0;
        zielfunktion(runde,n,21) = 1;
    end
    zielfunktion(runde,n,17) = zielfunktion(runde,n,15)/zielfunktion(runde,n,11);
    if zielfunktion(runde,n,17)>1; zielfunktion(runde,n,17)=-1; end
    zielfunktion(runde,n,18) = zielfunktion(runde,n,6)/sum(ismember([beckenpunkte_select2.no_pegel],[pegel(n).neighbors.alle_vorgaenger]));
    zielfunktion(runde,n,19) = zielfunktion(runde,n,7)/sum([beckenpunkte_select2(ismember([beckenpunkte_select2.no_pegel],[pegel(n).neighbors.alle_vorgaenger])).flaeche_see_hdam]);
    zielfunktion(runde,n,20) = zielfunktion(runde,n,9)/sum([beckenpunkte_select2(ismember([beckenpunkte_select2.no_pegel],[pegel(n).neighbors.alle_vorgaenger])).volumen_dam_hdam]);
    zielfunktion(runde,n,22) = max(pegel(n).inflow_wasim(:,3));
    
    vol_diff = pegel(n).inflow_wasim(:,3)-pegel(n).outflow{end}(:,2);
    lokal(runde,n,23) = sum(vol_diff(vol_diff>0))*3600; % zurückgehaltenes Volumen
end


%% Datenzusammenstellung für alle Pegelpunkte
for m = 1:length(pegel)
    daten_ges(runde,m,1) = pegel(m).gage.no;
    daten_ges(runde,m,2) = pegel(m).gage.area;
    daten_ges(runde,m,3) = pegel(m).gage.S_gesamt(end,1);
    daten_ges(runde,m,4) = pegel(m).gage.S_gesamt(end,2);
    daten_ges(runde,m,5) = pegel(m).gage.S_gesamt(end,3);
    
    zeitraum_inflowAW = 1:find(pegel(m).inflow_wasim(:,3)==max(pegel(m).inflow_wasim(:,3)),1);
    zeitraum_inflowA = 1:find(pegel(m).inflow{end}(:,2)==max(pegel(m).inflow{end}(:,2)),1);
    zeitraum_outflowA = 1:find(pegel(m).outflow{end}(:,2)==max(pegel(m).outflow{end}(:,2)),1);
    
    daten_ges(runde,m,7) = sum(pegel(m).inflow_wasim(zeitraum_inflowAW,3))*3600;
    daten_ges(runde,m,8) = sum(pegel(m).inflow{end}(zeitraum_inflowA,2))*3600;
    daten_ges(runde,m,9) = sum(pegel(m).outflow{end}(zeitraum_outflowA,2))*3600;
    daten_ges(runde,m,10) = length(zeitraum_inflowAW);
    daten_ges(runde,m,11) = length(zeitraum_inflowA);
    daten_ges(runde,m,12) = length(zeitraum_outflowA);
    daten_ges(runde,m,13) = max(pegel(m).inflow_wasim(:,3));
    daten_ges(runde,m,14) = max(pegel(m).inflow{end}(:,2));
    daten_ges(runde,m,15) = max(pegel(m).outflow{end}(:,2));
    daten_ges(runde,m,16) = pegel(m).becken_yn;
    daten_ges(runde,m,6) = (daten_ges(runde,m,13)-daten_ges(runde,m,15))/daten_ges(runde,m,13);
end

%% Datenzusammenstellung für alle Beckenpunkte
lokal(runde,max_becken_pos,1:45) = nan;

if mOpt ~= 13
for j = 1:size(becken_pos,1)
    lokal(runde,j,1) = becken_pos(j,1);
    lokal(runde,j,2) = becken_pos(j,2);
    spalte = find([pegel(becken_pos(j,1)).beckenparameter.beckennummer]==becken_pos(j,2));
    lokal(runde,j,3) = pegel(becken_pos(j,1)).gage.coordinates(1,1);
    lokal(runde,j,4) = pegel(becken_pos(j,1)).gage.coordinates(1,2);
    lokal(runde,j,5) = pegel(becken_pos(j,1)).beckenparameter(spalte).A_L;
    lokal(runde,j,6) = pegel(becken_pos(j,1)).beckenparameter(spalte).Q_maxF;
    lokal(runde,j,7) = pegel(becken_pos(j,1)).beckenparameter(spalte).Q_maxL;
    lokal(runde,j,8) = [pegel(becken_pos(j,1)).beckenparameter(spalte).S_max]';
    lokal(runde,j,9) = [pegel(becken_pos(j,1)).becken(spalte).S_max]';
    if lokal(runde,j,9) > lokal(runde,j,8) % ist > geplant
        lokal(runde,j,10) = 1;   % HWE notwendig
        lokal(runde,j,11) = 1; % Becken gefüllt
    elseif lokal(runde,j,9) < lokal(runde,j,8) % ist < geplant
        lokal(runde,j,10) = 0;
        lokal(runde,j,11) = 0;
    else
        lokal(runde,j,10) = 0;
        lokal(runde,j,11) = 1;
    end
    
    lokal(runde,j,12) = pegel(becken_pos(j,1)).beckenfuellung(spalte);
    lokal(runde,j,13) = pegel(becken_pos(j,1)).gage.area;
    lokal(runde,j,14) = pegel(becken_pos(j,1)).gage.S_gesamt(spalte,1);
    lokal(runde,j,15) = pegel(becken_pos(j,1)).gage.S_gesamt(spalte,2);
    lokal(runde,j,16) = pegel(becken_pos(j,1)).gage.S_gesamt(spalte,3);
    lokal(runde,j,17) = lokal(runde,j,15)/lokal(runde,j,13)*1e-3;
    lokal(runde,j,18) = pegel(becken_pos(j,1)).becken_vor_yn;
    
    lokal(runde,j,19) = max(pegel(becken_pos(j,1)).inflow_wasim(:,3));
    lokal(runde,j,20) = max(pegel(becken_pos(j,1)).inflow{spalte}(:,2));
    lokal(runde,j,21) = max(pegel(becken_pos(j,1)).outflow{spalte}(:,2));
    
    lokal(runde,j,22) = (lokal(runde,j,19)-lokal(runde,j,21))/lokal(runde,j,19);
    lokal(runde,j,23) = (lokal(runde,j,20)-lokal(runde,j,21))/lokal(runde,j,20);
    lokal(runde,j,24) = (lokal(runde,j,20)-lokal(runde,j,21))/lokal(runde,j,19);
    
    zeitraum_inflowAW = 1:find(pegel(becken_pos(j,1)).inflow_wasim(:,3)==max(pegel(becken_pos(j,1)).inflow_wasim(:,3)),1);
    zeitraum_inflowA = 1:find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)==max(pegel(becken_pos(j,1)).inflow{spalte}(:,2)),1);
    zeitraum_outflowA = 1:find(pegel(becken_pos(j,1)).outflow{spalte}(:,2)==max(pegel(becken_pos(j,1)).outflow{spalte}(:,2)),1);
    
    lokal(runde,j,25) = sum(pegel(becken_pos(j,1)).inflow_wasim(zeitraum_inflowAW,3))*3600;
    lokal(runde,j,26) = sum(pegel(becken_pos(j,1)).inflow{spalte}(zeitraum_inflowA,2))*3600;
    lokal(runde,j,27) = sum(pegel(becken_pos(j,1)).outflow{spalte}(zeitraum_outflowA,2))*3600;
    
    if lokal(runde,j,6)<lokal(runde,j,20)
        f = .9;
        if lokal(runde,j,6)<f*lokal(runde,j,20)
            clear diagonale
            diagonale = [repmat(lokal(runde,j,6),1,find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')-2),...
                lokal(runde,j,6):(f*lokal(runde,j,20)-lokal(runde,j,6))/(find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>f*lokal(runde,j,20),1,'last')-find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')+2):f*lokal(runde,j,20)];%,...
            diagonale = [diagonale,repmat(f*lokal(runde,j,20),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2))-length(diagonale))];
        else
            diagonale = repmat(lokal(runde,j,6),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        end
        diagonale = diagonale(1:length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        diff = pegel(becken_pos(j,1)).inflow{spalte}(:,2)-diagonale';
        diff(diff<0)=0;
        lokal(runde,j,28) = sum(diff)*3600;
        
        f = .8;
        if lokal(runde,j,6)<f*lokal(runde,j,20)
            clear diagonale
            diagonale = [repmat(lokal(runde,j,6),1,find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')-2),...
                lokal(runde,j,6):(f*lokal(runde,j,20)-lokal(runde,j,6))/(find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>f*lokal(runde,j,20),1,'last')-find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')+2):f*lokal(runde,j,20)];%,...
            diagonale = [diagonale,repmat(f*lokal(runde,j,20),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2))-length(diagonale))];
        else
            diagonale = repmat(lokal(runde,j,6),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        end
        diagonale = diagonale(1:length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        diff = pegel(becken_pos(j,1)).inflow{spalte}(:,2)-diagonale';
        diff(diff<0)=0;
        lokal(runde,j,29) = sum(diff)*3600;
        
        f = .7;
        if lokal(runde,j,6)<f*lokal(runde,j,20)
            clear diagonale
            diagonale = [repmat(lokal(runde,j,6),1,find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')-2),...
                lokal(runde,j,6):(f*lokal(runde,j,20)-lokal(runde,j,6))/(find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>f*lokal(runde,j,20),1,'last')-find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')+2):f*lokal(runde,j,20)];%,...
            diagonale = [diagonale,repmat(f*lokal(runde,j,20),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2))-length(diagonale))];
        else
            diagonale = repmat(lokal(runde,j,6),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        end
        diagonale = diagonale(1:length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        diff = pegel(becken_pos(j,1)).inflow{spalte}(:,2)-diagonale';
        diff(diff<0)=0;
        lokal(runde,j,30) = sum(diff)*3600;
        
        f = .6;
        if lokal(runde,j,6)<f*lokal(runde,j,20)
            clear diagonale
            diagonale = [repmat(lokal(runde,j,6),1,find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')-2),...
                lokal(runde,j,6):(f*lokal(runde,j,20)-lokal(runde,j,6))/(find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>f*lokal(runde,j,20),1,'last')-find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')+2):f*lokal(runde,j,20)];%,...
            diagonale = [diagonale,repmat(f*lokal(runde,j,20),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2))-length(diagonale))];
        else
            diagonale = repmat(lokal(runde,j,6),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        end
        diagonale = diagonale(1:length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        diff = pegel(becken_pos(j,1)).inflow{spalte}(:,2)-diagonale';
        diff(diff<0)=0;
        lokal(runde,j,31) = sum(diff)*3600;
        
        f = .5;
        if lokal(runde,j,6)<f*lokal(runde,j,20)
            clear diagonale
            diagonale = [repmat(lokal(runde,j,6),1,find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')-2),...
                lokal(runde,j,6):(f*lokal(runde,j,20)-lokal(runde,j,6))/(find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>f*lokal(runde,j,20),1,'last')-find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')+2):f*lokal(runde,j,20)];%,...
            diagonale = [diagonale,repmat(f*lokal(runde,j,20),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2))-length(diagonale))];
        else
            diagonale = repmat(lokal(runde,j,6),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        end
        diagonale = diagonale(1:length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        diff = pegel(becken_pos(j,1)).inflow{spalte}(:,2)-diagonale';
        diff(diff<0)=0;
        lokal(runde,j,32) = sum(diff)*3600;
        
        f = .4;
        if lokal(runde,j,6)<f*lokal(runde,j,20)
            clear diagonale
            diagonale = [repmat(lokal(runde,j,6),1,find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')-2),...
                lokal(runde,j,6):(f*lokal(runde,j,20)-lokal(runde,j,6))/(find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>f*lokal(runde,j,20),1,'last')-find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')+2):f*lokal(runde,j,20)];%,...
            diagonale = [diagonale,repmat(f*lokal(runde,j,20),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2))-length(diagonale))];
        else
            diagonale = repmat(lokal(runde,j,6),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        end
        diagonale = diagonale(1:length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        diff = pegel(becken_pos(j,1)).inflow{spalte}(:,2)-diagonale';
        diff(diff<0)=0;
        lokal(runde,j,33) = sum(diff)*3600;
        
        f = .3;
        if lokal(runde,j,6)<f*lokal(runde,j,20)
            clear diagonale
            diagonale = [repmat(lokal(runde,j,6),1,find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')-2),...
                lokal(runde,j,6):(f*lokal(runde,j,20)-lokal(runde,j,6))/(find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>f*lokal(runde,j,20),1,'last')-find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')+2):f*lokal(runde,j,20)];%,...
            diagonale = [diagonale,repmat(f*lokal(runde,j,20),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2))-length(diagonale))];
        else
            diagonale = repmat(lokal(runde,j,6),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        end
        diagonale = diagonale(1:length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        diff = pegel(becken_pos(j,1)).inflow{spalte}(:,2)-diagonale';
        diff(diff<0)=0;
        lokal(runde,j,34) = sum(diff)*3600;
        
        f = .2;
        if lokal(runde,j,6)<f*lokal(runde,j,20)
            clear diagonale
            diagonale = [repmat(lokal(runde,j,6),1,find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')-2),...
                lokal(runde,j,6):(f*lokal(runde,j,20)-lokal(runde,j,6))/(find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>f*lokal(runde,j,20),1,'last')-find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')+2):f*lokal(runde,j,20)];%,...
            diagonale = [diagonale,repmat(f*lokal(runde,j,20),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2))-length(diagonale))];
        else
            diagonale = repmat(lokal(runde,j,6),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        end
        diagonale = diagonale(1:length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        diff = pegel(becken_pos(j,1)).inflow{spalte}(:,2)-diagonale';
        diff(diff<0)=0;
        lokal(runde,j,35) = sum(diff)*3600;
        
        f = .1;
        if lokal(runde,j,6)<f*lokal(runde,j,20)
            clear diagonale
            diagonale = [repmat(lokal(runde,j,6),1,find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')-2),...
                lokal(runde,j,6):(f*lokal(runde,j,20)-lokal(runde,j,6))/(find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>f*lokal(runde,j,20),1,'last')-find(pegel(becken_pos(j,1)).inflow{spalte}(:,2)>lokal(runde,j,6),1,'first')+2):f*lokal(runde,j,20)];%,...
            diagonale = [diagonale,repmat(f*lokal(runde,j,20),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2))-length(diagonale))];
        else
            diagonale = repmat(lokal(runde,j,6),1,length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        end
        diagonale = diagonale(1:length(pegel(becken_pos(j,1)).inflow{spalte}(:,2)));
        diff = pegel(becken_pos(j,1)).inflow{spalte}(:,2)-diagonale';
        diff(diff<0)=0;
        lokal(runde,j,36) = sum(diff)*3600;
    else
        lokal(runde,j,28:36) = 0;
    end
    
    lokal(runde,j,37) = length(zeitraum_inflowAW);
    lokal(runde,j,38) = length(zeitraum_inflowA);
    lokal(runde,j,39) = length(zeitraum_outflowA);
    
    if lokal(runde,j,19) >= lokal(runde,j,6)
        [~,zeitraumA_B] = fun_zeitraum_B(pegel(becken_pos(j)).inflow_wasim(:,3),lokal(runde,j,6));
        lokal(runde,j,40) = (sum(pegel(becken_pos(j)).inflow_wasim(zeitraumA_B,3))-length(zeitraumA_B)*lokal(runde,j,6))*3600;
        lokal(runde,j,43) = length(zeitraumA_B);
    else
        lokal(runde,j,40) = 0;
        lokal(runde,j,43) = NaN;
    end
    
    if lokal(runde,j,20) >= lokal(runde,j,6)
        [~,zeitraumA_B] = fun_zeitraum_B(pegel(becken_pos(j)).inflow{spalte}(:,2),lokal(runde,j,6));
        lokal(runde,j,41) = (sum(pegel(becken_pos(j)).inflow{spalte}(zeitraumA_B,2))-length(zeitraumA_B)*lokal(runde,j,6))*3600;
        lokal(runde,j,44) = length(zeitraumA_B);
    else
        lokal(runde,j,41) = 0;
        lokal(runde,j,44) = NaN;
    end
    
    if lokal(runde,j,21) >= lokal(runde,j,6)
        [~,zeitraumA_B] = fun_zeitraum_B(pegel(becken_pos(j)).outflow{spalte}(:,2),lokal(runde,j,6));
        lokal(runde,j,42) = (sum(pegel(becken_pos(j)).outflow{spalte}(zeitraumA_B,2))-length(zeitraumA_B)*lokal(runde,j,6))*3600;
        lokal(runde,j,45) = length(zeitraumA_B);
    else
        lokal(runde,j,42) = 0;
        lokal(runde,j,45) = NaN;
    end
    
    % seit 03.05.2019:
    lokal(runde,j,46) = sum(pegel(becken_pos(j)).becken(spalte).S(:,2)>0);    % Einstaudauer
    vol_diff = pegel(becken_pos(j,1)).inflow{spalte}(:,2)-pegel(becken_pos(j,1)).outflow{spalte}(:,2);
    lokal(runde,j,47) = sum(vol_diff(vol_diff>0))*3600; % zurückgehaltenes Volumen
end
end

%% Datenzusammenstellung für alle Beobachtungspunkte
for k = 1:length(beobachtungspunkte)
    regional(runde,k,1) = beobachtungspunkte(k);
    becken_in_TG{runde,k} = [];
    for j = 1:size(becken_pos,1)
        if any(pegel(beobachtungspunkte(k)).gage.no==pegel(becken_pos(j,1)).neighbors.alle_nachfolger)==1
            becken_in_TG{runde,k} = [becken_in_TG{runde,k},becken_pos(j,2)];
        end
    end
    regional(runde,k,2) = length(becken_in_TG{runde,k});
    regional(runde,k,3) = pegel(beobachtungspunkte(k)).gage.area;
    regional(runde,k,4) = pegel(beobachtungspunkte(k)).gage.S_gesamt(end,1);
    regional(runde,k,5) = pegel(beobachtungspunkte(k)).gage.S_gesamt(end,2);
    regional(runde,k,6) = pegel(beobachtungspunkte(k)).gage.S_gesamt(end,3);
    regional(runde,k,7) = regional(runde,k,4)/regional(runde,k,3)*1e-3;
    regional(runde,k,8) = regional(runde,k,5)/regional(runde,k,3)*1e-3;
    regional(runde,k,9) = regional(runde,k,6)/regional(runde,k,3)*1e-3;
    
    if regional(runde,k,2)-1>0
        regional(runde,k,10) = sum([lokal(runde,find(ismember([lokal(runde,j,2)],[becken_in_TG{runde,k}])==1),18)],2)/(regional(runde,k,2)-1);
    else
        regional(runde,k,10) = 2;
    end
end
