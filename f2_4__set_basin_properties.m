function [pegel,startpegel] = f2_4__set_basin_properties(pegel,beckenpunkte_select2,becken_pos,Beckenkombination,E_no,E_opt,feste_Drosselweite,runde,existing_structures,Gebiet,mOpt)

bemessung = 3;  % HQ 1000

for j = 1:size(becken_pos,1)
    no_pegel = becken_pos(j,1);
    
    if mOpt~=13
        k = length(pegel(no_pegel).beckenparameter)+1;
        
        %%% Beckendefinition
        beckenpunkte_zeile = find([beckenpunkte_select2.no_becken]==becken_pos(j,2));
        pegel(no_pegel).beckenparameter(k).beckennummer = becken_pos(j,2);
        % Beckenvolumen in Betriebsfall (ohne HWE): < 50'000 m^3
        pegel(no_pegel).beckenparameter(k).S_max = beckenpunkte_select2(beckenpunkte_zeile).V_HRB(bemessung);
        
        pegel(no_pegel).beckenparameter(k).S0 = 0;
        % Einstauhöhe im Betriebsfall (ohne HWE): < 4 m
        pegel(no_pegel).beckenparameter(k).h_W = beckenpunkte_select2(beckenpunkte_zeile).h_HRB(bemessung);
        pegel(no_pegel).beckenparameter(k).h_HWE = beckenpunkte_select2(beckenpunkte_zeile).h_HWE;
        pegel(no_pegel).beckenparameter(k).h_D = beckenpunkte_select2(beckenpunkte_zeile).h_dam;
        pegel(no_pegel).beckenparameter(k).b_W = beckenpunkte_select2(beckenpunkte_zeile).b_wehr(bemessung);
        % Volumen
        pegel(no_pegel).beckenparameter(k).V_W = beckenpunkte_select2(beckenpunkte_zeile).V_HRB(bemessung);
        pegel(no_pegel).beckenparameter(k).V_HWE = beckenpunkte_select2(beckenpunkte_zeile).V_HWE;
        pegel(no_pegel).beckenparameter(k).V_D = beckenpunkte_select2(beckenpunkte_zeile).V_dam;
        % Öffnungsgröße:
        pegel(no_pegel).beckenparameter(k).D_min = beckenpunkte_select2(beckenpunkte_zeile).D_min;
        pegel(no_pegel).beckenparameter(k).D_max = beckenpunkte_select2(beckenpunkte_zeile).D_max;
        pegel(no_pegel).beckenparameter(k).ShAs = beckenpunkte_select2(beckenpunkte_zeile).ShAs;
        pegel(no_pegel).beckenparameter(k).geo = beckenpunkte_select2(beckenpunkte_zeile).elev_dam;
        pegel(no_pegel).beckenparameter(k).f = beckenpunkte_select2(beckenpunkte_zeile).h_dam - beckenpunkte_select2(beckenpunkte_zeile).h_HWE;
        
        pegel(no_pegel).beckenparametrisierung(k);
        
        %%% Drosselweite
        if feste_Drosselweite ==1
            % vorgegebene Drosselweite:
            zeile_becken = find([Beckenkombination.becken_real(:,runde)]==becken_pos(j,2));
            spalte = find([pegel(no_pegel).beckenparameter.beckennummer]==becken_pos(j,2));
            pegel(no_pegel).gage.AL_min(spalte) = Beckenkombination.drosselweite(zeile_becken,runde,E_opt);
            pegel(no_pegel).gage.AL_max(spalte) = Beckenkombination.drosselweite(zeile_becken,runde,E_opt);
        end
        
        % Beckenkombination analysieren und Pegel oberstrom
        % identifizieren (--> Identifikation möglicher Startpegel)
        pegel = becken_vorgaenger(pegel,no_pegel);
        
    else
        if isempty(pegel(no_pegel).neighbors.vorgaenger)==1
            pegel(no_pegel).becken_yn = 1;
            discharge_total = pegel(no_pegel).discharge_qdir(:,2)+pegel(no_pegel).discharge_qifl(:,2)+pegel(no_pegel).discharge_qbas(:,2);
            pegel(no_pegel).discharge_shares = [pegel(no_pegel).discharge_qdir(:,2)./discharge_total,...
                pegel(no_pegel).discharge_qifl(:,2)./discharge_total,pegel(no_pegel).discharge_qbas(:,2)./discharge_total];
        else
            pegel(no_pegel).becken_yn = 0;
        end
        pegel(no_pegel).beckenparameter.S_max = 0;
        pegel = becken_vorgaenger(pegel,no_pegel);
    end
    
end

%%% potential start gages for the routing
pegel_temp = findobj(pegel,'becken_yn',1,'-and','becken_vor_yn',0); % basin at the gage but not before
startpegel = pegel_temp(1).gage.no; % define one of the potential start gages as starting gage of the routing algorithm

%%% inclusion of existing structures (after start basin
% existing retention basins
if Gebiet == 2
    for p_include = 9:12
        pegel(p_include).beckenparameter.difference = [];
        pegel = fun_include_basin(p_include,existing_structures(E_no,p_include),pegel);
    end
end

